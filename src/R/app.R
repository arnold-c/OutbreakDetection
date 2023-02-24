library(shiny)
library(tidyverse)
library(rlang)
library(paletteer)


selection_tabs <- tabsetPanel(
  id = "selection_tabs",
  type = "hidden",
  tabPanel(
    "Total",
    # Input: Selector for tested cases ----
    sliderInput(
      inputId = "num_tested",
      label = "Choose the number of individuals tested:",
      value = 100,
      min = 0,
      max = 100
      ),
    # Input: Selector for Prior Probability ----
    sliderInput(
      inputId = "prior_prob",
      label = "Choose the prior probability (%) of a positive measles test:",
      value = 60,
      min = 0,
      max = 100
      ),
    # Input: Selector for GS ----
    sliderInput(
      inputId = "gs_spec",
      label = "Choose the specificity (%) of the Gold Standard test:",
      value = 97.3,
      min = 0,
      max = 100
      ),
    sliderInput(
      inputId = "gs_sens",
      label = "Choose the sensitivity (%) of the Gold Standard test:",
      value = 98.1,
      min = 0,
      max = 100
      ),
    # Input: Selector for RDT ----
    sliderInput(
      inputId = "rdt_spec",
      label = "Choose the specificity (%) of the Rapid Diagnostic test:",
      value = 96.2,
      min = 0,
      max = 100
      )
  ),
  tabPanel(
    "Positive", 
    # Input: Selector for positive cases ----
    sliderInput(
      inputId = "pos_tested",
      label = "Choose the number of positive cases identified by the Gold Standard:",
      value = 3,
      min = 0,
      max = 20
      ),
    # Input: Selector for Quantile ----
    sliderInput(
      inputId = "percentile_pos",
      label = "Choose the percentile of test distribution to return the desired number of true positive cases:",
      value = 50,
      min = 0,
      max = 100
      )
  )
)

# Define UI for dataset viewer app ----
ui <- fluidPage(

  # App title ----
  titlePanel("RDT - Gold Standard Equivalence"),

  # Sidebar layout with a input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Select whether to input number of gold standard positives or
      # total individuals tested
      radioButtons(
        inputId = "pos_or_total",
        label = "Do you want to input the number of gold standard positives or the total number of individuals tested?",
        selected = "Positive",
        choices = list(
          "Positive",
          "Total"
        )
      ),
      # Input: Selection tabs based on user input about postive vs total cases
      # being computed
      selection_tabs,


      # Input: Selector for RDT ----
      sliderInput(
        inputId = "rdt_sens",
        label = "Choose the sensitivity (%) of the Rapid Diagnostic test:",
        value = 90.0,
        min = 0,
        max = 100
        ),

      
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      
      p("
      Hiebert et al. 2021 evaluated the accuracy of 8 commerically available assays for measles-IgM detection (https://journals.asm.org/doi/10.1128/JCM.03161-20) and noted that the most accurate assay (Clin-Tech Microimmune ELISA) has a sensitivity of 98.1% and specificity of 97.3% (when treating equivocal as positive). This was one of three methods to exceed the Enzygnost EIA benchmark, so was chosen as the gold standard default values.
      "),
      br(),      
      p("
      Warrener et al. 2011 evaluated the performance of a lateral flow rapid point-of-care measles IgM test (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3165981/), which can be performed in < 30 minutes and with a serum or oral fluid sample, demonstarting a sensitivity of 90.0% and specificity of 96.2%.
      "),
      conditionalPanel(
        # Conditional: Show the plots and table for total case calculations
        condition = "input.pos_or_total == 'Total'",
        tabsetPanel(
          tabPanel(
            "Plots",
              # Output: Plot of results ----
              splitLayout(
                plotOutput("test_pos_plot"),
                plotOutput("ppv_plot")
              )
            ),
          tabPanel(
            # Output: Table of results ----
            "Table",
            tableOutput("table")
            )
        )
      ),
      conditionalPanel(
        # Conditional: Show the plots and table for positive case calculations
        condition = "input.pos_or_total == 'Positive'",
        tabPanel(
          "Negative Binomial",
          splitLayout(
            plotOutput("test_dist_dens_plot"),
            plotOutput("test_dist_cum_prob_plot")
            # tableOutput("test_dist_table")
          )
        )
      )
          
      
      
      )
    )
  )


# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {

  observeEvent(input$pos_or_total, {
    updateTabsetPanel(
      inputId = "selection_tabs", selected = input$pos_or_total
      )
  })

  # General values ----
  dis_pos <- reactive({
    input$num_tested * (input$prior_prob / 100)
  })

  dis_neg <- reactive({
    input$num_tested * (1 - (input$prior_prob / 100))
  })

  # General functions ----
  true_pos <- function(pos, test){
    sens <- paste0({{ test }}, "_sens")
    pos * (input[[sens]] / 100)
  }

  false_pos <- function(neg, test){
    spec <- paste0({{ test }}, "_spec")
    neg * (1 - (input[[spec]] / 100))
  }

  true_neg <- function(neg, test){
    spec <- paste0({{ test }}, "_spec")
    neg * (input[[spec]] / 100)
  }

  false_neg <- function(pos, test){
    sens <- paste0({{ test }}, "_sens")
    pos * (1 - (input[[sens]] / 100))
  }

  # Gold Standard Results ----
  gs_res  <- reactive({
    tibble(
      pos = true_pos(dis_pos(), test = "gs") + 
        false_pos(dis_neg(), test = "gs"),
      neg = true_neg(dis_neg(), test = "gs") + 
        false_neg(dis_pos(), test = "gs"),
      ppv = 100 * true_pos(dis_pos(), test = "gs") / 
        (true_pos(dis_pos(), test = "gs") + false_pos(dis_neg(), test = "gs")),
      npv = 100 * true_neg(dis_neg(), test = "gs") /
        (true_neg(dis_neg(), test = "gs") + false_neg(dis_pos(), test = "gs"))
    )
  })

  # Rapid Diagnostic Results ----
  rdt_res  <- reactive({
    tibble(
      pos = true_pos(dis_pos(), test = "rdt") + 
        false_pos(dis_neg(), test = "rdt"),
      neg = true_neg(dis_neg(), test = "rdt") + 
        false_neg(dis_pos(), test = "rdt"),
      ppv = 100 * true_pos(dis_pos(), test = "rdt") / 
        (true_pos(dis_pos(), test = "rdt") + false_pos(dis_neg(), test = "rdt")),
      npv = 100 * true_neg(dis_neg(), test = "rdt") /
        (true_neg(dis_neg(), test = "rdt") + false_neg(dis_pos(), test = "rdt"))

    )
  })

  # Generate a summary table of the test detection ----
  output$table  <- renderTable({
    tibble(
      "Disease positive" = dis_pos(),
      "Disease negative" = dis_neg(),
      "Gold Standard positive" = gs_res()$pos,
      "Gold Standard negative" = gs_res()$neg,
      "Gold Standard PPV (%)" = gs_res()$ppv,
      "Gold Standard NPV (%)" = gs_res()$npv,
      "Rapid Diagnostic positive" = rdt_res()$pos,
      "Rapid Diagnostic negative" = rdt_res()$neg,
      "Rapid Diagnostic PPV (%)" = rdt_res()$ppv,
      "Rapid Diagnostic NPV (%)" = rdt_res()$npv
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Result",
      values_to = "Value"
    )
  })

  # Generate general sens-spec plots ----
  sens_range <- reactive({seq(0, 1, 0.05)})
  spec_range <- reactive({seq(0, 1, 0.05)})
  gen_data <- reactive({
    tidyr::crossing(sens = sens_range(), spec = spec_range())
    })

  gen_plot_data <- reactive({
    purrr::map2_dfr(
      .x = gen_data()$sens,
      .y = gen_data()$spec,
      .f = function(.x, .y){
        true_pos <- dis_pos() * .x
        true_neg <- dis_neg() * .y
        false_pos <- dis_neg() * (1 - .y)
        false_neg <- dis_pos() * (1 - .x)

        test_pos <- true_pos + false_pos
        test_neg <- true_neg + false_neg

        ppv <- true_pos / test_pos
        npv <- true_neg / test_neg

        return(tibble(
          sens = .x,
          spec = .y,
          true_pos = true_pos,
          true_neg = true_neg,
          false_pos = false_pos,
          false_neg = false_neg,
          test_pos = test_pos,
          test_neg = test_neg,
          ppv = ppv,
          npv = npv
        ))
      }
    )
  })

  input_data <- reactive({
    tibble(
      gs_sens = input$gs_sens / 100,
      gs_spec = input$gs_spec / 100,
      rdt_sens = input$rdt_sens / 100,
      rdt_spec = input$rdt_spec / 100
    )
  })

  output$test_pos_plot <- renderPlot({

    gen_plot_data() %>%
      ggplot() +
      geom_contour_filled(aes(x = sens, y = spec, z = test_pos)) +
      geom_point(
        data = input_data(),
        aes(x = gs_sens, y = gs_spec), 
        shape = 21, size = 2, color = "black", fill = "gold"
        ) +
      geom_point(
        data = input_data(),
        aes(x = rdt_sens, y = rdt_spec),
        shape = 21, size = 2, color = "white", fill = "black"
        ) +
      scale_fill_paletteer_d("cartography::sand.pal", dynamic = TRUE) +
      labs(
        title = "Number of Positive Tests by Sensitivity and Specificity",
        subtitle = paste0(
          "Prior Probability = ", 
          input$prior_prob, 
          "%, Total Tests = ",
          input$num_tested,
          ", Disease Positives = ",
          dis_pos()
          ),
        x = "Sensitivity",
        y = "Specificity",
        fill = "Number of Positive Tests"
      ) +
      theme(legend.position = "bottom")
  })

  output$ppv_plot <- renderPlot({

    gen_plot_data() %>%
      ggplot() +
      geom_contour_filled(aes(x = sens, y = spec, z = ppv)) +
      geom_point(
        data = input_data(), 
        aes(x = gs_sens, y = gs_spec), 
        shape = 21, size = 2, color = "black", fill = "gold"
        ) +
      geom_point(
        data = input_data(),
        aes(x = rdt_sens, y = rdt_spec),
        shape = 21, size = 2, color = "white", fill = "black"
        ) +
      scale_fill_paletteer_d("cartography::harmo.pal", dynamic = TRUE) +
      labs(
        title = "PPV by Sensitivity and Specificity",
        subtitle = paste0(
          "Prior Probability = ",
          input$prior_prob,
          "%, Total Tests = ",
          input$num_tested,
          ", Disease Positives = ",
          dis_pos()
          ),
        x = "Sensitivity",
        y = "Specificity",
        fill = "PPV"
      ) +
      theme(legend.position = "bottom")
  })

  
  quantile_tested <- reactive({
    input$pos_tested + 
      qnbinom(p = input$percentile_pos/100, size = input$pos_tested, prob = input$rdt_sens/100)
    })

  failures_vec <- reactive({seq(0, quantile_tested() - input$pos_tested, 1)})

  # Determine how many failures should be included in the plots (dependent on
  # where the vertical line indicating the number of tests would lie on the axes)
  total_test_dist <- reactive({
    if (quantile_tested() > input$pos_tested + 10) {
      data <- tibble(failures = failures_vec())
    } else if(
        dnbinom(5, input$pos_tested, input$rdt_sens/100) < 0.05 & 
          quantile_tested() <= input$pos_tested + 5
        ) {
      data <- tibble(failures = 0:5)
    } else {
      data <- tibble(failures = 0:10)
    }
    
    data %>%
      mutate(
        total = failures + input$pos_tested,
        dens = dnbinom(failures, input$pos_tested, prob = input$rdt_sens/100),
        cum_prob = pnbinom(failures, input$pos_tested, prob = input$rdt_sens/100)
      ) 
  })

  output$test_dist_dens_plot <- renderPlot({
    total_test_dist() %>%
      ggplot() +
      geom_col(aes(x = total, y = dens, fill = dens)) +
      geom_vline(
        xintercept = quantile_tested(), 
        color = "grey80",
        size = 1.5,
        linetype = "dashed"
        ) +
      scale_x_continuous(
        breaks = seq(
          min(total_test_dist()$total), 
          max(total_test_dist()$total), 
          1
          )
      ) +
      scale_fill_gradient(low = "#fbe3c2", high = "#92351e") +
      labs(
        title = paste0(
          "PMF of Tests Required to Return ",
          input$pos_tested,
          " Positive Tests"
          ),
        subtitle = paste0(
          "Test Sensitivity = ",
          input$rdt_sens,
          "%, ",
          input$percentile_pos,
          "-th Percentile"
          ),
        x = "Total Tests",
        y = "Probability Density",
        caption = paste0(
          "The ",
          input$percentile_pos,
          "-th percentile of tests required to return at least ",
          input$pos_tested,
          " positive tests is shown as a vertical line"
          )
      ) +
      theme(legend.position = "none")
  })

  output$test_dist_cum_prob_plot <- renderPlot({
    total_test_dist() %>%
      ggplot() +
      geom_col(aes(x = total, y = cum_prob, fill = cum_prob)) +
      geom_hline(
        yintercept = input$percentile_pos/100,
        color = "grey80",
        size = 1.5,
        linetype = "dashed"
        ) +
      geom_vline(
        xintercept = quantile_tested(),
        color = "grey80",
        size = 1.5,
        linetype = "dashed"
        ) +
      scale_x_continuous(
        breaks = seq(
          min(total_test_dist()$total), 
          max(total_test_dist()$total), 
          1
          )
      ) +
      scale_fill_gradient(low = "#abc9c8", high = "#0a3351") +
      labs(
        title = paste0(
          "CDF of Tests Required to Return ",
          input$pos_tested,
          " Positive Tests"
          ),
        subtitle = paste0(
          "Test Sensitivity = ",
          input$rdt_sens,
          "%, ",
          input$percentile_pos,
          "-th Percentile"
          ),
        x = "Total Tests",
        y = "Cumulative Probability",
        caption = paste0(
          "The ",
          input$percentile_pos,
          "-th percentile of tests required to return at least ",
          input$pos_tested,
          " positive tests is shown as a vertical line,\nwith the horizontal line showing the associated cumulative probability (",
          input$percentile_pos,
          "%)"
          )
      ) +
      theme(legend.position = "none")
  })

}

shinyApp(ui = ui, server = server)