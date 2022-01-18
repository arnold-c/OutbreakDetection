library(shiny)
library(tidyverse)
library(rlang)
library(MetBrewer)


# Define UI for dataset viewer app ----
ui <- fluidPage(

  # App title ----
  titlePanel("RDT - Gold Standard Equivalence"),

  # Sidebar layout with a input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Selector for tested cases ----
      sliderInput(
        inputId = "num_tested",
        label = "Choose the number of individuals tested:",
        value = 100,
        min = 0,
        max = 100
        ),

      # Input: Selector for Prevalence ----
      sliderInput(
        inputId = "prev",
        label = "Choose the underlying prevalence (%) of measles:",
        value = 5,
        min = 0,
        max = 100
        ),

      # Input: Selector for GS ----
      sliderInput(
        inputId = "gs_sens",
        label = "Choose the sensitivity (%) of the Gold Standard test:",
        value = 98.1,
        min = 0,
        max = 100
        ),
      sliderInput(
        inputId = "gs_spec",
        label = "Choose the specificity (%) of the Gold Standard test:",
        value = 97.3,
        min = 0,
        max = 100
        ),

      # Input: Selector for RDT ----
      sliderInput(
        inputId = "rdt_sens",
        label = "Choose the sensitivity (%) of the Rapid Diagnostic test:",
        value = 90.0,
        min = 0,
        max = 100
        ),
      sliderInput(
        inputId = "rdt_spec",
        label = "Choose the specificity (%) of the Rapid Diagnostic test:",
        value = 96.2,
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
      tabsetPanel(
        # Output: Plot of results ----
        tabPanel(
          "Test Positive Plot",
          splitLayout(
            plotOutput("test_pos_plot"),
            plotOutput("ppv_plot")
          )
        ),

        # Output: Table of results ----
        tabPanel(
          "Table",
          tableOutput("table")
        )
      ),

    )
  )
)


# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {

  # General values ----
  dis_pos <- reactive({
    input$num_tested * (input$prev / 100)
  })

  dis_neg <- reactive({
    input$num_tested * (1 - (input$prev / 100))
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


  output$test_pos_plot <- renderPlot({

    input_data <- tibble(
      gs_sens = input$gs_sens / 100,
      gs_spec = input$gs_spec / 100,
      rdt_sens = input$rdt_sens / 100,
      rdt_spec = input$rdt_spec / 100
    )

    gen_plot_data() %>%
      ggplot() +
      geom_contour_filled(aes(x = sens, y = spec, z = test_pos)) +
      geom_point(
        data = input_data, 
        aes(x = gs_sens, y = gs_spec), 
        shape = 21, size = 2, color = "black", fill = "gold"
        ) +
      geom_point(
        data = input_data,
        aes(x = rdt_sens, y = rdt_spec),
        shape = 21, size = 2, color = "white", fill = "black"
        ) +
      scale_fill_manual(values = met.brewer("OKeeffe2", 10)) +
      labs(
        title = "Number of Positive Tests by Sensitivity and Specificity",
        subtitle = glue("Prevalence = {prev}, Total Tests = {tot}, Disease Positives = {dis_pos()}"),
        x = "Sensitivity",
        y = "Specificity",
        fill = "Number of Positive Tests"
      ) +
      theme(legend.position = "bottom")
  })

  output$ppv_plot <- renderPlot({
    
    input_data <- tibble(
      gs_sens = input$gs_sens / 100,
      gs_spec = input$gs_spec / 100,
      rdt_sens = input$rdt_sens / 100,
      rdt_spec = input$rdt_spec / 100
    )

    gen_plot_data() %>%
      ggplot() +
      geom_contour_filled(aes(x = sens, y = spec, z = ppv)) +
      geom_point(
        data = input_data, 
        aes(x = gs_sens, y = gs_spec), 
        shape = 21, size = 2, color = "black", fill = "gold"
        ) +
      geom_point(
        data = input_data,
        aes(x = rdt_sens, y = rdt_spec),
        shape = 21, size = 2, color = "white", fill = "black"
        ) +
      scale_fill_manual(values = met.brewer("Hokusai2", 10)) +
      labs(
        title = "PPV by Sensitivity and Specificity",
        subtitle = glue("Prevalence = {prev}, Total Tests = {tot}, Disease Positives = {dis_pos()}"),
        x = "Sensitivity",
        y = "Specificity",
        fill = "PPV"
      ) +
      theme(legend.position = "bottom")
  })

}

shinyApp(ui = ui, server = server)