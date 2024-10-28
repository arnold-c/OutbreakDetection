// Some definitions presupposed by pandoc's typst output.
#let blockquote(body) = [
  #set text( size: 0.92em )
  #block(inset: (left: 1.5em, top: 0.2em, bottom: 0.2em))[#body]
]

#let horizontalrule = [
  #line(start: (25%,0%), end: (75%,0%))
]

#let endnote(num, contents) = [
  #stack(dir: ltr, spacing: 3pt, super[#num], contents)
]

#show terms: it => {
  it.children
    .map(child => [
      #strong[#child.term]
      #block(inset: (left: 1.5em, top: -0.4em))[#child.description]
      ])
    .join()
}

// Some quarto-specific definitions.

#show raw.where(block: true): set block(
    fill: luma(230),
    width: 100%,
    inset: 8pt,
    radius: 2pt
  )

#let block_with_new_content(old_block, new_content) = {
  let d = (:)
  let fields = old_block.fields()
  fields.remove("body")
  if fields.at("below", default: none) != none {
    // TODO: this is a hack because below is a "synthesized element"
    // according to the experts in the typst discord...
    fields.below = fields.below.amount
  }
  return block.with(..fields)(new_content)
}

#let unescape-eval(str) = {
  return eval(str.replace("\\", ""))
}

#let empty(v) = {
  if type(v) == "string" {
    // two dollar signs here because we're technically inside
    // a Pandoc template :grimace:
    v.matches(regex("^\\s*$")).at(0, default: none) != none
  } else if type(v) == "content" {
    if v.at("text", default: none) != none {
      return empty(v.text)
    }
    for child in v.at("children", default: ()) {
      if not empty(child) {
        return false
      }
    }
    return true
  }

}

// Subfloats
// This is a technique that we adapted from https://github.com/tingerrr/subpar/
#let quartosubfloatcounter = counter("quartosubfloatcounter")

#let quarto_super(
  kind: str,
  caption: none,
  label: none,
  supplement: str,
  position: none,
  subrefnumbering: "1a",
  subcapnumbering: "(a)",
  body,
) = {
  context {
    let figcounter = counter(figure.where(kind: kind))
    let n-super = figcounter.get().first() + 1
    set figure.caption(position: position)
    [#figure(
      kind: kind,
      supplement: supplement,
      caption: caption,
      {
        show figure.where(kind: kind): set figure(numbering: _ => numbering(subrefnumbering, n-super, quartosubfloatcounter.get().first() + 1))
        show figure.where(kind: kind): set figure.caption(position: position)

        show figure: it => {
          let num = numbering(subcapnumbering, n-super, quartosubfloatcounter.get().first() + 1)
          show figure.caption: it => {
            num.slice(2) // I don't understand why the numbering contains output that it really shouldn't, but this fixes it shrug?
            [ ]
            it.body
          }

          quartosubfloatcounter.step()
          it
          counter(figure.where(kind: it.kind)).update(n => n - 1)
        }

        quartosubfloatcounter.update(0)
        body
      }
    )#label]
  }
}

// callout rendering
// this is a figure show rule because callouts are crossreferenceable
#show figure: it => {
  if type(it.kind) != "string" {
    return it
  }
  let kind_match = it.kind.matches(regex("^quarto-callout-(.*)")).at(0, default: none)
  if kind_match == none {
    return it
  }
  let kind = kind_match.captures.at(0, default: "other")
  kind = upper(kind.first()) + kind.slice(1)
  // now we pull apart the callout and reassemble it with the crossref name and counter

  // when we cleanup pandoc's emitted code to avoid spaces this will have to change
  let old_callout = it.body.children.at(1).body.children.at(1)
  let old_title_block = old_callout.body.children.at(0)
  let old_title = old_title_block.body.body.children.at(2)

  // TODO use custom separator if available
  let new_title = if empty(old_title) {
    [#kind #it.counter.display()]
  } else {
    [#kind #it.counter.display(): #old_title]
  }

  let new_title_block = block_with_new_content(
    old_title_block, 
    block_with_new_content(
      old_title_block.body, 
      old_title_block.body.body.children.at(0) +
      old_title_block.body.body.children.at(1) +
      new_title))

  block_with_new_content(old_callout,
    block(below: 0pt, new_title_block) +
    old_callout.body.children.at(1))
}

// 2023-10-09: #fa-icon("fa-info") is not working, so we'll eval "#fa-info()" instead
#let callout(body: [], title: "Callout", background_color: rgb("#dddddd"), icon: none, icon_color: black) = {
  block(
    breakable: false, 
    fill: background_color, 
    stroke: (paint: icon_color, thickness: 0.5pt, cap: "round"), 
    width: 100%, 
    radius: 2pt,
    block(
      inset: 1pt,
      width: 100%, 
      below: 0pt, 
      block(
        fill: background_color, 
        width: 100%, 
        inset: 8pt)[#text(icon_color, weight: 900)[#icon] #title]) +
      if(body != []){
        block(
          inset: 1pt, 
          width: 100%, 
          block(fill: white, width: 100%, inset: 8pt, body))
      }
    )
}

#let article(
  // Article's Title
  title: "Article Title",
  header-title: none,

  // A dictionary of authors.
  // Dictionary keys are authors' names.
  // Dictionary values are meta data of every author, including
  // label(s) of affiliation(s), email, contact address,
  // or a self-defined name (to avoid name conflicts).
  // Once the email or address exists, the author(s) will be labelled
  // as the corresponding author(s), and their address will show in footnotes.
  //
  // Example:
  // (
  //   "Author Name": (
  //     "affiliation": "affil-1",
  //     "email": "author.name@example.com", // Optional
  //     "address": "Mail address",  // Optional
  //     "name": "Alias Name" // Optional
  //   )
  // )
  authors: (),

  // A dictionary of affiliation.
  // Dictionary keys are affiliations' labels.
  // These labels show be constent with those used in authors' meta data.
  // Dictionary values are addresses of every affiliation.
  //
  // Example:
  // (
  //   "affil-1": "Institution Name, University Name, Road, Post Code, Country"
  // )
  affiliations: (),

  // The paper's abstract.
  abstract: [],

  // The paper's keywords.
  keywords: (),

  // The path to a bibliography file if you want to cite some external
  // works.
  bib: none,
  bib-title: "References",

  // Word count
  word-count: false,

  // Line numbers
  line-numbers: false,

  // Paper's content
  body
) = {
  // Set document properties
  set bibliography(title: bib-title)

  // Line numbers have not yet been implemented in a release version, but are coming soon
  // https://github.com/typst/typst/issues/352
  // https://github.com/typst/typst/pull/4516
  //if line-numbers {
  //  set par.line(numbering: "1")
  //  show figure: set par.line(numbering: none)
  //}

  set document(title: title, author: authors.keys())
  set page(numbering: "1", number-align: center)
  set text(font: ("Linux Libertine", "STIX Two Text", "serif"), lang: "en")
  show footnote.entry: it => [
    #set par(hanging-indent: 0.7em)
    #set align(left)
    #numbering(it.note.numbering, ..counter(footnote).at(it.note.location())) #it.note.body
  ]
  show figure.caption: emph
  set table(
    fill: (x, y) => {
      if y == 0 {gray}
    }
  )
  show table.cell: it => {
    if it.y == 0 {
      strong(it)
    } else {
      it
    }
  }

  // Title block
  align(center)[
    #block(text(size: 14pt, weight: "bold", title))

    #v(1em)

    // Authors and affiliations

    // Restore affiliations' keys for looking up later
    // to show superscript labels of affiliations for each author.
    #let inst_keys = affiliations.keys()

    // Authors' block
    #block([
      // Process the text for each author one by one
      #for (ai, au) in authors.keys().enumerate() {
        let au_meta = authors.at(au)
        // Don't put comma before the first author
        if ai != 0 {
          text([, ])
        }
        // Write auther's name
        let au_name = if au_meta.keys().contains("name") {au_meta.name} else {au}
        text([#au_name])

        // Get labels of author's affiliation
        let au_inst_id = au_meta.affiliation
        let au_inst_primary = ""
        // Test whether the author belongs to multiple affiliations
        if type(au_inst_id) == "array" {
          // If the author belongs to multiple affiliations,
          // record the first affiliation as the primary affiliation,
          au_inst_primary = affiliations.at(au_inst_id.first())
          // and convert each affiliation's label to index
          let au_inst_index = au_inst_id.map(id => inst_keys.position(key => key == id) + 1)
          // Output affiliation
          super([#(au_inst_index.map(id => [#id]).join([,]))])
        } else if (type(au_inst_id) == "string") {
          // If the author belongs to only one affiliation,
          // set this as the primary affiliation
          au_inst_primary = affiliations.at(au_inst_id)
          // convert the affiliation's label to index
          let au_inst_index = inst_keys.position(key => key == au_inst_id) + 1
          // Output affiliation
          super([#au_inst_index])
        }

        // Corresponding author
        if au_meta.keys().contains("corresponding") and (au_meta.corresponding == "true" or au_meta.corresponding == true) {
          [#super[,]#footnote(numbering: "*")[
            Corresponding author. #au_name. Address:
            #if not au_meta.keys().contains("address") or au_meta.address == "" {
              [#au_inst_primary.]
            }
            #if au_meta.keys().contains("email") {
              [Email: #link("mailto:" + au_meta.email.replace("\\", "")).]
            }
          ]]
        }

        if au_meta.keys().contains("equal-contributor") and (au_meta.equal-contributor == "true" or au_meta.equal-contributor == true){
          if ai == 0 {
            [#super[,]#footnote(numbering: "*")[Equal contributors.]<ec_footnote>]
          } else {
            [#super[,]#footnote(numbering: "*", <ec_footnote>)]
          }
        }
      }
    ])

    #v(1em)

    // Affiliation block
    #align(left)[#block([
      #set par(leading: 0.4em)
      #for (ik, key) in inst_keys.enumerate() {
        text(size: 0.8em, [#super([#(ik+1)]) #(affiliations.at(key))])
        linebreak()
      }
    ])]
  ]

  if header-title == "true" {
    header-title = title
  }

  set page(
    header: [
      #set text(8pt)
      #align(right)[#header-title]
    ],
  )
  if word-count {
    import "@preview/wordometer:0.1.2": word-count, total-words
    show: word-count.with(exclude: (heading, table, figure.caption))
  }

  // Abstract and keyword block
  if abstract != [] {
    pagebreak()

    block([
        #if word-count {
        import "@preview/wordometer:0.1.2": word-count, word-count-of, total-words
        text(weight: "bold", [Word count: ])
        text([#word-count-of(exclude: (heading))[#abstract].words])
      }
      #heading([Abstract])
      #abstract

      #if keywords.len() > 0 {
        linebreak()
        text(weight: "bold", [Key words: ])
        text([#keywords.join([; ]).])
      }
    ])

    v(1em)
  }

  // Display contents

  pagebreak()

  show heading.where(level: 1): it => block(above: 1.5em, below: 1.5em)[
    #set pad(bottom: 2em, top: 1em)
    #it.body
  ]

  set par(first-line-indent: 0em)

  if word-count {
      import "@preview/wordometer:0.1.2": word-count, word-count-of, total-words
      text(weight: "bold", [Word count: ])
      text([#word-count-of(exclude: (heading, table, figure.caption))[#body].words])
  }

  body

}


 // title: "Title",
 // header-title: "Header Title",
 // authors: (
 //   "Author 1": (
 //     "affiliation": ("affil-1", "affil-2"),
 //     "email": "a1-email",
 //   ),
 //   "Author 2": (
 //     "affiliation": ("affil-1", "affil-3"),
 //     "email": "a2-email",
 //   ),
 // ),
 // affiliations: (
 //   "affil-1": "Affiliation 1",
 //   "affil-2": "Affiliation 2",
 //   "affil-3": "Affiliation 3",
 // ),
 // abstract: [
 //   == Background
 //   Abstract Background
 //
 //   == Methods
 //   Abstract Methods
 //
 //   == Result
 //   Abstract Results
 //
 //   == Conclusions
 //   Abstract Conclusions
 // ],
 // keywords: (),
 // bib: "refs.bib",
 // bib-title: "Refs",
 // word-count: false,
 // line-numbers: false,

#show: body => article(
      title: "The Need to Develop a Holistic Infectious Disease Surveillance System",
        header-title: "true",
        authors: (
                    "Callum R.K. Arnold": (
            affiliation: ("PSU-Bio", "CIDD"),
            corresponding: "true",
            
            email: "cfa5228\@psu.edu",
      ),                    "Alex C. Kong": (
            affiliation: ("Hopkins-IH"),
            
            
            
      ),                    "Amy K. Winter": (
            affiliation: ("UGA"),
            
            
            
      ),                    "William J. Moss": (
            affiliation: ("Hopkins-IH", "Hopkins-Epi"),
            
            
            
      ),                    "Bryan N. Patenaude": (
            affiliation: ("Hopkins-IH"),
            
            
            
      ),                    "Matthew J. Ferrari": (
            affiliation: ("PSU-Bio", "CIDD"),
            
            
            
      ),        ),
        affiliations: (
          "PSU-Bio": "Department of Biology, Pennsylvania State University, University Park, PA, USA 16802",
          "CIDD": "Center for Infectious Disease Dynamics, Pennsylvania State University, University Park, PA, USA 16802",
          "Hopkins-IH": "Department of International Health, Johns Hopkins Bloomberg School of Public Health, Baltimore, MD, USA 21205",
          "UGA": "Department of Epidemiology, College of Public Health, University of Georgia, Athens, GA, USA 30602",
          "Hopkins-Epi": "Department of Epidemiology, Johns Hopkins Bloomberg School of Public Health, Baltimore, MD, USA 21205"    ),
          bib: "OD.bib",
            keywords: ("Rapid-Diagnostic Tests","ELISA","Infectious Disease Surveillance","Outbreak Detection"),
        abstract: [== Background
<background>
Infectious disease surveillance and outbreak detection systems often utilize diagnostic testing to validate case identification. The metrics of sensitivity, specificity, and positive predictive value are commonly discussed when evaluating the performance of diagnostic tests, and to a lesser degree, the performance of outbreak detection systems. However, the interaction of the two levels’ (the test and the alert system) metrics, is typically overlooked. Here, we describe how equivalent regions of detection accuracy can exist over a range of diagnostic test characteristics, examining the sensitivity to background noise structure and magnitude.

== Methods
<methods>
We generated a stochastic SEIR model with importation to simulate true measles and non-measles sources of febrile rash (noise) daily incidence. We generated time series of febrile rash (i.e., measles clinical case definition) by summing the daily incidence of measles and either independent Poisson noise or non-measles dynamical noise (consistent with rubella virus). For each time series we assumed a fraction of all cases were seen at a healthcare clinic, and a subset of those were diagnostically confirmed using a test with sensitivity and specificity consistent with either a rapid diagnostic test (RDT) or an enzyme-linked immunosorbent assay (ELISA). From the resulting time series of test-positive cases, we define an outbreak alert as the exceedance of a threshold by the 7-day rolling average of observed (test positive) cases. For each threshold level, we calculated percentages of alerts that were aligned with an outbreak (analogous to the positive predictive value), the percentage of outbreaks detected (analogous to the sensitivity), and combined these two measures into an accuracy metric for outbreak detection. We selected the optimal threshold as the value that maximizes accuracy. We show how the optimal threshold and resulting accuracy depend on the diagnostic test, reporting rate, and the type and magnitude of the non-measles noise.

== Results
<results>
The optimal threshold for each test increased monotonically as the percentage of clinic visits who were tested increased. With Poisson-only noise, similar outbreak detection accuracies could be achieved with imperfect RDT-like tests as with ELISA-like diagnostic tests (c. 93%), given moderately high testing rates. With larger delays (14 days) between the ELISA test administration and result date, RDTs could outperform the ELISA. Similar numbers of unavoidable cases and outbreak alert delays could be achieved between the test types. With dynamical noise, however, the accuracy of ELISA scenarios was far superior to those achieved with RDTs (c.~93% vs.~73%). For dynamical noise, RDT-based scenarios typically favored more sensitive alert threshold than ELISA-based scenarios (at a given reporting rate), observed with lower numbers of unavoidable cases and detection delays.

== Conclusions
<conclusions>
The performance of an outbreak detection system is highly sensitive to the structure and the magnitude of background noise. Under the assumption that the noise is relatively static over time, RDTs can perform as well as ELISA in a surveillance system. However, when the noise is temporally correlated, as from a separate SEIR process, imperfect tests cannot overcome their accuracy limitations through higher testing rates.

],
        word-count: true,
        line-numbers: true,
    body,
)

= Background
<background>
At the heart of an outbreak detection system is a surveillance program built upon case detection, often utilizing individual diagnostic tests as necessary components of epidemiological investigations before the declaration of an outbreak @murrayInfectiousDiseaseSurveillance2017@zhouEliminationTropicalDisease2013@pahoIntegratedApproachCommunicable2000@worldhealthorganizationMeaslesOutbreakGuide2022@craggOutbreakResponse2018. For diseases with non-specific symptoms, accurate measurement tools are often necessary to confidently and correctly ascribe changes in symptom prevalence within a population to a particular disease, and therefore detect outbreaks of specific pathogens. As a result, it has been commonplace for surveillance systems to be developed around high-accuracy tests, such as Polymerase Chain Reaction (PCR) tests and immunoglobulin (Ig) tests, when financially and logistically feasible @gastanaduyMeasles2019@commissionerCoronavirusCOVID19Update2020@grasslyComparisonMolecularTesting2020@ezhilanSARSCoVMERSCoVSARSCoV22021@worldhealthorganizationCholera2023@essentialprogrammeonimmunizationepiimmunizationClinicalSpecimensLaboratory2018. Depending on the disease in question, either sensitivity (the ability to correctly detect a true positive individual) or specificity (the ability to correctly discount a true negative) will be prioritized, as they often are at odds with each other @westreichDiagnosticTestingScreening2019@shrefflerDiagnosticTestingAccuracy2024@parikhUnderstandingUsingSensitivity2008. This balance is commonly defined within the Target Product Profile (TPP) of a test @worldhealthorganizationTargetProductProfiles, which is a set of minimum characteristics that should be met for production and widespread use, helping to guide research and development. For example, in the wake of the 2013 Ebola outbreak in Guinea a TPP was developed that listed the minimum acceptable sensitivity of 95% and specificity of 99% @chuaCaseImprovedDiagnostic2015. Recognizing that Ebola is not the major cause of fever and other non-specific symptoms in the region, it is arguably more important to prioritize the specificity of the disease, although the authors note that the severity of the infection requires a high level of sensitivity as the consequences of a missed case are dire at an individual and population level @chuaCaseImprovedDiagnostic2015.

Much like the accuracy of an individual test, outbreak detection systems face the same issue regarding the prioritization of sensitive or specific alerts @germanSensitivityPredictiveValue2000@worldhealthorganizationOperationalThresholds2014@lewisTimelyDetectionMeningococcal2001. For many disease systems, particularly in resource constrained environments where the burden of infectious diseases is highest @gbd2019childandadolescentcommunicablediseasecollaboratorsUnfinishedAgendaCommunicable2023@roserBurdenDisease2023, cases are counted and if a pre-determined threshold is breached, be that weekly, monthly, or some combination of the two, an alert is triggered that may launch a further investigation and/or a response @worldhealthorganizationMeaslesOutbreakGuide2022@worldhealthorganizationOperationalThresholds2014. In effect, this discretizes a distinctly continuous phenomenon (observed cases) into a binary measure, outbreak or no outbreak, for decision making purposes. Reactive interventions, such as vaccination campaigns and non-pharmaceutical based interventions, designed to reduce transmission or limit and suppress outbreaks, early action has the potential to avert the most cases @atkinsAnticipatingFutureLearning2020@taoLogisticalConstraintsLead@graisTimeEssenceExploring2008@ferrariTimeStillEssence2014@worldhealthorganizationConfirmingInvestigatingManaging2009@minettiLessonsChallengesMeasles2013. While this framing would point towards a sensitive (i.e., early alert) surveillance system being optimal, each action comes with both direct and indirect financial and opportunity costs stemming from unnecessary activities that limit resources for future response capabilities. Just as the balance of sensitivity and specificity of a test for an individual must be carefully evaluated, so must the balance at the outbreak level.

The concept of using incidence-based alert triggers to define the discrete event of an "outbreak" with characteristics analogous to individual tests has been well documented in the case of meningitis, measles, and malaria @worldhealthorganizationMeaslesOutbreakGuide2022@lewisTimelyDetectionMeningococcal2001@worldhealthorganizationConfirmingInvestigatingManaging2009@trotterResponseThresholdsEpidemic2015@cooperReactiveVaccinationControl2019@zalwangoEvaluationMalariaOutbreak2024@kanindaEffectivenessIncidenceThresholds2000. However, an overlooked, yet critical, aspect of an outbreak detection system is the interplay between the individual test and outbreak alert characteristics. With their success within malaria surveillance systems, and particularly since the COVID-19 pandemic, rapid diagnostic tests (RDTs) have garnered wider acceptance, and their potential for use in other disease systems has been gaining interest @warrenerEvaluationRapidDiagnostic2023. Despite concerns of their lower diagnostic accuracy slowing their adoption until recently @millerAddressingBarriersDevelopment2015, RDTs generally have reduced cold-chain requirements @brownRapidDiagnosticTests2020 and faster speed of result that has been show to outweigh the cost of false positive/negative results in some settings @warrenerEvaluationRapidDiagnostic2023@mcmorrowMalariaRapidDiagnostic2011@larremoreTestSensitivitySecondary2021@middletonModelingTransmissionMitigation2023.

In this paper we examine how the use of imperfect diagnostic tests affects the performance of outbreak detection for measles outbreaks in the context of a febrile rash surveillance system that includes both measles and non-measles cases. Because measles symptoms are non-specific, it is important to account for non-measles sources of febrile rash e.g., rubella, parvovirus, varicella, etc., producing the potential for false positive results in the context of imperfect tests. Currently, measles outbreaks are declared on the basis of either suspected measles cases (i.e., an individual with fever and maculopapular rash @essentialprogrammeonimmunizationepiimmunizationvaccinesandbiologicalsivbMeasles2018) alone, cases confirmed by enzyme-linked immunosorbent assay (ELISA) to detect the presence of measles-specific IgM antibodies, or a combination of the two, depending on the elimination status of the region, with countries nearer elimination increasing the use of ELISA diagnostic confirmation for suspected cases @essentialprogrammeonimmunizationepiimmunizationvaccinesandbiologicalsivbMeasles2018. Each of these detection systems have their flaws. Although clinical case definition is very fast and requires limited resources, it is highly sensitive, and in the face of high "background noise" from non-measles sources of febrile rash, can lead to low positive predictive value (PPV), i.e., the probability that an alert accurately captures an outbreak @hutchinsEvaluationMeaslesClinical2004. And while ELISA confirmation is the standard diagnostic test for measles surveillance and has higher specificity @essentialprogrammeonimmunizationepiimmunizationvaccinesandbiologicalsivbMeasles2018, the training and facility requirements generally mean that samples must be transported from the point of care to a separate laboratory, which incurs both costs and delays #emph[#strong[\[ADD REFERENCE HERE\]];];. In resource-poor settings these delays may be days to weeks. A rapid diagnostic test offers the opportunity to provide a compromise between diagnostic accuracy and timeliness, and recent developments show encouraging signs in the field as well as in theory @warrenerEvaluationRapidDiagnostic2023@brownRapidDiagnosticTests2020.

By examining the combination of alert threshold and individual test characteristic in a modeling study that explicitly incorporates dynamical background noise, we aim to illustrate the need to develop a TPP for the whole detection system, not just one component. To evaluate the alert system performance, we develop a set of outbreak definition criteria and surveillance metrics, drawing inspiration from acceptance sampling, ecological surveillance systems, and epidemiological surveillance systems guidelines and reviews @yemshanovAcceptanceSamplingCosteffective2020@christensenHerdlevelInterpretationTest2000@muratoEvaluationSamplingMethods2020@sternAutomatedOutbreakDetection1999@calbaSurveillanceSystemsEvaluation2015@guidelinesworkinggroupUpdatedGuidelinesEvaluating2001. Using these metrics we overcome issues encountered by early warning systems that rely on dynamical values such as $R_(upright("effective"))$ in defining outbreaks @jombartRealtimeMonitoringCOVID192021@proverbioPerformanceEarlyWarning2022@stolermanUsingDigitalTraces2023@brettAnticipatingEpidemicTransitions2018@salmonMonitoringCountTime2016, for example, characterizing the end of an epidemic period is important in a time series where multiple outbreaks will occur.

= Methods
<methods>
== Model Structure
<model-structure>
We constructed a stochastic compartmental non-age structured SEIR model of measles and simulated using a modified Tau-leaping algorithm with a time step of 1 day, utilizing binomial draws to ensure compartment sizes remained positive valued @chatterjeeBinomialDistributionBased2005@gillespieApproximateAcceleratedStochastic2001. To stochastically reintroduce infections and cause recurring epidemics, commuter-style imports were added that are proportional to the size of the population and $R_0$ @keelingModelingInfectiousDiseases2008, and the transmission parameter ($beta_t$) is sinusoidal with a period of one year and 20% seasonal amplitude. $R_0$ was set to 16, with a latent period of 10 days and infectious period of 8 days. The population was initialized with 500,000 individuals with Ghana-like birth and vaccination rates, and the final results were scaled up to the approximate 2022 population size of Ghana (33 million) @worldbankGhana. The full table of parameters can be found in @tbl-model-parameters. All simulations and analysis was completed in Julia version 1.10.5 @bezansonJuliaFreshApproach2017, with all code stored at #link("https://github.com/arnold-c/OutbreakDetection");.

#figure([
#let parameters = csv("parameters.csv")
#let import_rate = $(1.06*μ*R_0)/(√(N))$
#let parameter_labels = ( "Parameters", $R_0$, $"Latent period ("#sym.sigma")"$, $"Infectious period ("#sym.gamma")"$, "Seasonal amplitude", $"Birth rate ("#sym.mu")"$, $"Vaccination rate at birth ("#sym.rho")"$)

#table(
  columns: 3,
  ..for ((label), (.., measles, nonmeasles)) in parameter_labels.zip(parameters) {
      (label, measles, nonmeasles)
  },
  [Importation rate], table.cell(colspan: 2, align: center, $(1.06*μ*R_0)/(√(N))$),
  [Population size (N)], table.cell(colspan: 2, align: center, "500,000, scaled to 33M"),
  [Initial proportion susceptible], table.cell(colspan: 2, align: center, "0.05"),
  [Initial proportion exposed], table.cell(colspan: 2, align: center, "0.0"),
  [Initial proportion infected], table.cell(colspan: 2, align: center, "0.0"),
  [Initial proportion recovered], table.cell(colspan: 2, align: center, "0.95"),
)
], caption: figure.caption(
position: bottom, 
[
Compartmental model parameters
]), 
kind: "quarto-float-tbl", 
supplement: "Table", 
)
<tbl-model-parameters>


To examine the sensitivity of the detection system to background noise, we layered the measles incidence time series with one of four noise time series structures: Poisson-only noise; or dynamical noise with rubella-like parameters that could be in- or out-of-phase, or with independent seasonality to the measles dynamics. Only dynamical in-phase noise and Poisson-only noise are presented in this manuscript. For each noise structure, five different magnitudes of noise were simulated and evaluated: the total annual incidence of noise in each system was either the same magnitude as the average annual measles incidence, two times, four times, six times, or eight times the average annual measles incidence. 100 time series of 100 years were simulated for each scenario, before summarizing the distributions of outbreak detection methods.

#emph[Mention that the same model was used to define the noise as the outbreaks.]

== Defining Outbreaks
<defining-outbreaks>
It is common to use expert review to define outbreaks when examining empirical data, but this is not feasible in a modeling study where tens of thousands of years are being simulated. To account for this, many studies only simulate a single outbreak within a time series (repeating this short stochastic simulation multiple times to ensemble results), define an outbreak as a period where $R_(upright("effective"))$ \> 1, or use a threshold of \> 2 standard deviations (s.d.) over the mean seasonal incidence observed in empirical data (or from a 'burn-in' period of the simulation) @sternAutomatedOutbreakDetection1999@jombartRealtimeMonitoringCOVID192021@stolermanUsingDigitalTraces2023@salmonMonitoringCountTime2016@teklehaimanotAlertThresholdAlgorithms2004@leclereAutomatedDetectionHospital2017. Each method has its uses, but to evaluate the performance of an outbreak detection system in an endemic region where multiple sequential epidemics are expected it is important to clearly define the bounds of the outbreak, which can only be achieved by 2 s.d. \> mean ($R_(upright("effective"))$ will be less than 1 after an outbreak’s peak, but still within what can be reasonably defined as the outbreak’s bounds). This, however, assumes strong seasonal forcing and regular periodicity of incidence to produce a smooth enough baseline, which is not present as countries near measles elimination status @grahamMeaslesCanonicalPath2019. Here we define a true measles outbreak as a region of the time series that meets the following three criteria:

- The daily measles incidence must be greater than, or equal to, 5 cases
- The daily measles incidence must remain above 5 cases for greater than, or equal to, 30 consecutive days
- The total measles incidence must be great than, or equal to, 500 cases within the bounds of the outbreak

Each of these is a necessary, but not sufficient, condition to the definition of an outbreak; all must be met. For example, a region may produce 1000 total cases within the lower and upper bounds defined by the daily incidence being $gt.eq$ 5 cases, but if this period only lasts for 20 days it would not be considered an outbreak as it would not be possible to mount a response to. The incidence of non-measles febrile rash (i.e., noise) does not affect the outbreak status of a region but may affect the alert status triggered by the testing protocol.

== Triggering Alerts
<triggering-alerts>
#emph[Need to mention that an alert is any time the threshold is exceeded, regardless of the number of consecutive days of exceedance.] #emph[Clearly state that 1 day of exceedance is and 30 days of consecutive exceedance are both treated as a single alert i.e., the same.]

The use of alert thresholds is common in locations burdened by measles, where access to high performance computing is limited making nowcasting-style approaches to outbreak detection impractical. For this reason, we define an "alert" as any single day where the 7-day moving average of the daily incidence is greater than, or equal to, a pre-specified alert threshold. To examine the interaction of test and surveillance systems we vary the alert threshold, between 1 and 15 cases per day, before computing the evaluation metrics and identifying the optimal threshold for each combination of individual test, testing rate, noise structure, and noise magnitude. We also examined an alert method that required either the daily incidence or the 7-day moving average of daily incidence to exceed the alert threshold (see supplemental results).

Each day, 60% of the measles and non-measles febrile rash cases visit the clinic for treatment, and a percentage of these clinic visits are tested, as all clinic visits are deemed to be suspected measles cases because they meet the clinical case definition. This percentage of clinic visits that are tested is varied between 10% and 60%, in 10% increments, for all combinations of diagnostic test (except clinical case definition) and alert threshold, defining the "testing scenario". Each testing scenario uses one of the following tests:

- An RDT equivalent with 85% sensitivity and specificity, and 0-day lag in result return. That is, 85% of true measles cases will be correctly labelled as positive, and 15% of non-measles febrile rash individuals that are tested will be incorrectly labelled as positive for measles. This acts as a lower bound of acceptability for a new measles RDT.
- An RDT equivalent with 90% sensitivity and specificity, and 0-day lag in result return.
- An ELISA equivalent test with 100% sensitivity and specificity with a 0-day test result delay.
- An ELISA equivalent test with 100% sensitivity and specificity with a 14-day test result delay.

== Optimal Thresholds
<optimal-thresholds>
For each of the previously defined tests, we calculate the number of test positive cases (that will include false positive and negative cases resulting from imperfect diagnostic tests), which is used to categorize the time series by alert status in conjunction with the alert threshold. The overlap with true outbreak status, or lack thereof, is used to classify the following metrics, and has been illustrated in @fig-outbreak-schematic:

- The percentage of alerts that are 'correct'. This is analogous to the PPV of the alert system. Note that it is possible to trigger multiple alerts within a single outbreak, and each would be considered correct.
- The percentage of outbreaks that are detected. This is analogous to the sensitivity of the alert system.
- The detection delay. If an alert precedes the start of the outbreak, as long as it successfully "captures" the outbreak within its bounds, it is considered to be a 'correct' alert, resulting in a negative detection delay i.e., an early warning triggered by false positives.
- The number of unavoidable and avoidable outbreak cases and deaths. Avoidable cases are defined as those that occur within an outbreak after a correct alert is first triggered i.e., cases that could theoretically be prevented with a perfectly effective and timely response. Unavoidable cases are the inverse: those that occur before a correct alert, or those that occur in an undetected outbreak. Using a fitted age distribution of measles cases in Ghana @mintaProgressMeaslesElimination2023 and the age-specific case fatality ratio for Ghana in 2022 @sbarraEstimatingNationallevelMeasles2023, we scaled the values for each one-year age cohort, per annum. In practice, not all cases deemed avoidable are (due to imperfect and delays in responses), but to minimize the sensitivity of the results to the response implementation and operational constraints we are counting them as such.

#figure([
#box(image("manuscript_files/figure-typst/fig-outbreak-schematic-output-1.png"))
], caption: figure.caption(
position: bottom, 
[
A schematic of the outbreak definition and alert detection system. A) Noise time series. B) Measles incidence time series. C) Observed time series resulting from testing noise & measles cases that visit the healthcare facility. The orange bands/vertical lines represent regions of the measles time series that meet the outbreak definition criteria. The green bands/vertical lines represent regions of the observed (measles - noise) time series that breach the alert threshold (the horizontal dashed line), and constitute an alert.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-outbreak-schematic>


For each combination of diagnostic test and testing rate, the optimal alert threshold is calculated by selecting the threshold that produces the highest accuracy. Each of the 100 simulations produces an accuracy (the mean of the percentages of alerts that are correct and the percentage of the outbreaks that are detected), and the median accuracy of this distribution is used to determine the optimal alert threshold. Each testing scenario (the combination of the diagnostic test, testing rate, and alert threshold), is then compared (using the metrics above) at the optimal threshold that is specific for the scenario. This allows for conclusions to be made about the surveillance system as a whole, instead of just single components.

= Results
<results>
The threshold that maximized outbreak detection accuracy given a testing scenario depends heavily on the noise structure and magnitude. For example, when the ratio of average noise incidence was 8 times higher than the average measles incidence, the optimal threshold for an RDT with 90% sensitivity/specificity and 40% testing of clinic visits was 4 cases per day with the dynamical in-phase noise, whereas it is 5 cases per day when Poisson-only noise was simulated (@tbl-optimal-thresholds). This corresponds to an accuracy of 72% and 93%, respectively (@fig-accuracy). Much higher accuracy can be achieved with imperfect tests under Poisson noise because large spikes of non-measles febrile rash that would lead to many false positives and trigger erroneous alerts do not occur. With dynamical noise, this possibility can and does occur and, when the ratio of non-measles to measles cases is sufficiently high, drastically impacts the accuracy of an alert system. The results from ELISA tests are identical between the different noise structures: with a perfect test there are no false positive results. For dynamical noise, RDTs never perform as well as an ELISA, even under scenarios that include a 14-day lag between test procurement and result receipt. However, with Poisson noise, once testing reaches 20% of clinic visits, equivalence can be observed at approximately 91% accuracy (@fig-accuracy).

#figure([
#let optimal_thresholds = csv("optimal-thresholds.csv")

#table(
  columns: 9,
  fill: (x, y) => {
    if y == 0 {gray}
    if y == 1 {gray}
  },
  align: center,
  [], table.cell(colspan: 2, align: center, "Test Characteristic"), table.cell(colspan: 6, align: center, "Testing Rate"),
  ..optimal_thresholds.flatten()
)
], caption: figure.caption(
position: bottom, 
[
Optimal threshold for each testing scenario, when the average noise incidence is 8 times higher than the average measles incidence. A) the noise structure is dynamical, and the seasonality is in-phase with the measles incidence. B) the noise structure is Poisson only
]), 
kind: "quarto-float-tbl", 
supplement: "Table", 
)
<tbl-optimal-thresholds>


Moving across a single row of @tbl-optimal-thresholds, as the testing rate increases, so does the optimal threshold, because more individuals are testing positive. However, the changes in accuracy are not monotonically increasing, as one might anticipate. Even with a perfect test, there are sequential cells that increase the testing rate, yet accuracy decreases (@fig-accuracy). For example, the ELISA equivalent with a 0-day test result lag outbreak detection accuracy decreases from 91% to 87% when clinic testing rates increase from 30% to 40% of clinic visits. The reason behind this unintuitive result stems from the use of integer valued thresholds. While computing the moving average test positive incidence helps to smooth out the alert system, an integer valued threshold can result in step-changes of accuracy between two threshold values, and the expected increase in the alert system’s PPV from a higher threshold value is outweighed by the loss in the alert system’s sensitivity to detecting outbreaks. Even with a perfect test, the alert system needs to discriminate between endemic/imported cases and epidemic cases; increasing the testing rate will result in higher numbers of test positive individuals, and a lower threshold can result in an overly sensitive alert system, triggering for measles infections, but not those within an outbreak. A higher threshold will face the opposite issue; not triggering for smaller outbreaks, which is more likely to be an issue at lower testing rates. Examining the number of unavoidable cases (@fig-unavoidable) and the detection delays (@fig-delay), we can see that the decrease in accuracy results from a decrease in the sensitivity of the alert system: the detection delay increases by 6 days, and the unavoidable cases by c.~2300. Although the testing rate increases, so does the optimal threshold, from 3 test positives to 4 test positives per day. The subsequent change in testing rate (40 - 50%) is associated with no change in the optimal threshold for the ELISA, and as expected, the number of unavoidable cases and detection delay decrease, resulting in a more sensitive alert system, with an higher system accuracy indicating that the increase more than offset the decrease in the system’s PPV.

#figure([
#box(image("manuscript_files/figure-typst/fig-accuracy-output-1.png"))
], caption: figure.caption(
position: bottom, 
[
The accuracy of outbreak detection systems under different testing rates and noise structures. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-accuracy>


#figure([
#box(image("manuscript_files/figure-typst/fig-unavoidable-output-1.png"))
], caption: figure.caption(
position: bottom, 
[
The number of unavoidable cases of outbreak detection systems under different testing rates and noise structures. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-unavoidable>


#figure([
#box(image("manuscript_files/figure-typst/fig-delay-output-1.png"))
], caption: figure.caption(
position: bottom, 
[
The detection delay of outbreak detection systems under different testing rates and noise structures. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-delay>


For Poisson noise, similar detection delays are observed for all test and noise magnitudes, with variation by testing rate (mean of -3.7 to 36.1 days). Under dynamical noise, there are clearer differences in the performance of ELISA and RDTs, with the separation of outcomes occurring later than observed for detection accuracy (8 time noise magnitude vs.~2 times, respectively) (@fig-accuracy). With large amounts of dynamical noise (8 times the measles incidence), the mean detection delay of the 90% and 85% RDTs range from -17.5 days to 3.2 days, and from -25.2 days to -3.4 days, respectively. Negative delays indicate that alerts are being triggered before the start of the outbreak and is correlated with the proportion of the time series that is under alert, with larger negative delays associated with more and/or longer alert periods (@fig-alert-proportion, Supplemental Figure 2).

#figure([
#box(image("manuscript_files/figure-typst/fig-alert-proportion-output-1.png"))
], caption: figure.caption(
position: bottom, 
[
The difference between proportion of the time series in alert for outbreak detection systems under different testing rates and noise structures. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-alert-proportion>


= Discussion
<discussion>
The performance of an outbreak detection system is highly sensitive to the structure and level of background noise in the simulation. Despite the mean daily noise incidence set to equivalent values between the dynamical and Poisson-only simulations, drastically different results are observed.

Under the assumption that non-measles febrile rash is relatively static in time (Poisson), RDTs can perform as well, if not better than ELISA tests at moderate to high testing rates, and at a fraction of the cost. However, if it is expected that the noise is non-stationary, imperfect tests cannot overcome their accuracy limitations through higher testing rates, saturating at c.~74% accuracy, relative to ELISA’s 93%. This occurs because, despite the same average incidence of noise in each (comparable) scenario, the relative proportion of measles to noise varies differently throughout the time series, exacerbating the effects of imperfect diagnostic tests that will produce higher rates of false positives and negatives than ELISA-like diagnostics. Imperfect tests are also more susceptible to incorrect results due to the varying prevalence throughout the time series, so the PPV of the system will not be static throughout the time series.

- #emph[Make point about balance of false positive/negative outbreak detection]
- #emph[Analysis isn’t true optimization:]
  - #emph[Requires explicit decisions about preference for speed / false alerts vs higher PPV]
  - #emph[Surveillance is counting for action (WHO quote)]

The purpose of routine surveillance is to characterize the infection landscape. With strong public health infrastructure and infectious disease surveillance programs, it is possible to develop a strong understanding of the shape of febrile rash cases, regardless of source. With this information, countries can tailor their future activities to rely more or less heavily upon RDTs, depending on the dynamics of the target disease and its relationship to background noise, favoring RDTs when there are low levels of noise and ELISAs during large rubella outbreaks, for example.

== Limitations and Strengths
<limitations-and-strengths>
To our knowledge, this is one of the first simulation studies to examine the relationship between individual test characteristics and the wider surveillance program. By explicitly modeling the interaction between the two, we make a case that surveillance systems should take a holistic approach; ignoring one component can lead to drastically different, and suboptimal, results. Additionally, by defining outbreak bounds concretely we have been able to calculate metrics of outbreak detection performance that draw parallels to those used when evaluating individual diagnostic tests, allowing for intuitive and simple implementation of this method in resource-constrained environments, something that is not possible with most outbreak detection and early warning system simulations in the literature. An evaluation of all outbreak detection algorithms is beyond the scope of this work, but a more computationally expensive approach based on nowcasting incidence may help overcome the shortcomings of RDTs in high-noise scenarios.

For computational simplicity, this paper did not include demography in the model structure. And while a simulation-based approach allows for complete determination of true infection status i.e., measles vs non-measles febrile rash cases, and therefore an accurate accounting of the outbreak and alert bounds, these simulations do not specifically represent any real-world setting. The evaluation of empirical data does provide this opportunity, but at the cost of knowing the true infection status of individuals and confounding of multiple variables, limiting analysis to only those who are observed (i.e., not those in the community who do not visit a healthcare center), and removing the possibility to explore the sensitivity of the results to parameters of interest to a surveillance program e.g., testing rate, and the test itself.

Additionally, is has been well documented that the performance of an individual test is highly sensitive to its timing within a person’s infection cycle @gastanaduyMeasles2019@larremoreTestSensitivitySecondary2021@middletonModelingTransmissionMitigation2023@kisslerViralDynamicsAcute2021@ratnamPerformanceIndirectImmunoglobulin2000, so it is possible that different conclusions would be drawn if temporal information about the test acquisition was included in the simulation.

Finally, the optimal threshold for a testing scenario depends heavily on the costs ascribed to incorrect actions, be that failing to detect an outbreak or incorrectly mounting a response for an outbreak that doesn’t exist. In the simulations we have weighted them equally, but it is likely that they should not be deemed equivalent: missing an outbreak leads to many thousands of cases and associated DALYs, whereas an unnecessary alert would generally launch an initial low-cost investigation for full determination of the outbreak status. This is particularly important in countries with vast heterogeneity in transmission: different weightings should be applied to higher vs.~lower priority/risk regions to account for discrepancies in consequences of incorrect decisions. For example, a high priority zone could still benefit from a false alert as many high-risk healthcare regions within a country are targeted for Supplemental Immunization Activities, so a false alert would just hasten the (proactive) vaccination response.

Given these limitations, the explicit values (i.e., optimal thresholds, accuracies etc.) should be interpreted with caution, and the exact results observed in the real-world will likely be highly dependent on unseen factors, such as the proportion of measles and non-measles sources of febrile rash that seek healthcare. However, the general patterns should hold, and more importantly, the analysis framework provides a consistent and holistic approach to evaluating the trade-off between individual level tests and the alert system enacted to detect outbreaks. 

#pagebreak()
= Funding
<funding>
- #emph[Something about GAVI/Gates]

This work was supported by funding from the Office of the Provost and the Clinical and Translational Science Institute, Huck Life Sciences Institute, and Social Science Research Institutes at the Pennsylvania State University. The project described was supported by the National Center for Advancing Translational Sciences, National Institutes of Health, through Grant UL1 TR002014. The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH. The funding sources had no role in the collection, analysis, interpretation, or writing of the report.

= Acknowledgements
<acknowledgements>
== Author Contributions
<author-contributions>
#emph[Conceptualization:] CA, MJF

#emph[Data curation:] MJF, CA

#emph[Formal analysis:] CA, MJF

#emph[Funding acquisition:] MJF, WM, AW

#emph[Investigation:] CA, MJF

#emph[Methodology:] CA, MJF

#emph[Project administration:] MJF

#emph[Software:] CA

#emph[Supervision:] MJF, WM, AW, BP

#emph[Validation:] CA, MJF

#emph[Visualization:] CA

#emph[Writing - original draft:] CA, MJF

#emph[Writing - review and editing:] all authors.

== Conflicts of Interest and Financial Disclosures
<conflicts-of-interest-and-financial-disclosures>
The authors declare no conflicts of interest.

== Data Access, Responsibility, and Analysis
<data-access-responsibility-and-analysis>
Callum Arnold and Dr.~Matthew J. Ferrari had full access to all the data in the study and take responsibility for the integrity of the data and the accuracy of the data analysis. Callum Arnold and Dr.~Matthew J. Ferrari (Department of Biology, Pennsylvania State University) conducted the data analysis.

== Data Availability
<data-availability>
All code and data for the simulations can be found at #link("https://github.com/arnold-c/OutbreakDetection");.

#pagebreak()
#block[
] <refs>


 

#set bibliography(style: "springer-vancouver")


#bibliography("OD.bib")

