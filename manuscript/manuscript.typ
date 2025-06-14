#import "template.typ": article
#import "metadata.typ": meta

#show: article.with(
  ..meta,
  abstract: [
    == Background
    Outbreak detection frequently relies on imperfect individual-level case diagnosis.
    Both outbreaks and cases are discrete events that can be misclassified and uncertainty at the case level may impact the performance of outbreak alert and detection systems.
    Here, we describe how the performance of outbreak detection depends on individual-level diagnostic test characteristics and population-level epidemiology and describe settings where imperfect individual-level tests can achieve consistent performance comparable to “perfect” diagnostic tests.

    == Methodology
    We generated a stochastic SEIR model to simulate daily incidence of measles (i.e., true) and non-measles (i.e., noise) febrile rash illness.
    We modeled non-measles sources as either independent static (Poisson) noise, or dynamical noise consistent with an independent SEIR process (e.g., rubella).
    Defining outbreak alerts as the exceedance of a threshold by the 7-day rolling average of observed test positives, we optimized the threshold that maximized outbreak detection accuracy across set of noise structures and magnitudes, diagnostic test accuracy (consistent with either a perfect test, or proposed rapid diagnostic tests), and testing rates.

    == Conclusions

    The optimal threshold for each diagnostic test typically increased monotonically with testing rate.
    With static noise, outbreak detection with RDT-like and perfect tests achieved accuracies of 90%, with comparable delays to outbreak detection.
    With dynamical noise, the accuracy of perfect test scenarios was superior to those achieved with RDTs ($approx$ 90% vs. $approx$ 80%).
    Outbreak detection accuracy declined as dynamical noise increased and leads to permanent alert status with RDT-like tests at very high noise.
    The performance of an outbreak detection system is highly sensitive to the structure and the magnitude of the background noise.
    Depending on the epidemiological context, outbreak detection using RDTs can perform as well as perfect tests.
  ],
  author-summary: [
    To respond to outbreaks of infectious diseases, we first need to detect them.
    This detection is inherently flawed, in part, due to imperfect diagnostic tests used to indicate whether individuals are positive or negative for a disease.
    We evaluated the impact of imperfect diagnostic tests for infectious disease on the accuracy and timeliness of outbreak detection in the context of a set of background infections that could be mistaken for the disease of interest, and consequently cause false positive test results.
    We find that when outbreak detection performance is highly dependent on the structure and magnitude of the background "noise" infections.
    When the rate of background infections far exceeds that of the target infection, and are dynamical, such that there are large peaks and troughs of "noise infections", imperfect diagnostic tests are not able to accurately distinguish the "signal" (target infections) from the background "noise".
    If the background "noise" infections are either less cyclical in their dynamics, or do not outnumber true infections by a great deal, imperfect diagnostic tests can perform well.
  ],
  // word-count: true,
  line-numbers: true,
  line-spacing: 2,
)

= Background

Diagnostics are medical devices and techniques that can be used to detect the presence of a specific pathogen in a host @DiagnosticsGlobal.
This may include #emph[in vivo] measures, such as x-ray imagery, or #emph[in vitro] tests to quantify the presence of the pathogen itself, e.g., polymerase chain reaction (PCR), or the host's response to the pathogen e.g., enzyme immunoassay/enzyme-linked immunosorbent assay (EIA/ELISA) @DiagnosticsGlobal @yangPCRbasedDiagnosticsInfectious2004 @alhajjEnzymeLinkedImmunosorbent2024.
Any given diagnostic will vary in its ability to correctly identify the presence of the pathogen, which is described by its sensitivity and specificity.
The sensitivity of a diagnostic is the ability to correctly identify a positive result, conditional on a positive individual being tested i.e., a true positive result @westreichDiagnosticTestingScreening2019 @shrefflerDiagnosticTestingAccuracy2024 @parikhUnderstandingUsingSensitivity2008.
The specificity is the opposite: the ability to correctly determine a true negative result, conditional on a negative individual being tested @westreichDiagnosticTestingScreening2019 @shrefflerDiagnosticTestingAccuracy2024 @parikhUnderstandingUsingSensitivity2008.
Due to the translation of quantitative measures e.g., immunoglobulin M (IgM) antibody titers, into a binary outcome (positive/negative), the sensitivity and specificity of a diagnostic are often at odds with one another.
For example, using a low optical density value to define the threshold for detection for an ELISA will produce a diagnostic that is highly sensitive, as it only requires a small host response to the pathogen and many resulting antibody titers will exceed this value.
However, this may lead to low specificity due to an increase in spurious false positive results from non-infected individuals.
To account for these differences, the target product profile (TPP) of a diagnostic provides a minimum set of characteristics that should be met, helping to guide the development and use @worldhealthorganizationTargetProductProfiles.

The choice to prioritize sensitivity or specificity will be pathogen and context specific.
When the cost of a false negative result is disproportionately high relative to a false positive, such as for Ebola @chuaCaseImprovedDiagnostic2015, highly specific tests may be preferred.
This balance will, however, vary as the prevalence of infection in a population varies.
Higher presence of infection in a population will increase the positive predictive value (PPV) of the test i.e., the probability that a positive test reflects a positive individual, that unlike the sensitivity of the test, is not conditioned upon the infection status of the tested individual @westreichDiagnosticTestingScreening2019 @shrefflerDiagnosticTestingAccuracy2024.
Regions of high disease burden may therefore prioritize test sensitivity, in contrast to a lower burden location's preference for test specificity and PPV, all else being equal.

At the heart of an outbreak detection system is a surveillance program that enumerates the baseline rate of case incidence and defines an outbreak as a time period with anomalously high incidence relative to that baseline @murrayInfectiousDiseaseSurveillance2017 @zhouEliminationTropicalDisease2013 @pahoIntegratedApproachCommunicable2000 @craggOutbreakResponse2018.
As many disease symptoms reflect generic host responses to infection e.g., febrile rash, and infection with a given pathogen can give rise to a wide range of disease symptoms and severity across individuals, accurate methods of case identification are required.
Given the imperfect nature of diagnostic classification, any result for an individual is uncertain.
Accumulating multiple individual test results to produce population-level counts will propagate this uncertainty, and may result in over- or under-counts due to a preponderance of the diagnostic to produce false positive and negative individual test results, respectively.
When the prevalence of a surveillance program's target disease is low relative to the prevalence of other sources of clinically-compatible cases (as might be expected at the start of an outbreak), the PPV of an individual diagnostic will decrease, increasing the number of false positives, making it harder to identify true anomalies in disease incidence.
As a result, it has been commonplace for infectious disease surveillance systems to be developed around high-accuracy tests, such as PCR and ELISA tests, when financially and logistically feasible @gastanaduyMeasles2019 @commissionerCoronavirusCOVID19Update2020@grasslyComparisonMolecularTesting2020@ezhilanSARSCoVMERSCoVSARSCoV22021 @worldhealthorganizationCholera2023 @essentialprogrammeonimmunizationepiimmunizationClinicalSpecimensLaboratory2018.

Outbreak detection systems, like diagnostic tests, must prioritize the sensitivity or specificity of an alert to detect an outcome (the outbreak) @germanSensitivityPredictiveValue2000@worldhealthorganizationOperationalThresholds2014 @lewisTimelyDetectionMeningococcal2001.
For many disease systems, particularly in resource constrained environments where the burden of infectious diseases is typically highest @gbd2019childandadolescentcommunicablediseasecollaboratorsUnfinishedAgendaCommunicable2023 @roserBurdenDisease2023, cases are counted and if a pre-determined threshold is breached #sym.dash.em be that weekly, monthly #sym.dash.em or some combination of the two, an alert is triggered that may launch a further investigation and/or a response @worldhealthorganizationMeaslesOutbreakGuide2022 @worldhealthorganizationOperationalThresholds2014.
In effect, this discretizes a distinctly continuous phenomenon (observed cases) into a binary measure (outbreak or no outbreak) for decision making purposes.
For reactive management approaches, such as vaccination campaigns and non-pharmaceutical based interventions that are designed to reduce transmission or limit and suppress outbreaks, early action has the potential to avert the most cases @atkinsAnticipatingFutureLearning2020@taoLogisticalConstraintsLead @graisTimeEssenceExploring2008 @ferrariTimeStillEssence2014 @worldhealthorganizationConfirmingInvestigatingManaging2009 @minettiLessonsChallengesMeasles2013.
While this framing would point towards a sensitive (i.e., early alert) surveillance system being optimal, each action comes with both direct and indirect financial and opportunity costs stemming from unnecessary activities that limit resource availability for future responses.
Much like the need to carefully evaluate the balance of an individual diagnostic test's sensitivity and specificity, it is essential to consider these characteristics at the outbreak level.

The concept of using incidence-based alert triggers to detect the discrete event of an outbreak with characteristics analogous to individual tests has been well documented in the case of meningitis, measles, and malaria @worldhealthorganizationMeaslesOutbreakGuide2022 @lewisTimelyDetectionMeningococcal2001 @worldhealthorganizationConfirmingInvestigatingManaging2009 @trotterResponseThresholdsEpidemic2015 @cooperReactiveVaccinationControl2019 @zalwangoEvaluationMalariaOutbreak2024 @kanindaEffectivenessIncidenceThresholds2000.
However, an overlooked, yet critical, aspect of an outbreak detection system is the interplay between the individual test and outbreak alert characteristics.
With their success within malaria surveillance systems, and particularly since the COVID-19 pandemic, rapid diagnostic tests (RDTs) have garnered wider acceptance, and their potential for use in other disease systems has been gaining interest @warrenerEvaluationRapidDiagnostic2023.
Despite concerns about their lower diagnostic accuracy slowing their adoption until recently @millerAddressingBarriersDevelopment2015, the reduced cold-chain requirements @brownRapidDiagnosticTests2020, reduced training and laboratory requirements and costs @essentialprogrammeonimmunizationepiimmunizationClinicalSpecimensLaboratory2018 @worldhealthorganizationMeaslesOutbreakGuide2022 @brownRapidDiagnosticTests2020, and faster speed of result provided by RDTs has been show to outweigh the cost of false positive/negative results in some settings @warrenerEvaluationRapidDiagnostic2023 @mcmorrowMalariaRapidDiagnostic2011 @larremoreTestSensitivitySecondary2021 @middletonModelingTransmissionMitigation2024.


In this paper, we examine how the use of imperfect diagnostic tests affects the performance of outbreak detection in the context of measles where RDTs are being developed with promising results @20240613_tpp_measles_rubell_FV_EN @brownRapidDiagnosticTests2020 @warrenerEvaluationRapidDiagnostic2023 @shonhaiInvestigationMeaslesOutbreak2015 (though not exclusively @seninMeaslesIgMRapid2024).
We evaluate the scenarios under which equivalence in outbreak detection can be achieved, where altering testing rates can offset the reduction in diagnostic discrimination of imperfect tests relative to perfect tests, and meaningful improvements can be attained with respect to specific metrics e.g., speed of response.
By examining the combination of the alert threshold and individual test characteristics in a modeling study that explicitly incorporates dynamical background noise, we illustrate the need to develop TPPs for surveillance programs as a whole.


= Methods
== Model Structure
We constructed a stochastic compartmental non-age structured Susceptible-Exposed-Infected-Recovered (SEIR) model of measles, and simulated disease transmission using a modified Tau-leaping algorithm with a time step of 1 day @gillespieApproximateAcceleratedStochastic2001.
We utilized binomial draws to ensure compartment sizes remained positive valued @chatterjeeBinomialDistributionBased2005.
We assumed that the transmission rate ($beta_t$) is sinusoidal with a period of one year and 20% seasonal amplitude.
$R_0$ was set to 16, with a latent period of 10 days and infectious period of 8 days @guerraBasicReproductionNumber2017 @gastanaduyMeasles2019.
The population was initialized with 500,000 individuals with Ghana-like birth and vaccination rates, and the final results were scaled up to the approximate 2022 population size of Ghana (33 million) @worldbankGhana.
Ghana was chosen to reflect a setting with a high-performing measles vaccination program that has not yet achieved elimination status (c. 80% coverage for two doses of measles-containing vaccine), and must remain vigilant to outbreaks @WHOImmunizationData @masreshaTrackingMeaslesRubella2024.
We assumed commuter-style imports at each time step to avoid extinction; the number of imports each day were drawn from a Poisson distribution with mean proportional to the size of the population and $R_0$ @keelingModelingInfectiousDiseases2008.
The full table of parameters can be found in @tbl_od-model-parameters.
All simulations and analyses were completed in Julia version 1.11.5 @bezansonJuliaFreshApproach2017, with all code stored at #link("https://github.com/arnold-c/OutbreakDetection").

#let table_math(inset: 6pt, size: 14pt, content) = table.cell(
  inset: inset,
  text(size: size, content),
)

#let import_rate = $(1.06*μ*R_0)/(√(N))$

#figure(
  table(
    columns: 3, align: horizon,
    [Parameters], [Measles], [Dynamical noise],
    [R0], [16], [5],
    [Latent period (s)], [10 days], [7 days],
    [Infectious period (g)], [8 days], [14 days],
    [Seasonal amplitude], [0.2], [0.2],
    [Vaccination rate at birth (r)], [80%], [(5-85)%],
    [Birth/death rate (m)], table.cell(
      colspan: 2,
      align: center,
      "27 per 1000 per annum",
    ),
    [Importation rate], table.cell(
      colspan: 2,
      align: center,
      table_math[$(1.06*μ*R_0)/(√(N))$],
    ),
    [Population size (N)], table.cell(
      colspan: 2,
      align: center,
      "500,000, scaled to 33M",
    ),
    [Initial proportion susceptible], table.cell(
      colspan: 2,
      align: center,
      "0.05",
    ),
    [Initial proportion exposed], table.cell(colspan: 2, align: center, "0.0"),
    [Initial proportion infected], table.cell(colspan: 2, align: center, "0.0"),
    [Initial proportion recovered], table.cell(
      colspan: 2,
      align: center,
      "0.95",
    ),
  ),
  caption: [Compartmental model parameters],
)
<tbl_od-model-parameters>

To examine the sensitivity of the detection system to background noise, we generated a time series of symptomatic febrile rash by combining the measles incidence time series with a noise time series.
The noise time series was modeled as either Poisson-only noise (subsequently referred to as static noise), to represent the incidence of non-specific febrile rash due to any of a number of possible etiologies, or dynamical noise modeled as a rubella SEIR process.
For static noise, the time series of non-measles febrile rash cases each day was constructed by independent draws from a Poisson distribution.
For dynamical noise, we generated time series of cases from an SEIR model that matched the measles model in structure, but had $R_0 = 5$, mean latent period of 7 days, and mean infectious period of 14 days.
We also added additional static noise drawn from a Poisson distribution with mean equal to 15% of the average daily rubella incidence to account for non-rubella sources of febrile rash (@tbl_od-model-parameters) @papadopoulosEstimatesBasicReproduction2022 @RubellaCDCYellow.
The seasonality for the rubella noise was simulated to be in-phase with measles.

For each noise structure, we simulated five magnitudes of noise ($Lambda$), representing the average daily noise incidence.
$Lambda$ was calculated as a multiple ($c$) of the average daily measles incidence ($angle.l Delta I_M angle.r$): $Lambda = c dot.op angle.l Delta I_M angle.r upright("where") c in { 1, 2, 4, 6, 8 }$.
Noise magnitudes will be denoted as $Lambda (c)$ for the rest of the manuscript e.g., $Lambda (8)$ to denote scenarios where the average noise incidence is 8 times that of the average measles incidence.
For the static noise scenarios, independent draws from a Poisson distribution with mean $c dot.op angle.l Delta I_M angle.r$ were simulated to produce the noise time series i.e., $Lambda (c) = upright("Pois")(c dot.op angle.l Delta I_M angle.r)$.
For the dynamical noise scenarios, the rubella vaccination rate at birth was set to 85.38%, 73.83%, 50.88%, 27.89%, or 4.92% to produce equivalent values of $Lambda$ (to within 2 decimal places): $Lambda (c) = angle.l Delta I_R angle.r + upright("Pois")(0.15 dot.op angle.l Delta I_R angle.r)$.
We simulated 100 time series of 100 years for each scenario before summarizing the distributions of outbreak detection methods.

== Defining Outbreaks
It is common to use expert review to define outbreaks when examining empirical data, but this is not feasible in a modeling study where tens of thousands of years are being simulated.
Previous simulation studies define an outbreak as a period where $R_"t" > 1$ with the aim of detecting an outbreak during the grow period @jombartRealtimeMonitoringCOVID192021 @stolermanUsingDigitalTraces2023, or use a threshold of > 2 standard deviations (s.d.) over the mean seasonal incidence observed in empirical data (or from a 'burn-in' period of the simulation) @sternAutomatedOutbreakDetection1999 @salmonMonitoringCountTime2016 @teklehaimanotAlertThresholdAlgorithms2004 @leclereAutomatedDetectionHospital2017.

Here we simulate time series of 100 years and we define a measles outbreak as a region of the time series that meets the following three criteria:

- The daily measles incidence must be greater than or equal to 5 cases
- The daily measles incidence must remain above 5 cases for greater than or equal to 30 consecutive days
- The total measles incidence must be great than or equal to 500 cases within the bounds of the outbreak

Only events meeting all 3 criteria are classified as outbreaks.
The incidence of non-measles febrile rash (i.e., noise) does not affect the outbreak status of a region but may affect the alert status triggered by the testing protocol.

Each day, a percentage (P) of clinically-compatible cases of febrile rash are tested; P is fixed in a given scenario to a value between 10% and 100%, in 10% increments.
Each "testing scenario" combines a testing rate (P) with one of the following tests:

- An imperfect test with 85% sensitivity and specificity, and 0-day lag in result return. That is, 85% of true measles cases will be correctly labeled as positive, and 15% of non-measles febrile rash individuals that are tested will be incorrectly labeled as positive for measles. This acts as a lower bound of acceptability for a hypothetical measles RDT @20240613_tpp_measles_rubell_FV_EN
- An imperfect test with 90% sensitivity and specificity, and 0-day lag in result return @brownRapidDiagnosticTests2020
- A perfect test with 100% sensitivity and specificity, and a 0-day test result delay. This is more accurate than is observed for current ELISA tests @hiebertEvaluationDiagnosticAccuracy2021, but it used to evaluate the theoretical best-case scenario
- A perfect test with 100% sensitivity and specificity, and a 14-day test result delay that represents a best-case test under more realistic reporting delays in result return

For each time series of true measles cases, we define outbreaks as the range of time that meets the definition above (@fig-outbreak-schematic\a).
We then add non-measles noise (@fig-outbreak-schematic\b) and test according to the testing scenario, which yields 5 time series of test positive cases (@fig-outbreak-schematic\c): one time series of all clinically compatible cases and 4 reflecting the testing scenarios.

#figure(
  image("manuscript_files/plots/schematic-plot.svg"),
  caption: [
    A schematic of the outbreak definition and alert detection system. A) Measles incidence time series. B) Noise incidence time series. C) Observed time series of test positive cases according to a given testing scenario. The orange bands present in all 3 panels represent regions of the measles time series that meet the outbreak definition criteria. In panel C, the dark blue bands represent regions of the test positive time series that breach the alert threshold (the horizontal dashed line), and constitute an alert.
  ],
)
<fig-outbreak-schematic>

== Triggering Alerts
We define an "alert" as any consecutive string of 1 or more days where the 7-day (trailing) moving average of the test positive cases is greater than or equal to a pre-specified alert threshold, T.
For each time series of test positive cases, we calculate the percentage of alerts that are "correct", defined as any overlap of 1 or more days between the outbreak and alert periods (@fig-outbreak-schematic\c).
This is analogous to the PPV of the alert system and will be referred to as such for the rest of the manuscript.
It is possible to have multiple alerts within a single outbreak if the 7-day moving average of test positive cases drops below, and then re-crosses, the threshold, T, and we count each as correct.
For all outbreaks in the measles time series, we calculate the percentage that contain at least 1 alert within the outbreak’s start and end dates (@fig-outbreak-schematic\c).
We refer to this as the sensitivity of the alert system.
We also calculate the detection delay as the time from the start of an outbreak to the start of its first alert.
If the alert period starts before the outbreak and continues past the start date of the outbreak, this would be considered a correct alert with a negative delay i.e., an early warning triggered by false positive test results.
Finally, for each time series we calculate the number of unavoidable and avoidable outbreak cases.
Unavoidable cases are those that occur before a correct alert, or those that occur in an undetected outbreak.
Avoidable cases are those that occur within an outbreak and after the first alert; we do not quantity avoidable cases here as the value depends critically on the explicit details of the response, which we do not model.

We define the accuracy of the surveillance system for a given time series as the arithmetic mean of the system’s PPV and sensitivity.
To examine the interaction of the test with the surveillance system's characteristics (i.e., testing rate, noise structure and magnitude), we optimized the alert threshold, T, for a given "testing scenario".
Each of the 100 simulations per scenario produces an accuracy, and we identified the optimal alert threshold, T#sub([O]), as the value that produced the highest mean accuracy for a given scenario.
To identify T#sub([O]) we implemented the TikTak multistart optimization algorithm @arnoudBenchmarkingGlobalOptimizers2023, using 100 initial values (alert thresholds) selected from a Sobol' low-discrepancy sequence @sobolDistributionPointsCube1967 initialized with lower and upper bounds of 0.0 and 50.0, respectively.
In brief, the Sobol' sequence is a deterministic, quasi-random sequence of numbers that maximizes the uniformity of the explored parameter space by approximately iteratively bisecting the parameter space @sobolDistributionPointsCube1967 @lemieuxQuasiMonteCarlo2009.
After 100 initial alert thresholds are generated, the accuracy is evaluated and the 10 alert thresholds (points) with the highest accuracy are retained.
The 10 retained alert thresholds are sorted in descending order of accuracy, creating the sequence of Sobol' points ($upright(bold(s#sub[1])) dots upright(bold(s#sub[10]))$) that are used to calculate the seed points for local optimization that subsequently performed using the BOBYQA derivative-free algorithm @powellBOBYQAAlgorithmBound.
For each of the 10 local optimizations, the starting seed is computed as the weighted combination of the Sobol' point $upright(bold(s#sub[i]))$ and the alert threshold that produced the maximum accuracy so far, with increasing weight provided to the alert threshold that maximized accuracy; more information can be found in Appendix B.6 of @arnoudBenchmarkingGlobalOptimizers2023.
The TikTak algorithm is implemented in the #link("https://github.com/tpapp/MultistartOptimization.jl")[MultistartOptimization.jl] package @pappTpappMultistartOptimizationjl2025, with local optimization (BOBYQA) implemented in the #link("https://github.com/jump-dev/NLopt.jl")[NLOpt.jl] package @johnsonNLoptNonlinearoptimizationPackage2025.

We then compare testing scenarios at their respective optimal alert threshold.
This allows for conclusions to be made about the surveillance system as a whole, rather than just single components.
We also present results for optimizations based upon the harmonic mean (F-1 score) of the system's PPV and sensitivity in the Supplement (@fig-f1-alert-threshold, @fig-f1-accuracy, @fig-f1-delay, @fig-f1-alert-proportion).


#pagebreak()

= Author Contributions
#emph[Conceptualization:] all authors

#emph[Data curation:] CA, MJF

#emph[Formal analysis:] CA

#emph[Funding acquisition:] ACK, AW, BP, MJF, WJM

#emph[Investigation:] CA

#emph[Methodology:] CA, MJF

#emph[Project administration:] MJF

#emph[Software:] CA

#emph[Supervision:] AW, BP, MJF, WJM

#emph[Validation:] CA

#emph[Visualization:] CA

#emph[Writing - original draft:] CA

#emph[Writing - review and editing:] all authors

== Conflicts of Interest and Financial Disclosures
The authors declare no conflicts of interest.

== Data Access, Responsibility, and Analysis
Callum Arnold and Dr. Matthew J. Ferrari had full access to all the data in the study and take responsibility for the integrity of the data and the accuracy of the data analysis. Callum Arnold (Department of Biology, Pennsylvania State University) conducted the data analysis.

== Data Availability
All code and data for the simulations can be found at #link("https://github.com/arnold-c/OutbreakDetection")

#pagebreak()

#set bibliography(style: "elsevier-vancouver")

#bibliography("./OD.bib")
