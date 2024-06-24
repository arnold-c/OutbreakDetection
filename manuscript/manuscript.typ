#set page(paper: "us-letter", numbering: "1")
#set par(justify: true)
#set text(
  font: "Linux Libertine",
  size: 10pt,
)
#show figure.caption: emph
#set table(
  fill: (x, y) => {
    if y == 0 {gray}
  }
)
#show table.cell: it => {
  if it.y == 0 {
    strong(it)
  } else {
    it
  }
}
#set math.equation(numbering: "1.")


#let title = [The Need to Develop a Holistic Infectious Disease Surveillance System]
#let authors = [Callum R.K. Arnold#super[a, b], Alex C. Kong#super[c], Amy K. Winter#super[d], William J. Moss#super[c, e], Bryan N. Patenaude#super[c], Matthew J. Ferrari#super[a, b]]

#align(center, text(14pt)[
  *#title*
])

#align(center, text(10pt)[
  #authors
])

#linebreak()

#super[a] Department of Biology, Pennsylvania State University, University Park, PA, USA 16802

#super[b] Center for Infectious Disease Dynamics, Pennsylvania State University, University Park, PA, USA 16802

#super[c] Department of International Health, Johns Hopkins Bloomberg School of Public Health, Baltimore, MD, USA 21205

#super[d] Department of Epidemiology, College of Public Health, University of Georgia, Athens, GA, USA 30602

#super[e] Department of Epidemiology, Johns Hopkins Bloomberg School of Public Health, Baltimore, MD, USA 21205

#linebreak()

Corresponding author: Callum R.K. Arnold, _cfa5228\@psu.edu_, Center for Infectious Disease Dynamics, Pennsylvania State University, University Park, PA, USA 16802

#pagebreak()

#set page(
  header: align(right)[#title]
)

= Abstract
== Background

Infectious disease surveillance and outbreak detection systems often utilize diagnostic testing to validate case identification.
The metrics of sensitivity, specificity, and positive predictive value are commonly discussed when evaluating the performance of diagnostic tests, and to a lesser degree, the performance of outbreak detection systems.
However, the interaction of the two levels’ (the test and the alert system) metrics, is typically overlooked.
Here, we describe how equivalent regions of detection accuracy can exist over a range of diagnostic test characteristics, examining the sensitivity to background noise structure.

== Methods

We generated a stochastic SEIR model with importation to simulate true measles and non-measles sources of febrile rash (noise) daily incidence.
We simulated the noise to have in-phase, out-of-phase, non-seasonal dynamics, with respect to the measles simulation.
We also examined noise where each day were independent Poisson draws.
For rapid diagnostic test (RDT) and enzyme-linked immunosorbent assay (ELISA) like tests, we varied the reporting rate, and the observed daily incidence threshold required to trigger an alert, and computed the: percentages of alerts that were correct (analogous to the positive predictive value), the percentage of outbreaks detected (analogous to the sensitivity), the delay between outbreak initiation and alert, and the number of (un)avoidable outbreak cases.
Summarizing the percentage of alerts correct and outbreaks detected into an accuracy metric, we selected the “optimal threshold” for each combination of diagnostic test and reporting rate, presenting the summary statistics above.

== Results

The optimal threshold for each test increased monotonically as the percentage of clinic visited who were tested increased, though there are adjacent reporting rates that have identical thresholds due to the integer requirement for the daily threshold.
With Poisson-only noise, similar outbreak detection accuracies could be achieved with imperfect RDT-like tests as with ELISA-like diagnostic tests (c.
93%), given moderately high testing rates.
With larger delays (14 days) between the ELISA test procurement and result date, RDTs could outperform the ELISA.
Similar numbers of unavoidable cases and outbreak alert delays could be achieved between the test types.
With dynamical, in-phase, noise, however, the accuracy of ELISA scenarios was far superior to those achieved with RDTs (c.
93% vs. 73%). For dynamical noise, RDT-based scenarios typically favored more sensitive alert threshold than ELISA-based scenarios (at a given reporting rate), observed with lower numbers of unavoidable cases and detection delays.

== Conclusions

The performance of an outbreak detection system is highly sensitive to the structure and the level of background noise.
Under the assumption that the noise is relatively static over time (as influenza-like illness was early in the COVID-19 pandemic), RDTs can perform as well as ELISA in a surveillance system.
However, when the noise is non-stationary, imperfect tests cannot overcome their accuracy limitations through higher testing rates. 

#pagebreak()

= Background

At the heart of an outbreak detection system is a surveillance program built upon case detection, often utilizing individual diagnostic tests as necessary components of epidemiological investigations before the declaration of an outbreak @murrayInfectiousDiseaseSurveillance2017 @zhouEliminationTropicalDisease2013 @pahoIntegratedApproachCommunicable2000 @worldhealthorganizationMeaslesOutbreakGuide2022 @craggOutbreakResponse2018.
When diseases have non-specific symptoms, without accurate measurement tools we cannot have confidence in our ability to correctly ascribe changes in symptoms within a population to a particular disease, and therefore detect outbreaks of specific pathogens.
As a result, it has been commonplace for surveillance systems to be developed around high-accuracy tests, such as Polymerase Chain Reaction (PCR) tests and immunoglobulin (Ig) tests, when financially and logistically feasible @gastanaduyMeasles2019 @commissionerCoronavirusCOVID19Update2020 @grasslyComparisonMolecularTesting2020 @ezhilanSARSCoVMERSCoVSARSCoV22021 @worldhealthorganizationCholera2023 @essentialprogrammeonimmunizationepiimmunizationClinicalSpecimensLaboratory2018.
Depending on the disease in question, either sensitivity (the ability to correctly detect a true positive individual) or specificity (the ability to correctly discount a true negative) will be prioritized, as they often are at odds with each other @westreichDiagnosticTestingScreening2019 @shrefflerDiagnosticTestingAccuracy2024 @parikhUnderstandingUsingSensitivity2008.
This balance is commonly defined within the Target Product Profile (TPP) of a test @worldhealthorganizationTargetProductProfiles, which is a set of minimum characteristics that should be met for production and widespread use, helping to guide research and development.
For example, in the wake of the 2013 Ebola outbreak in Guinea a TPP was developed that listed the minimum acceptable sensitivity of 95% and specificity of 99% @chuaCaseImprovedDiagnostic2015.
Recognizing that Ebola is not the major cause of fever and other non-specific symptoms in the region, it is arguably more important to prioritize the specificity of the disease, although the authors note that the severity of the infection requires a high level of sensitivity as the consequences of a missed case are dire at an individual and population level @chuaCaseImprovedDiagnostic2015.

Much like the accuracy of an individual test, outbreak detection systems face the same issue regarding the prioritization of sensitive vs specific alerts @germanSensitivityPredictiveValue2000 @worldhealthorganizationOperationalThresholds2014 @lewisTimelyDetectionMeningococcal2001.
For many disease systems, particularly in resource constrained environments where the burden of infectious diseases is highest @gbd2019childandadolescentcommunicablediseasecollaboratorsUnfinishedAgendaCommunicable2023 @roserBurdenDisease2023, cases are counted and if they breach a threshold, be that weekly, monthly, or some combination of the two, an alert is triggered that may launch a further investigation and/or a response @worldhealthorganizationMeaslesOutbreakGuide2022 @worldhealthorganizationOperationalThresholds2014.
In effect, this discretizes a distinctly continuous phenomenon into a binary measure, outbreak or no outbreak, for decision making purposes.
When considering reactive vaccination campaigns and non-pharmaceutical based interventions designed to reduce transmission and limit and suppress outbreaks through direct and indirect protection, early action has the potential to avert the most cases, and in some instances, can be the most cost effective as it provides the opportunity to incorporate later information within an “adaptive management” framework @atkinsAnticipatingFutureLearning2020 @taoLogisticalConstraintsLead @graisTimeEssenceExploring2008 @ferrariTimeStillEssence2014 @worldhealthorganizationConfirmingInvestigatingManaging2009 @minettiLessonsChallengesMeasles2013.
While this framing would point towards a sensitive surveillance system being optimal, with the benefits clearly defined in terms of morbidity and mortality avoided, each action comes with both direct and indirect financial and opportunity costs stemming from unnecessary activities that remove resources from future response activities.
Just as the balance of sensitivity and specificity of a test must be carefully evaluated, so must the balance at the outbreak level.

The concept of using incidence-based alert triggers to define the discrete event of an “outbreak” with characteristics analogous to individual tests has been well documented in the case of meningitis, measles, and malaria @worldhealthorganizationMeaslesOutbreakGuide2022 @lewisTimelyDetectionMeningococcal2001 @worldhealthorganizationConfirmingInvestigatingManaging2009 @trotterResponseThresholdsEpidemic2015 @cooperReactiveVaccinationControl2019 @zalwangoEvaluationMalariaOutbreak2024 @kanindaEffectivenessIncidenceThresholds2000.
However, an overlooked, yet critical, aspect of an outbreak detection system is the interplay between the individual test and outbreak alert characteristics.
With their success within malaria surveillance systems, and particularly since the COVID-19 pandemic, rapid diagnostic tests (RDTs) have garnered wider acceptance, and their potential for use in other disease systems has been gaining interest @warrenerEvaluationRapidDiagnostic2023.
Despite concerns of their lower diagnostic accuracy slowing their adoption until recently @millerAddressingBarriersDevelopment2015, the reduced cold-chain requirements @brownRapidDiagnosticTests2020, and gains in speed of result and the impact on the timeliness of an alert system may outweigh the costs in some settings @warrenerEvaluationRapidDiagnostic2023 @mcmorrowMalariaRapidDiagnostic2011 @larremoreTestSensitivitySecondary2021 @middletonModelingTransmissionMitigation2023.

In this paper we examine how the interaction between imperfect diagnostic tests affects the optimal outbreak alert method for the detection of measles outbreaks, in the context of a febrile rash surveillance system that detects and tests suspected measles cases i.e., an individual with fever and maculopapular (non-vesicular) rash (or who the healthcare worker suspects has measles) @essentialprogrammeonimmunizationepiimmunizationvaccinesandbiologicalsivbMeasles2018.
Because measles symptoms are non-specific, it is important to account for non-measles sources of febrile rash e.g., rubella, parvovirus, varicella, etc., producing the potential for false positive results in the context of imperfect tests.
Currently, measles outbreaks are declared on the basis of either suspected measles cases alone, the use of an enzyme-linked immunosorbent assay (ELISA) to detect the presence of measles-specific IgM antibodies, or a combination of the two, depending on the elimination status of the region, with countries nearer elimination increasing the use of ELISA diagnostic confirmation for suspected cases @essentialprogrammeonimmunizationepiimmunizationvaccinesandbiologicalsivbMeasles2018.
Each of these detection systems have their flaws.
Although clinical case definition is very fast and requires limited resources, providing a highly sensitive alert system, in the face of high “background noise” from non-measles sources of febrile rash relative to true measles prevalence, it can lead to great uncertainty in the alert system and low positive predictive value (PPV), i.e., the probability that an alert accurately captures an outbreak @hutchinsEvaluationMeaslesClinical2004.
And while ELISA confirmation is the standard diagnostic test for measles surveillance @essentialprogrammeonimmunizationepiimmunizationvaccinesandbiologicalsivbMeasles2018, it is expensive due to its cold-chain and transport costs and delays resulting from the requirement to transfer samples from point of care to a central laboratory, leading to a slower and less sensitive (but more specific) alert system.
A rapid diagnostic test offers the opportunity to provide a compromise between diagnostic accuracy and timeliness, and recent developments show encouraging signs in the field as well as in theory @warrenerEvaluationRapidDiagnostic2023 @brownRapidDiagnosticTests2020.

By examining the combination of alert threshold and individual test characteristic in a modeling study that explicitly incorporates dynamical background noise, we aim to illustrate the need to develop a TPP for the whole detection system, not just one component.
To evaluate the alert system performance, we develop a set of outbreak definition criteria and surveillance metrics, drawing inspiration from acceptance sampling, ecological surveillance systems, and epidemiological surveillance systems guidelines and reviews @yemshanovAcceptanceSamplingCosteffective2020 @christensenHerdlevelInterpretationTest2000 @muratoEvaluationSamplingMethods2020 @sternAutomatedOutbreakDetection1999 @calbaSurveillanceSystemsEvaluation2015 @guidelinesworkinggroupUpdatedGuidelinesEvaluating2001.
Using these metrics we overcome issues encountered by early warning systems that rely on dynamical values such as $R#sub[effective]$ in defining outbreaks @jombartRealtimeMonitoringCOVID192021 @proverbioPerformanceEarlyWarning2022 @stolermanUsingDigitalTraces2023 @brettAnticipatingEpidemicTransitions2018 @salmonMonitoringCountTime2016, for example characterizing the end of an epidemic period is important in a timeseries where multiple outbreaks will occur.

= Methods
== Model Structure

We constructed a stochastic compartmental non-age structured SEIR model of measles and simulated using a modified Tau-leaping algorithm with a time step of 1 day, utilizing binomial draws to ensure compartment sizes remained positive valued @chatterjeeBinomialDistributionBased2005 @gillespieApproximateAcceleratedStochastic2001.
To stochastically reintroduce infections and cause recurring epidemics, commuter-style imports were added that are proportional to the size of the population and $R#sub[0]$ @keelingModelingInfectiousDiseases2008, and the transmission parameter (#sym.beta#sub[t]) is sinusoidal with a period of one year.
$R#sub[0]$ was set to 16, with a latent period of 10 days and infectious period of 8 days. The population was initialized with 500_000 individuals with Ghana-like birth and vaccination rates, and the final results were scaled up to the 2022 population size of Ghana (33 million). The full table of parameters can be found in @table_model-parameters. All simulations and analysis was completed in Julia version 1.10.3 @bezansonJuliaFreshApproach2017, with all code stored at _ https://github.com/arnold-c/OutbreakDetection _.

#let parameters = csv("parameters.csv")
#let import_rate = $(1.06*μ*R_0)/(√(N))$
#let parameter_labels = ( "Parameters", $R#sub[0]$, $"Latent period ("#sym.sigma")"$, $"Infectious period ("#sym.gamma")"$, "Seasonal amplitude", $"Birth rate ("#sym.mu")"$, $"Vaccination rate at birth ("#sym.rho")"$)

#figure(
	table(
		columns: 3,
		..for ((label), (.., measles, nonmeasles)) in parameter_labels.zip(parameters) {
			(label, measles, nonmeasles)
		},
		[Importation rate], table.cell(colspan: 2, align: center, $(1.06*μ*R_0)/(√(N))$),
		[Population size (N)], table.cell(colspan: 2, align: center, "500_000, scaled to 33M"),
		[Initial proportion susceptible], table.cell(colspan: 2, align: center, "0.05"),
		[Initial proportion exposed], table.cell(colspan: 2, align: center, "0.0"),
		[Initial proportion infected], table.cell(colspan: 2, align: center, "0.0"),
		[Initial proportion recovered], table.cell(colspan: 2, align: center, "0.95"),
	),
	caption: "Compartmental model parameters"
)<table_model-parameters>

To examine the sensitivity of the detection system to background noise, we layered the measles incidence timeseries with one of four noise time series: Poisson-only noise, or dynamical noise with rubella-like parameters that could be in- or out-of-phase or with independent seasonality to the measles dynamics.
The total annual incidence of noise in each system was approximately 8 times the average annual measles incidence, with the dynamical noise systems also containing low levels of Poisson noise, reflecting Ghana-like dynamics.
Despite the average incidence of noise in each scenario being the same, the relative proportion of measles to noise varies differently throughout the timeseries, exacerbating the effects of imperfect diagnostic tests that will produce higher rates of false positives and negatives than ELISA-like diagnostics.
Imperfect tests are also more susceptible to incorrect results because of the varying prevalence throughout the timeseries, so the PPV of the system will not be static throughout the timeseries.
To account for the stochasticity of the model, we simulated 100 time series for each scenario, each lasting 100 years, before summarizing the distributions of outbreak detection methods.

== Defining Outbreaks

It is common to use expert review to define outbreaks when examining empirical data, but this is not feasible in a modeling study where tens of thousands of years are being simulated.
To account for this, many studies only simulate a single outbreak within a timeseries (repeating this short stochastic simulation multiple times to ensemble results), define an outbreak as a period where $R#sub[effective]$ > 1, or use a threshold of > 2 standard deviations (s.d.) over the mean seasonal incidence observed in empirical data @sternAutomatedOutbreakDetection1999 @jombartRealtimeMonitoringCOVID192021 @stolermanUsingDigitalTraces2023 @salmonMonitoringCountTime2016 @teklehaimanotAlertThresholdAlgorithms2004 @leclereAutomatedDetectionHospital2017.
Each method has its uses, but to evaluate the performance of an outbreak detection system in an endemic region where multiple sequential epidemics are expected it is important to clearly define the bounds of the outbreak, which can only be achieved by one of these methods (2 s.d. > mean).
This, however, assumes strong seasonal forcing and regular periodicity of incidence to produce a smooth enough baseline, which is not present as countries near measles elimination status @grahamMeaslesCanonicalPath2019.
Here we define a true measles outbreak as a region of the timeseries that meets the following three criteria:

+ The daily measles incidence must be greater than, or equal to, 5 cases
+ The daily measles incidence must remain above 5 cases for greater than, or equal to, 30 consecutive days
+ The total measles incidence must be great than, or equal to, 500 cases within the bounds of the outbreak

Each of these is a necessary, but not sufficient, condition to the definition of an outbreak; all must be met.
For example, a region may produce 1000 total cases within the lower and upper bounds defined by the daily incidence being #sym.gt.eq 5 cases, but if this period only lasts for 20 days it would not be considered an outbreak as it would not be possible to mount a response to.
The incidence of non-measles febrile rash (i.e., noise) does not affect the outbreak status of a region but does affect the alert status triggered by the testing protocol.

== Triggering Alerts

The use of alert thresholds is common in locations burdened by measles, where access to high performance computing is limited making nowcasting-style approaches to outbreak detection impractical.
For this reason, we define an “alert” as any single day where the 7-day moving average of the daily incidence is greater than, or equal to, a threshold of X.
To examine the interaction of test and surveillance systems we vary the alert threshold, X, between 1 and 15 cases per day, before computing the evaluation metrics and identifying the optimal threshold for each individual test.
We also examined an alert method that required either the daily incidence or the 7-day moving average of daily incidence to exceed the alert threshold X (see supplemental results).

Each day, 60% of the measles and non-measles febrile rash cases visit the clinic for treatment, and Y% percentage of these clinic visits are tested, as all clinic visits are deemed to be suspected measles cases because they meet the clinical case definition.
This percentage of clinic visits, Y, that are tested is varied between 10% and 60%, in 10% increments, for all combinations of diagnostic test (apart from the clinical case definition) and alert threshold, defining the “testing scenario”.
While a reporting rate of 6-48% is high, we are only simulating a population of 500_000 individuals i.e., approximately 1.5% of Ghana’s population, which produces absolute incidence values that are consistent with observed reporting on the national scale. Each testing scenario uses one of the following tests:

+ No test, where every suspected measles case is counted as positive for measles and used to trigger an outbreak alert, which is equivalent to the use of the clinical case definition. This is implemented as a diagnostic test with 100% sensitivity, 0% specificity, a 0-day lag in returning the result, and all individuals who visit the clinic are tested (not varied between 10-60%, like the other diagnostic tests described).
+ An RDT equivalent with 85% sensitivity and specificity, and 0-day lag in result return. That is, 85% of true measles cases will be correctly labelled as positive, and 15% of non-measles febrile rash individuals that are tested will be incorrectly labelled as positive for measles. This acts as a lower bound of acceptability for a new measles RDT.
+ An ELISA equivalent test with 100% sensitivity and specificity with a 0-day test result delay.
+ An ELISA equivalent test with 100% sensitivity and specificity with a 7-day test result delay.
+ An ELISA equivalent test with 100% sensitivity and specificity with a 14-day test result delay.

== Optimal Thresholds

For each of the previously defined tests, we calculate the number of test positive cases (that will include false positive and negative cases resulting from imperfect diagnostic tests), which is used to categorize the time series by alert status in conjunction with the alert threshold X.
The overlap with true outbreak status, or lack thereof, is used to classify the following metrics, and has been illustrated in @figure_outbreak-schematic:

+ The percentage of alerts that are correct. This is analogous to the PPV of the alert system. Note that it is possible to trigger multiple alerts within a single outbreak, and each would be considered correct.
+ The percentage of outbreaks that are detected. This is analogous to the sensitivity of the alert system.
+ The detection delay. Note that if an alert precedes the start of the outbreak, as long as it successfully “captures” the outbreak within its bounds, it is considered to be correct, resulting in a negative detection delay i.e., an early warning triggered by false positives.
+ The number of unavoidable and avoidable outbreak cases and deaths. Avoidable cases are defined as those that occur within an outbreak after a correct alert is first triggered i.e., cases that could theoretically be prevented with a perfectly effective and timely response. Unavoidable cases are the inverse: those that occur before a correct alert, or those that occur in an undetected outbreak. Using a fitted age distribution of measles cases in Ghana @mintaProgressMeaslesElimination2023 and the age-specific case fatality ratio for Ghana in 2022 @sbarraEstimatingNationallevelMeasles2023, we scaled the values for each one-year age cohort, per annum. In practice, not all cases deemed avoidable are (due to imperfect and delays in responses), but to minimize the sensitivity of the results to the response implementation and operational constraints we are counting them as such.

#figure(
	image("./schematic.png", width: 85%),
	caption: "A schematic of the outbreak definition and alert detection system"
)<figure_outbreak-schematic>

For each combination of diagnostic test and testing rate, the optimal alert threshold is calculated by selecting the threshold that produces the highest accuracy.
Each of the 100 simulations produces an accuracy (the mean of the percentages of alerts that are correct and the percentage of the outbreaks that are detected), and the median accuracy of this distribution is used to determine the optimal alert threshold.
Each testing scenario (the combination of the diagnostic test, testing rate, and alert threshold), is then compared (using the metrics above) at the optimal threshold that is specific for the scenario.
This allows for conclusions to be made about the surveillance system as a whole, instead of just single components.

== Cost-Benefit Analysis

To quantify the relative benefit of each testing scenario at their optimal threshold, we calculated the Incremental Cost-Effectiveness Ratios (ICERs) using 1) the number of diagnostic tests used in each scenario, 2) the estimated cost of a single diagnostic test (accounting for transport and other logistical costs), and 3) the cost of responding to all alerts within a year.
If a testing scenario produced clinical benefits e.g., fewer unavoidable cases and smaller alert delay, and cost less than \$2,203.56 per DALY averted (the 2022 GDP per capita for Ghana @worldbankGhana), it was deemed to be cost effective.

$ "ICER" = (("Total cost of alternative") - ("Total cost of reference")) / (("Unavoidable DALYs for reference") - ("Unavoidable DALYs for alternative")) $<icer>

Disability-Adjusted Life Years (DALYs) are calculated by summing the discounted Years of Life Lost (YLLs) to premature mortality (Equation 2) with the Years Lost due to Disability (YLDs) (Equation 3).
YLLs included a discounting factor of 3% per year, and used Ghana’s 2019 life table from the WHO @worldhealthorganizationLifeTablesLife
Disability weights of 0.051 and 0.133 were used for moderate and severe measles cases, respectively, with a duration of illness of 10 days, assuming 50% of cases were severe and the remaining moderate @vosGlobalBurden3692020.

$ "YLL"#sub[p] = sum_(n=0)^("LE"#sub[p]) (1/1.003^n) $

#align(center, [
YLL#sub[p] = total discounted years of life lost for an individual who dies at age p \
LE#sub[p] = expected years of remaining life for an individual who dies at age p  \
n = the nth year of additional life an individual was expected to live
])

#linebreak()

$ "YLD" = sum_(s) "I"#sub[s] * "DW"#sub[s] * "d"#sub[s] $

#align(center, [
	I#sub[s] = number of incident cases for disease severity s \
	DW#sub[s] = disability weighing for disease severity s \
	d#sub[s] = duration of illness for disease severity s
])

= Results
== Optimal Thresholds, Accuracy, and Unavoidable Cases

The threshold that maximized outbreak detection accuracy given a testing scenario depends heavily on the noise structure.
For example, the optimal threshold for an RDT with 90% sensitivity/specificity and 40% testing of clinic visits was 4 cases per day with the dynamical in-phase noise, whereas it is 5 cases per day when Poisson-only noise was simulated (Table 2).
This corresponds to an accuracy of 72% and 93%, respectively (Table 3).
Much higher accuracy can be achieved with imperfect tests under Poisson noise because large spikes of non-measles febrile rash do not occur that would lead to many false positives, triggering an erroneous alert.
With dynamical noise, this possibility can and does occur, and when the ratio of measles to non-measles cases is sufficiently low, drastically impacts the accuracy of an alert system.
This is why the values observed among ELISA tests are identical between the different noise structures: with a perfect test there are no false positive results.
For dynamical noise, RDTs never perform as well as an ELISA, even scenarios that include a 14-day lag between test procurement and result receipt.
However, with Poisson noise, once testing reaches 20% of clinic visits, equivalence can be observed at approximately 91% accuracy (Table 3).


//#let optimal_thresholds = csv("optimal-thresholds.csv")
//
//#table(
//	columns: 10,
//	..optimal_thresholds.flatten()
//)

#table(
	columns: 2,
	[Test], [Second columns],
	[], []
)

_Table 2) Optimal threshold for each testing scenario. A) the noise structure is dynamical, and the seasonality is in-phase with the measles incidence. B) the noise structure is Poisson only_

Moving across a single row of Table 2, as the testing rate increases, so does the optimal threshold, because more individuals are testing positive.
However, the changes in accuracy are not monotonically increasing, as one might anticipate.
Even with a perfect test, there are sequential cells that increase the testing rate, yet accuracy decreases. For example, the ELISA equivalent with a 0-day test result lag goes from an outbreak detection accuracy of 91% down to 87% when testing increases from 30% to 40% of clinic visits.
The reason behind this unintuitive result stems from the use of integer valued thresholds.
While computing the moving average incidence helps to smooth out the alert system, requiring the daily alert to be an integer value results in the optimal threshold for alert accuracy being stuck at a lower level than a truly optimal intermediate floating point number, because the increase in the alert system’s PPV from a higher threshold is outweighed by the larger decrease in alert sensitivity, or the threshold is forced to increase.
Examining the number of unavoidable cases (Table 4) and the detection delays (Table 5), we can see that the decrease in accuracy results from a decrease in the sensitivity of the alert system: the detection delay increases by 6 days, and the unavoidable cases by c. 2300.
Although the testing rate increases, so does the optimal threshold, from 3 test positives to 4 test positives per day.

_Table 4) Unavoidable cases per annum of each testing scenario at their specific optimal thresholds, scaled up to Ghana’s 2022 population. A) the noise structure is dynamical, and the seasonality is in-phase with the measles incidence. B) the noise structure is Poisson only._

	Test Characteristic	% of Clinic Visits Tested
	Testing Scenario	Test lag	10	20	30	40	50	60	100
Dynamical noise: in-phase	Clinical Case Definition	0							20643
	RDT Equivalent 0.85	0	436	3070	4510	8102	3824	4666
	RDT Equivalent 0.90	0	488	3240	4618	6463	2671	5092
	ELISA Equivalent	0	770	5980	8893	11172	4529	6144
	ELISA Equivalent	3	992	6670	9660	11949	5164	6858
	ELISA Equivalent	7	1293	7572	10619	12870	6007	7761
	ELISA Equivalent	14	2015	9277	12363	14643	7641	9495
Poisson Noise	Clinical Case Definition	0							53419
	RDT Equivalent 0.85	0	766	6592	9111	9107	11865	11765
	RDT Equivalent 0.90	0	770	3178	9736	12808	5205	8111
	ELISA Equivalent	0	770	5980	8893	11172	4529	6144
	ELISA Equivalent	3	992	6670	9660	11949	5164	6858
	ELISA Equivalent	7	1293	7572	10619	12870	6007	7761
	ELISA Equivalent	14	2015	9277	12363	14643	7641	9495

_Table 5) Outbreak alert delay (days) of each testing scenario at their specific optimal thresholds. A) the noise structure is dynamical, and the seasonality is in-phase with the measles incidence. B) the noise structure is Poisson only._

	Test Characteristic	% of Clinic Visits Tested
	Testing Scenario	Test lag	10	20	30	40	50	60	100
Dynamical noise: in-phase	Clinical Case Definition	0							-4674.81
	RDT Equivalent 0.85	0	-28.07	-15.53	-11.84	-5.76	-18.47	-13.11
	RDT Equivalent 0.90	0	-19.85	-7.56	-3.55	0.37	-13.10	-5.74
	ELISA Equivalent	0	-3.69	22.61	30.64	36.08	17.92	23.05
	ELISA Equivalent	3	-2.64	24.55	33.04	38.83	19.55	24.94
	ELISA Equivalent	7	-0.46	27.85	36.51	42.09	22.67	28.23
	ELISA Equivalent	14	3.69	33.64	42.38	48.18	28.21	34.07
Poisson Noise	Clinical Case Definition	0							-2095.63
	RDT Equivalent 0.85	0	-3.75	24.49	31.45	31.85	38.17	37.87
	RDT Equivalent 0.90	0	-3.69	12.94	32.92	40.14	20.26	28.74
	ELISA Equivalent	0	-3.69	22.61	30.64	36.08	17.92	23.05
	ELISA Equivalent	3	-2.64	24.55	33.04	38.83	19.55	24.94
	ELISA Equivalent	7	-0.46	27.85	36.51	42.09	22.67	28.23
	ELISA Equivalent	14	3.69	33.64	42.38	48.18	28.21	34.07


== Cost-Benefit Analysis

…

= Discussion

The performance of an outbreak detection system is highly sensitive to the structure and level of background noise in the simulation.
Despite the mean daily noise incidence set to equivalent values between the dynamical and Poisson-only simulations, drastically different results are observed.
Under the assumption that non-measles febrile rash is relatively static in time (Poisson), RDTs can perform as well, if not better than ELISA tests at moderate to high testing rates, and at a fraction of the cost.
However, if it is expected that the noise is non-stationary, imperfect tests cannot overcome their accuracy limitations through higher testing rates, saturating at c. 74% accuracy, relative to ELISA’s 93%.

The purpose of routine surveillance is to characterize the infection landscape.
With strong public health infrastructure and infectious disease surveillance programs, it is possible to develop a strong understanding of the shape of febrile rash cases, regardless of source.
With this information, countries can tailor their future activities to rely more or less heavily upon RDTs, depending on the dynamics of the target disease and its relationship to background noise, favoring RDTs when there are low levels of noise and ELISAs during large rubella outbreaks, for example.

== Limitations and Strengths

To our knowledge, this is one of the first simulation studies to examine the relationship between individual test characteristics and the wider surveillance program.
By explicitly modeling the interaction between the two, we make a case that surveillance systems should take a holistic approach; ignoring one component can lead to drastically different, and suboptimal, results.
Additionally, by defining outbreak bounds concretely we have been able to calculate metrics of outbreak detection performance that draw parallels to those used when evaluating individual diagnostic tests, allowing for intuitive and simple implementation of this method in resource-constrained environments, something that is not possible with most outbreak detection and early warning system simulations in the literature.
An evaluation of all outbreak detection algorithms is beyond the scope of this work, but a more computationally expensive approach based on nowcasting incidence may help overcome the shortcomings of RDTs in high-noise scenarios.

For computational simplicity, this paper did not include demography in the model structure.
And while a simulation-based approach allows for complete determination of true infection status i.e., measle vs non-measles febrile rash cases, and therefore an accurate accounting of the outbreak and alert bounds, these simulations do not specifically represent any real-world setting.
The evaluation of empirical data does provide this opportunity, but at the cost of knowing the true infection status of individuals and confounding of multiple variables, limiting analysis to only those who are observed (i.e., not those in the community who do not visit a healthcare center), and removing the possibility to explore the sensitivity of the results to parameters of interest to a surveillance program e.g., testing rate, and the test itself.

Additionally, is has been well documented that the performance of an individual test is highly sensitive to its timing within a person’s infection cycle @gastanaduyMeasles2019 @larremoreTestSensitivitySecondary2021 @middletonModelingTransmissionMitigation2023 @kisslerViralDynamicsAcute2021 @ratnamPerformanceIndirectImmunoglobulin2000, so it is possible that different conclusions would be drawn if temporal information about the test acquisition was included in the simulation.

Finally, the optimal threshold for a testing scenario depends heavily on the costs ascribed to incorrect actions, be that failing to detect an outbreak or incorrectly mounting a response for an outbreak that doesn’t exist.
In the simulations we have weighted them equally, but it is likely that they should not be deemed equivalent: missing an outbreak leads to many thousands of cases and associated DALYs, whereas an unnecessary alert would generally launch an initial low-cost investigation for full determination of the outbreak status.
This is particularly important in countries with vast heterogeneity in transmission: different weightings should be applied to higher vs. lower priority/risk regions to account for discrepancies in consequences of incorrect decisions.
For example, a high priority zone could still benefit from a false alert as many high-risk healthcare regions within a country are targeted for Supplemental Immunization Activities, so a false alert would just hasten the (proactive) vaccination response.

Given these limitations, the explicit values (i.e., optimal thresholds, accuracies etc.) should be interpreted with caution, and the exact results observed in the real-world will likely be highly dependent on unseen factors, such as the proportion of measles and non-measles sources of febrile rash that seek healthcare.
However, the general patterns should hold, and more importantly, the analysis framework provides a consistent and holistic approach to evaluating the tradeoff between individual level tests and the alert system enacted to detect outbreaks. 

#pagebreak()

= Funding

This work was supported by funding from Gavi, the Vaccine Alliance.

= Acknowledgments

== Author Contributions

*Conceptualization:* CA, MJF

*Data curation:* MJF, CA

*Formal analysis:* CA, MJF

*Funding acquisition:* MJF, WM, AW

*Investigation:* CA, MJF

*Methodology:* CA, MJF

*Project administration:* MJF

*Software:* CA

*Supervision:* MJF, WM, AW, BP

*Validation:* CA, MJF

*Visualization:* CA

*Writing - original draft:* CA, MJF

*Writing - review and editing:* all authors.

== Conflicts of Interest and Financial Disclosures

The authors declare no conflicts of interest.

== Data Access, Responsibility, and Analysis

Callum Arnold and Dr. Matthew J. Ferrari had full access to all the data in the study and take responsibility for the integrity of the data and the accuracy of the data analysis. Callum Arnold and Dr. Matthew J. Ferrari (Department of Biology, Pennsylvania State University).

== Data Availability

The datasets generated during and/or analyzed during the current study are not publicly available due to containing personally identifiable information but are available from the corresponding author on reasonable request.

#pagebreak()

#bibliography("OD.bib")
