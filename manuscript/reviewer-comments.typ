#show heading.where(level: 1): it => align(center, it)
#show link: underline

// Function to format responses to reviewer comments
#let reviewer(content) = {
  block(
    inset: (left: 2em),
    text(fill: red, emph(content)),
  )
}

= Reviewer Comments

== General Response

We would like to thank the reviewers for their insightful comments; many of them resulted in substantial reflection of our modeling approach and assumptions and have led to improvements to the paper.
As a result of these comments we have refined the code base to not only increase its speed and modularity for clarity, but have also made a few changes to the logic.
These changes are detailed at the #link(<logic-changes>)[end of this document].
We have also make the following changes that were previously bugs in the code.
Similarly, these have been detailed at the #link(<bug-fixes>)[end of this document].

Throughout the rest of the document, all inset red text restates the reviewer comments and questions; our responses and descriptions are detailed in black text.

== Journal Requirements

#reviewer[
  1.
  We ask that a manuscript source file is provided at Revision.
  Please upload your manuscript file as a .doc, .docx, .rtf or .tex.
  If you are providing a .tex file, please upload it under the item type 'LaTeX Source File' and leave your .pdf version as the item type 'Manuscript'.
]

Completed #sym.dash.em uploaded as .docx.

#reviewer[
  2.
  Please upload all main figures as separate Figure files in .tif or .eps format.
  For more information about how to convert and format your figure files please see our guidelines: https://journals.plos.org/ploscompbiol/s/figures
]

Completed.

#reviewer[
  3.
  We have noticed that you have uploaded Supporting Information files, but you have not included a list of legends.
  Please add a full list of legends for your Supporting Information files after the references list.
]

Completed.

#reviewer[
  4.
  Please revise your current Competing Interest statement to the standard "The authors have declared that no competing interests exist."
]

Completed.

#reviewer[
  Note: If the reviewer comments include a recommendation to cite specific previously published works, please review and evaluate these publications to determine whether they are relevant and should be cited.
  There is no requirement to cite these works unless the editor has indicated otherwise.
]


== Reviewer \#1

#reviewer[
  Dear Authors,

  Thank you for the opportunity in reviewing this manuscript.
  The work is a consistent and detailed model of detectability of outbreak given different accuracy of tests, which is a crucial part in outbreak detection.

  Though interesting and worth publishing, I have minor comments that should be addressed.
]

=== Comment 1.1 (Line 100, Figure 1)

#reviewer[
  Figure 1 shows that perfect test with 0-day and 14-days delay performs exactly equal, which is not wrong as the figure only shows proportion of individual tested against alert threshold.
  It would be interesting to add or to point from that figure to a time series plot of each of the 4 different scenarios showing how the delay in testing would affect the timing of alert, which is of importance in outbreak detection algorithm.
  As supplementary material is okay.
]


Thank you for your suggestion.
Figure 1 illustrates that perfect tests with 0-day and 14-day delays achieve their optimal outbreak detection at the same alert threshold (for a given proportion of individuals tested),
and Figure 2 shows a slightly diminished detection accuracy with the 14-day delay tests.
We have included a comparison of two test-positive time series, resulting from the same measles and rubella simulation, in the supplement (S10 Figure).


=== Comment 1.2 (Figure 2 and Line 102)

#reviewer[
  The same figure 2 and for the text that follows it, around line 102.
  Would be great to have another version of Figure 3 where is clear shown that perfect tests and imperfect test separately.
  The difference in perfect test with different delay do not yield any clear or statistical significant difference of outbreak alert, not sure if I followed the explanation given to this.
  Though is clear from figure 3 that imperfect testing can be problematic triggering the surveillance system in unnecessary levels of alert and by chance or seasonality of the disease, such as is for measles.
]


Thank you for your comment.
As you mention, there is not a significant difference in outbreak detection accuracy (or in the other metrics).
As such, we feel that showing clear and substantial overlap is helpful in illustrating the point, and producing separate figures that only highlight the differences between perfect tests with 0-day and 14-day delays may lead the reader to believe there more of a difference in outcomes than present.
Similarly, given the 10-panel nature of the figures at present, separating the perfect and imperfect test results into (20) separate facets may be visually overwhelming, and make it harder to compare the relative performance of the tests.

Regarding the wording, we have expanded upon the explanation to clarify the reasoning behind test result lags affecting outbreak detection accuracy (lines 113-116).
We have also added a reference to Mina et al.'s perspective in lines 197-199.


=== Comment 1.3 (Line 173, Citations)

#reviewer[
  The phrase is not completely correct, Larremore et al.
  2021 that you cites is an example of modeling characteristics and its consequence for a wider surveillance program, as also Hay et al.
  2021.
  In a more general sense and targeted to non modeling expert, the perspective work of Mina et al.
  2020, entitled "Rethinking Covid-19 Test Sensitivity — A Strategy for Containment" addresses this, would recommend to reshape the phrase and cite those works.
]

Thank you for this comment - we have clarified this opening line to state that the novelty of our work is that it examines the interaction between the test and the background noise.


=== Comment 1.4 (Simulation Duration)

#reviewer[
  Minor comment, Why the simulated time series is as long as 100 years? In daily resolution? This is extremely long and somewhat unrealistic.
  As measles settings around the world varies a lot would be interesting to address some short-term outbreaks and low level of incidence, as some settings are close to elimination.
  An interesting development would be to understand how effectively different this 4 kind of testing scenarios can affect the elimination of measles.
  No action needed and not mandatory to be addressed, though would interesting to read some speculation on those scenarios in light of the work developed at the present manuscript.
  Regarding this, what are the time unit for figure 5?
]


We simulated 100 years per time series to ensure multiple outbreaks could be generated within a single simulation run, and to minimize the influence of the initial conditions on the timing of those outbreaks.
While it would be possible to simulate shorter, and more, time series, producing the same total number of years, this would likely decrease the variability in when outbreaks were initiated within a time series, reducing the number of opportunities for overlap between rubella and measles outbreaks.

With respect to simulating a daily timescale, we chose to use the Tau-leaping algorithm as it is a good approximation of the Gillespie stochastic simulation algorithm, while running at substantially higher speeds, allowing for a greater number of years to be simulated.
A daily resolution allows for a greater exploration of test result delays than an alternative model formulation such as the TSIR model, as well as simplifying the integration of two disease simulations with different latent and infectious periods.
Daily simulations also provide greater flexibility for future extensions that may explore the time-varying nature of test sensitivity, as an individual progresses through their infectious period (see cited work by Middleton and Larremore, 2024, for example).


== Reviewer \#2

=== Overview

#reviewer[
  This paper presents a modelling study that evaluates how accurately measles outbreaks can be detected in syndromic surveillance time-series data, over a range of scenarios that consider a variety of diagnostic tests, testing rates, and epidemiological contexts.
  As the authors highlight in the introduction, rapid diagnostics tests for measles are being developed and evaluated for surveillance purposes.
  Here, the authors identify scenarios where using rapid diagnostic tests with imperfect sensitivity and specificity can yield similar outbreak detection performance to that achieved by using a perfect test.
  They also quantify how the structure and magnitude of background noise in the surveillance time-series data affects the outbreak detection performance for each type of diagnostic test.

  The manuscript is well written, the modelling framework and analyses are described in detail, and the findings are both sensible and interesting.
  I have provided a number of comments below, primarily regarding clarifications and seeking further details, and concerning the chosen measures of outbreak detection performance.
]

=== Major Comments

==== Comment 2.1 (Alert Definition and Negative Delays) <comment-2-1>

#reviewer[
  The conditions that trigger an alert are defined in the Methods section ("Triggering Alerts", p17):

  #quote[We define an "alert" as any consecutive string of 1 or more days where the 7-day (trailing) moving average of the test positive cases is greater than or equal to a pre-specified alert threshold, T.]

  And because an alert period may begin before an outbreak and span part of the outbreak, early alarms are considered to be correct with a negative delay:

  #quote[If the alert period starts before the outbreak and continues past the start date of the outbreak, this would be considered a correct alert with a negative delay i.e., an early warning triggered by false positive test results.]

  I appreciate why the authors chose to define alerts in this way.
  However, I think there's merit in also presenting a brief analysis that distinguishes between alerts with negative delays (i.e., false alarms that persist long enough to become considered correct) and alerts that are only triggered once the surveillance data includes evidence of an outbreak.
  As the authors emphasise (lines 158-161), maintaining a perpetual alert status is not feasible.
]

Thank you for your comment #sym.dash.en as a result of this we re-examined our assumptions and decided to adjust our algorithm that matches alerts to outbreaks.
In our reformulated matching algorithm, we only allow one outbreak to be matched to each alert.
Previously this was not the case, and a single alert could be successfully matched to multiple outbreaks if it (the alert) is long enough.
We have also implemented a new filtering algorithm to the alert threshold optimization to allow for only positive-delay matches to be successfully counted within the optimization phase, although the main results do not have this restriction and this is only an option for the user to implement.
Upon re-running the analysis with the new single outbreak per alert matching strategy, all negative delays disappeared.
We think that this matching is more intuitive for the reasons highlighted above: a perpetual alert status is not feasible and would certainly be investigated before multiple outbreaks start, finish, and restart.

==== Comment 2.2 (Framework Flexibility and Trade-offs)

#reviewer[
  Regarding trade-offs between true and false alert rates, the authors state in the Discussion section (lines 168-171):

  #quote[These trade-offs must be explicitly acknowledged when designing surveillance systems, and we present a framework to account for the deep interconnectedness of individual and population-level uncertainties that arise from necessary categorizations.]

  How flexible is this framework in allowing a potential user to define their preferred trade-off(s) and identify optimal detection thresholds accordingly? It seems to me that the methods as described should apply to any definition for surveillance accuracy, but I don't see any statements to this effect in the manuscript (they could accompany the text in lines 214-222, for example).
]

This framework is very flexible, however, at present the code does not allow the user to supply a different objective function.
Since the original submission we have substantially reworked the package code to make it more modular and clearer to navigate.
As such, redefining the objective function should not be necessary as the surveillance accuracy function can be overloaded to accept a new Sum Type that creates a new balance between the sensitivity and PPV metrics that the current objective function generate.
#link(
  "https://github.com/arnold-c/OutbreakDetection/blob/f6034987dde1c6cafb07a3e14358dd289c69f87b/OutbreakDetectionCore/src/detection/detection-metric-functions.jl#L28-L34",
)[These lines of the package code illustrate what would need to be defined in the end-user's code.
  We have added this detail to lines 249-254.
]

#reviewer[
  On a related note, when using imperfect tests in the present of high levels of dynamical noise, the authors write (lines 212-214):

  #quote[In these situations, the added complexity from large numbers of false positive test results likely warrants a different decision criteria than a binary detection threshold.]

  Are there relevant real-world examples of decision criteria that the authors could suggest here?
]

Médecins Sans Frontières (MSF) routinely responds to outbreaks in the Democratic Republic of the Congo (DRC).
This work also involves participating in and conducting sustained surveillance activities in partnership with the DRC Ministry of Health.
Because test coverage with ELISA tests is low relative to measles incidence, MSF often relies on clinical suspicion (which could be viewed as the close to the equivalent of a diagnostic test with 100% sensitivity and 0% specificity).
Depending on the source of the "test" results, in combination with the estimated vaccination coverage and recency of a Supplemental Immunization Activity (SIA) in the region, the thresholds used to determine a "suspected outbreak" vary from #emph["5 suspected cases reported by a single geographic unit in a one month period"] to #emph["Within a geographic unit: Number of cases or weekly incidence higher than in previous (non-epidemic) years, or the same as in an epidemic year OR If no data from previous years: increase in the number of cases in the last 3 or 4 weeks"] to #emph["> 2 confirmed cases (IgM+) in a one-month period"] (#link("https://medicalguidelines.msf.org/en/viewport/mme/english/3-3-confirming-the-outbreak-32407880.html")[MSF decision criteria here).
  We have added a summary of this to lines 237-241.
]

==== Comment 2.3 (Results Section Details)

#reviewer[
  With the methods section coming at the end of the paper, I think some additional details are needed at the start of the results section to make it easier for the reader to correctly interpret and evaluate the presented results.

  (a) Lines 84-85 report the range of optimal alert thresholds as "between 0.39 and 16.80 test positive cases per day" without specifying a time period.
  The abstract refers to a "7-day rolling average", and this is mentioned two paragraphs later in the results (line 99).
  But there is no mention of a 7-day window in the main text until after lines 84-85.
]

Added.

#reviewer[
  (b) Optimal alert thresholds were defined as maximising surveillance accuracy, which was defined as "the arithmetic mean of the system's PPV and sensitivity" (Methods, lines 321-322).
  I would like to see this definition included in the first paragraph of the results section, so that the reader can interpret the results without having to jump to the end of the manuscript.
]

Added.

==== Comment 2.4 (Figure 1 Interpretation)

#reviewer[
  Regarding Figure 1, the sharp drops and rises in the alert threshold that are shown in the bottom-right panel (Dynamical Noise, $Lambda$(8)) are due to the interplay between false-positives and high testing rates, as the authors describe in the text (lines 98-101).
  The subsequent rebound in alert threshold for the imperfect test (85%) when 90% of individuals are tested is not explained, and it's not obvious to me whether the false-positives and/or false-negatives could be the cause.
  Can the authors offer any insights?
]


In this revised version of the analysis this artifact is no longer present.
In the prior version, the rebound appears to occur due to a combination of simulation choices that when combined, result in a flat multi-modal minima.
One is the choice in allowing multiple outbreaks to be "correctly" matched to a single alert, not penalizing very few but long alert durations, as discussed above.

The others relate to how we handled the calculation of test positive individuals.
Specifically, for computational efficiency, we originally calculated both the number of individuals tested, and the number of test positives, by multiplying the integer values by the respective proportions (test sensitivity/specificity for test positives) before rounding to convert back to an integer value.
This has been updated to use binomial sampling #sym.dash.em extensive refactoring substantially increased the simulation speed, allowing for this change.
Previously, we also rounded the moving average of test positive values to generate integer values.
We no longer perform this rounding.
Both of these prior decisions resulted in flatter minima for the optimizer as they force alert thresholds close to each other (e.g., 8.4 and 7.8 alerts) to produce more similar results than they otherwise would.
Details of the binomial sampling have been added to lines 318-319 and line 334-337.


==== Comment 2.5 (Test Result Lags Clarity)

#reviewer[
  I found it difficult to understand the remarks about test result lags decreasing the surveillance accuracy (lines 102-105) and I suspect it's due to a single choice of word.
  When the authors state "This will occur if an alert occurs within the duration of the lag (e.g., 14 days) of the end of the outbreak", am I correct in thinking that refers to an alert occurring 1 to 14 days *after* the end of the outbreak? I initially interpreted "within the duration of the lag (e.g., 14 days) *of* the end of the outbreak" as referring to a window of 14 days in either direction (before or after the end of the outbreak).
]

Thank you for highlighting this.
We have added the word #emph[after], as well further details on lines 113-116 to clarify this process.

==== Comment 2.6 (Seasonality Assumptions)

#reviewer[
  The seasonality for the rubella noise was simulated to be in-phase with measles (lines 257-258).
  This is an interesting choice, presumably then the background signal (noise) is proportional to the outbreak signal? I suspect this makes detection easier than if the seasonality was chosen at random, which might lead to many more false positives?

  I had a look at the code and despite not being particularly familiar with Julia, I found it well structured and commented.
  I noticed that `dynamical_noise_correlation` could be either `"in-phase"`, `"out-of-phase"`, or `"none"`, but that the second and third options were commented out in `scripts/ensemble-sim.jl`.
  I'm curious why the authors chose not to also explore these scenarios and report the results in this manuscript — but please note that I am *not* suggesting that they should do so in this manuscript.
]

Thank you for examining the code and for this question.
While developing this project we did explore the effect of seasonality on the results presented.
Ultimately, there was a relatively minor impact on the quantitative results, and no difference on the qualitative conclusions, so we decided not to include those additional analyses in the manuscript for brevity and clarity #sym.dash.em the main results already required a substantial amount of explanation to build up to, and seasonality was secondary to this.


=== Minor Comments

==== Comment 2.7 (Nowcasting Approach)

#reviewer[
  In the Discussion section the authors mention a more computationally expensive approach based on now-casting case incidence (lines 184-186).
  These nowcasts would presumably aim to characterise the dynamics noise process, as the background against which an outbreak detection algorithm would look for signals of an outbreak, right?
]


Extending this work to incorporate nowcasting could be as simple as estimating the reporting rates and correcting the observed test positive numbers (without making corrections for noise-induced false positives), before applying the alert threshold.
Equally, it could be much more sophisticated and try to characterize the background dynamical noise process, as you suggest, and distinguish it from the target disease.
If successful, this second approach would potentially lead to a far more accurate alert system with imperfect tests, but it becomes less necessary as test sensitivity and specificity improves.


==== Comment 2.8 (Tau-leaping Algorithm)

#reviewer[
  At the very start of the Methods section (lines 232-234), the authors mention that they used "a modified Tau-leaping algorithm", but don't explain what modifications were made.
]


Thank you for pointing this out.
We have updated lines 266-267 to clarify that we use binomial draws to ensure small compartments remain positive valued, as using draws from a Poisson distribution can result in negative compartments if the original size is small (as is the case for the $E$ and $I$ compartments).


==== Comment 2.9 (Importation Rate)

#reviewer[
  Regarding the importation rate defined in Table 1 (page 14), what does the mu parameter represent? And the importation rate defined here doesn't appear to be proportional to the population size N (as stated in the text above).
]


#sym.mu is the birth/death rate.
We have update Table 1 to reflect this #sym.dash.em thank you for catching this typo.
We have update the description of imports in lines 273-275 to clarify that it is proportional to the inverse of $N$.


==== Comment 2.10 (Testing Selection Process)

#reviewer[
  For the given percentage P of cases that are tested (lines 286-288), are individuals selected for testing using a binomial process, or is a fixed fraction of all individuals selected?
]

As part of the refactoring described above, we have changed this to move from using a fixed fraction to binomial draws for both the number of individuals tested, as well as the number of test positive results produced.

==== Comment 2.11 (Test Result Delay)

#reviewer[
  A perfect test with a 14-day test result delay represents a best-case test under "more realistic reporting delays in result return".
  What are the major contributing factors to such a long delay? Is it primarily the delay between symptom onset and production of antibodies to a detectable level?
]


In general, the turnaround time for the ELISA from specimen collection to central lab to centralized report is a substantial component of the reporting lag, and in many regions of the Democratic Republic of the Congo, for example, could be as long as 3+ weeks.
As such, we chose to focus on this delay when defining test result lag.
However, you raise an excellent point regarding IgM detectability, which follows much of the same reasoning we detail in lines 219-232 of the discussion.
IgM detectability is at its peak at approximately 3+ days post rash onset (https://www.cdc.gov/mumps/media/pdfs/2025/02/MMRV-Testing-for-Clinicians_Jan2025.pdf), which occurs approximately 4 days into an individual's infectious period (that typically lasts about 8 days).
As such, in most individuals, they are only detectable via IgM for approximately 25% of their infectious period.
Future work should aim to incorporate this aspect into the analysis by either explicitly modeling viral kinetics, or by using an alternative approach to adjust the test sensitivity and specificity throughout an individual's infectious period.


== Logic Changes <logic-changes>

=== Calculation of moving averages

`./OutbreakDetectionCore/src/utilities/calculate-moving-average.jl`

The core calculations are the same, but instead of rounding, converting, and returning an Int64 integer for the moving average (of the test positives, principally), return a Float64 floating point value.

=== Parameter sampling

In the previous version, there were few parameters that were sampled in each of the simulations of the ensemble.
In the refactored version, most are sampled.

`./OutbreakDetectionCore/src/diagnostic-testing/calculate-num-tested.jl`

`./OutbreakDetectionCore/src/diagnostic-testing/calculate-num-positive.jl`

`./OutbreakDetectionCore/src/types/dynamics-parameters.jl`

Instead of multiplying the integer number of individuals by the appropriate proportion, rounding and then converting to an integer, we now sample from a Binomial distribution to calculate the number of individuals tested and test positive.
This was changed to sampling to reduce the impact of small numbers all being rounded to 0 equally.

=== Fixed noise vaccination coverage

`./OutbreakDetectionCore/src/threshold-optimization/evaluate-missing-optimizations.jl`

`./OutbreakDetectionCore/src/noise/noise-parameters-optimization.jl`

As part of the refactoring, the new implementation uses MultistartOptimization to optimize the vaccination level in the dynamical noise simulations to achieve $Lambda(c)$ noise level.
The benefit of the refactor is that if any other scenarios need to be run, or a change to the target or noise disease etc. then the optimization will handle those changes.
With the other changes detailed, the maximum dynamical noise scaling that can be achieved is $Lambda(7)$, not $Lambda(8)$ that was used in the original manuscript.

=== Endemic noise initialization

`./OutbreakDetectionCore/src/noise/noise-recreation.jl`

During the refactoring, a calculation of the endemic state for the dynamical noise was added to seed the noise simulations in an endemic state.
This was done to avoid always starting with a large outbreak and wasn't present in the original implementation.

=== SEIR model Poisson transitions

`./OutbreakDetectionCore/src/simulation/seir-model.jl`

In the refactored version, we use the `_smart_transition()` function to use Poisson samples when the size of the compartments is large relative to the rate, and Binomial when small.
This is done for computational efficiency, and makes a negligible difference to the number of individuals moving between compartments.

=== Alert-outbreak matching

`./OutbreakDetectionCore/src/detection/match-alert-outbreak-thresholds.jl`

As mentioned in a response to a #link(<comment-2-1>)[comment 2.1], in the newest analysis we have updated the alert-outbreak matching algorithm.
Previously, a single alert could be correctly matched to multiple outbreaks.
Now, the default matching pattern no longer allows this and only matches to the first outbreak it overlaps with.
As a result, infrequent and long alerts generated by very low alert thresholds are penalized as only the first outbreak is counted as detected.
Under the prior implementation, all outbreaks would be "detected", resulting in a higher sensitivity than would otherwise be expected, with no penalization to the PPV of the system.
In the extreme, a single long alert that spans the whole time series would match all outbreaks, and the PPV would be 100% as there is only one alert and it has been matched to at least one outbreak.

In one respect, this accounts for the reality that once an alert has triggered and been investigated, it is likely not going to result in additional follow up if there is not a cessation in alert status.

=== Frequency Dependent Beta value

`./OutbreakDetectionCore/src/types/dynamics-parameters.jl`

The `DynamicsParameterSpecification` constructor used to normalize the value of `beta_mean` by the initial population size `N`, rather than using a strictly frequency-dependent approach where the normalization is included in the simulation of infections.
This has been changed for consistency reasons.

== Bug Fixes <bug-fixes>

=== Noise dynamics parameters

`./OutbreakDetectionCore/src/noise/noise-dynamics-parameters.jl`

In the previous version, the noise dynamics parameters struct was creating using the target disease's `beta_mean` values i.e,. measles' `beta_mean`, not the value calculated for rubella.

=== Beta value calculation

`./OutbreakDetectionCore/src/simulation/transmission-functions.jl`

The previous version calculated `beta` using the SIR equation, not SEIR equation that includes the movement out of the latent state.


=== SEIR model beta value index

`./OutbreakDetectionCore/src/simulation/seir-model.jl`

In the refactored version, the SEIR model loop uses the previous value of the beta to determine how many individuals are infected at the current time step (also using the previous values of the states).
This is correct as the otherwise the individuals are being affected by the state of a dynamical system that doesn't (quite) exist yet.

=== SEIR model imported infections

`./OutbreakDetectionCore/src/simulation/seir-model.jl`

In the refactored version, we calculate the number of importer individuals as a binomial sample from the number of remaining susceptible individuals (after removing those from direct contact infections).
The previous version sampled the number of imported individuals from the total population size.
