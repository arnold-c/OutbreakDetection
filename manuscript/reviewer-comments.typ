#show heading.where(level: 1): it => align(center, it)

// Function to format the planned response to reviewer comments
#let responseplan(content) = {
  block(
    inset: (left: 2em),
    text(fill: luma(45%), emph(content)),
  )
}

// Function to format responses to reviewer comments
#let response(content) = {
  block(
    inset: (left: 2em),
    text(fill: red, emph(content)),
  )
}

= Reviewer Comments

== Journal Requirements

1. We ask that a manuscript source file is provided at Revision. Please upload your manuscript file as a .doc, .docx, .rtf or .tex. If you are providing a .tex file, please upload it under the item type 'LaTeX Source File' and leave your .pdf version as the item type 'Manuscript'.

2. Please upload all main figures as separate Figure files in .tif or .eps format. For more information about how to convert and format your figure files please see our guidelines: https://journals.plos.org/ploscompbiol/s/figures

3. We have noticed that you have uploaded Supporting Information files, but you have not included a list of legends. Please add a full list of legends for your Supporting Information files after the references list.

4. Please revise your current Competing Interest statement to the standard "The authors have declared that no competing interests exist."

_Note: If the reviewer comments include a recommendation to cite specific previously published works, please review and evaluate these publications to determine whether they are relevant and should be cited. There is no requirement to cite these works unless the editor has indicated otherwise._

== Reviewer \#1

Dear Authors,

Thank you for the opportunity in reviewing this manuscript. The work is a consistent and detailed model of detectability of outbreak given different accuracy of tests, which is a crucial part in outbreak detection.

Though interesting and worth publishing, I have minor comments that should be addressed.

=== Comment 1.1 (Line 100, Figure 1)

Figure 1 shows that perfect test with 0-day and 14-days delay performs exactly equal, which is not wrong as the figure only shows proportion of individual tested against alert threshold. It would be interesting to add or to point from that figure to a time series plot of each of the 4 different scenarios showing how the delay in testing would affect the timing of alert, which is of importance in outbreak detection algorithm. As supplementary material is okay.

=== Comment 1.2 (Figure 2 and Line 102)

The same figure 2 and for the text that follows it, around line 102. Would be great to have another version of Figure 3 where is clear shown that perfect tests and imperfect test separately. The difference in perfect test with different delay do not yield any clear or statistical significant difference of outbreak alert, not sure if I followed the explanation given to this. Though is clear from figure 3 that imperfect testing can be problematic triggering the surveillance system in unnecessary levels of alert and by chance or seasonality of the disease, such as is for measles.

=== Comment 1.3 (Line 173, Citations)

The phrase is not completely correct, Larremore et al. 2021 that you cites is an example of modeling characteristics and its consequence for a wider surveillance program, as also Hay et al. 2021. In a more general sense and targeted to non modeling expert, the perspective work of Mina et al. 2020, entitled "Rethinking Covid-19 Test Sensitivity — A Strategy for Containment" addresses this, would recommend to reshape the phrase and cite those works.

=== Comment 1.4 (Simulation Duration)

Minor comment, Why the simulated time series is as long as 100 years? In daily resolution? This is extremely long and somewhat unrealistic. As measles settings around the world varies a lot would be interesting to address some short-term outbreaks and low level of incidence, as some settings are close to elimination. An interesting development would be to understand how effectively different this 4 kind of testing scenarios can affect the elimination of measles. No action needed and not mandatory to be addressed, though would interesting to read some speculation on those scenarios in light of the work developed at the present manuscript. Regarding this, what are the time unit for figure 5?

#responseplan[More about producing something with multiple outbreaks to remove influence of initial conditions - cuold have just simulated more smaller outbreaks.]

== Reviewer \#2

=== Overview

This paper presents a modelling study that evaluates how accurately measles outbreaks can be detected in syndromic surveillance time-series data, over a range of scenarios that consider a variety of diagnostic tests, testing rates, and epidemiological contexts. As the authors highlight in the introduction, rapid diagnostics tests for measles are being developed and evaluated for surveillance purposes. Here, the authors identify scenarios where using rapid diagnostic tests with imperfect sensitivity and specificity can yield similar outbreak detection performance to that achieved by using a perfect test. They also quantify how the structure and magnitude of background noise in the surveillance time-series data affects the outbreak detection performance for each type of diagnostic test.

The manuscript is well written, the modelling framework and analyses are described in detail, and the findings are both sensible and interesting. I have provided a number of comments below, primarily regarding clarifications and seeking further details, and concerning the chosen measures of outbreak detection performance.

=== Major Comments

==== Comment 2.1 (Alert Definition and Negative Delays)

The conditions that trigger an alert are defined in the Methods section ("Triggering Alerts", p17):

#quote[We define an "alert" as any consecutive string of 1 or more days where the 7-day (trailing) moving average of the test positive cases is greater than or equal to a pre-specified alert threshold, T.]

And because an alert period may begin before an outbreak and span part of the outbreak, early alarms are considered to be correct with a negative delay:

#quote[If the alert period starts before the outbreak and continues past the start date of the outbreak, this would be considered a correct alert with a negative delay i.e., an early warning triggered by false positive test results.]

I appreciate why the authors chose to define alerts in this way. However, I think there's merit in also presenting a brief analysis that distinguishes between alerts with negative delays (i.e., false alarms that persist long enough to become considered correct) and alerts that are only triggered once the surveillance data includes evidence of an outbreak. As the authors emphasise (lines 158-161), maintaining a perpetual alert status is not feasible.

#responseplan[What's the performance loss of retrospectively removing negative delays from true positives? Also think about how to explain using a different objective function and the changes to the results]

==== Comment 2.2 (Framework Flexibility and Trade-offs)

Regarding trade-offs between true and false alert rates, the authors state in the Discussion section (lines 168-171):

#quote[These trade-offs must be explicitly acknowledged when designing surveillance systems, and we present a framework to account for the deep interconnectedness of individual and population-level uncertainties that arise from necessary categorizations.]

How flexible is this framework in allowing a potential user to define their preferred trade-off(s) and identify optimal detection thresholds accordingly? It seems to me that the methods as described should apply to any definition for surveillance accuracy, but I don't see any statements to this effect in the manuscript (they could accompany the text in lines 214-222, for example).

#responseplan[Very flexible - just provide another measure of accuracy. Add comment in lines 214-222.]

#responseplan[Within reason, can get a different answer by changing the objective function. Framework done to illustrate this. We had to define some objective, and we chose to determine negative delays as correct. Objective to be pre-defined before doing the analysis - don't want to do post-hoc alterations to the objective to get a particular result.]


On a related note, when using imperfect tests in the present of high levels of dynamical noise, the authors write (lines 212-214):

#quote[In these situations, the added complexity from large numbers of false positive test results likely warrants a different decision criteria than a binary detection threshold.]

Are there relevant real-world examples of decision criteria that the authors could suggest here?

#responseplan[MSF multiple thresholds for measles outbreak response in DRC?]

==== Comment 2.3 (Results Section Details)

With the methods section coming at the end of the paper, I think some additional details are needed at the start of the results section to make it easier for the reader to correctly interpret and evaluate the presented results.

(a) Lines 84-85 report the range of optimal alert thresholds as "between 0.39 and 16.80 test positive cases per day" without specifying a time period. The abstract refers to a "7-day rolling average", and this is mentioned two paragraphs later in the results (line 99). But there is no mention of a 7-day window in the main text until after lines 84-85.

#responseplan[Add 7-day rolling average to lines 84-85]

(b) Optimal alert thresholds were defined as maximising surveillance accuracy, which was defined as "the arithmetic mean of the system's PPV and sensitivity" (Methods, lines 321-322). I would like to see this definition included in the first paragraph of the results section, so that the reader can interpret the results without having to jump to the end of the manuscript.

#responseplan[Add accuracy definition to line 81]

==== Comment 2.4 (Figure 1 Interpretation)

Regarding Figure 1, the sharp drops and rises in the alert threshold that are shown in the bottom-right panel (Dynamical Noise, Λ(8)) are due to the interplay between false-positives and high testing rates, as the authors describe in the text (lines 98-101). The subsequent rebound in alert threshold for the imperfect test (85%) when 90% of individuals are tested is not explained, and it's not obvious to me whether the false-positives and/or false-negatives could be the cause. Can the authors offer any insights?

#responseplan[I believe the rebound is because of an increase in false positive as a result of a larger number of individuals being tested]

==== Comment 2.5 (Test Result Lags Clarity)

I found it difficult to understand the remarks about test result lags decreasing the surveillance accuracy (lines 102-105) and I suspect it's due to a single choice of word. When the authors state "This will occur if an alert occurs within the duration of the lag (e.g., 14 days) of the end of the outbreak", am I correct in thinking that refers to an alert occurring 1 to 14 days *after* the end of the outbreak? I initially interpreted "within the duration of the lag (e.g., 14 days) *of* the end of the outbreak" as referring to a window of 14 days in either direction (before or after the end of the outbreak).

#responseplan[Correct - this is a typo to be corrected]

==== Comment 2.6 (Seasonality Assumptions)

The seasonality for the rubella noise was simulated to be in-phase with measles (lines 257-258). This is an interesting choice, presumably then the background signal (noise) is proportional to the outbreak signal? I suspect this makes detection easier than if the seasonality was chosen at random, which might lead to many more false positives?

I had a look at the code and despite not being particularly familiar with Julia, I found it well structured and commented. I noticed that `dynamical_noise_correlation` could be either `"in-phase"`, `"out-of-phase"`, or `"none"`, but that the second and third options were commented out in `scripts/ensemble-sim.jl`. I'm curious why the authors chose not to also explore these scenarios and report the results in this manuscript — but please note that I am *not* suggesting that they should do so in this manuscript.

#responseplan[Initially had simulated and explored this, but it was too complicated to describe within the scope of this paper given the number of other axes of analysis that need to be explained. Didn't really make a lot of difference to the results, particularly to moderately accurate tests (>=90% sens/spec), ultimately. At lower accuracy tests (80-90%), having no seasonality resulted in slightly higher accuracy outbreak detection]

=== Minor Comments

==== Comment 2.7 (Nowcasting Approach)

In the Discussion section the authors mention a more computationally expensive approach based on now-casting case incidence (lines 184-186). These nowcasts would presumably aim to characterise the dynamics noise process, as the background against which an outbreak detection algorithm would look for signals of an outbreak, right?

#responseplan[Nowcasting could be as simple as estimating the reporting rates and correcting the observed test positive numbers (with some noise), before applying the alert threshold to the estimated number of test positives, or could be more sophisticated and try to infer the background noise shape and magnitude. If successful in the latter approach, it would potentially provide more accurate alert systems.]

==== Comment 2.8 (Tau-leaping Algorithm)

At the very start of the Methods section (lines 232-234), the authors mention that they used "a modified Tau-leaping algorithm", but don't explain what modifications were made.

#responseplan[Update lines 232 to mention Binomial for small compartments that could go extinct (maybe update to use new approach in CSDNoise)]

==== Comment 2.9 (Importation Rate)

Regarding the importation rate defined in Table 1 (page 14), what does the mu parameter represent? And the importation rate defined here doesn't appear to be proportional to the population size N (as stated in the text above).

#responseplan[#sym.mu is the birth/death rate. Can update lines 243-244 to reflect that it's inversely proportional to $N$]

==== Comment 2.10 (Testing Selection Process)

For the given percentage P of cases that are tested (lines 286-288), are individuals selected for testing using a binomial process, or is a fixed fraction of all individuals selected?

#responseplan[Currently fixed, but should update this to use binomial process like CSDNoise, as well as the calculation of test positive/negative]

==== Comment 2.11 (Test Result Delay)

A perfect test with a 14-day test result delay represents a best-case test under "more realistic reporting delays in result return". What are the major contributing factors to such a long delay? Is it primarily the delay between symptom onset and production of antibodies to a detectable level?

#responseplan[More that the turnaround time for an ELISA to go from collection to central lab to centralized report. Should comment about mean time to symptom onset (X days? - incubation period) vs detectability (Y days?)]
