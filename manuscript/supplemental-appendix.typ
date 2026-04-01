#import "template.typ": article
#import "metadata.typ": meta

#show: article.with(..meta)

#show figure.where(kind: image): set figure(supplement: "Supplemental Figure")
#show figure.where(kind: table): set figure(supplement: "Supplemental Table")
#set math.equation(numbering: "1", supplement: "Supplemental Equation")

== Tables


== Figures

=== Arithmetic Mean Optimizations


=== F-1 Score Optimizations


== Multiple Outbreaks per Alert Optimizations

In the main manuscript, the matching process between alerts and outbreaks requires an alert to only be "attached" to a single outbreak.
This means that a single long alert resulting from a low alert threshold will only be deemed to have detected one outbreak, even if its duration is long enough to span multiple outbreaks.
In effect, this penalizes more sensitive alert systems that may prioritize lower alert thresholds but would produce fewer alerts, therefore resulting in a larger number of "missed outbreaks".
Below, we show the results of an alternative matching scheme where this restriction is not applied i.e., a single long alert can be ascribed to multiple outbreaks in the calculation of alert accuracy.

The general patterns follow the main results, where increasing static noise has relatively limited influence on outbreak detection accuracy for perfect and imperfect tests, and decreasing detection accuracy for imperfect tests with increasing levels of dynamic noise.
Unlike the main results, under this alternate matching strategy, imperfect tests did not show a monotonic increase in the optimal alert threshold as the proportion of individuals tested increased (@fig-multi-outbreak-matching-alert-threshold).
Instead, under both static and dynamical noise scenarios, as testing and noise magnitude increased, eventually the number of false positives overwhelmed the number of true positives being detected, resulting in a system that produced the "optimal" alert accuracy by reducing the threshold to, or close to, 0 (@fig-multi-outbreak-matching-alert-threshold).
Consequently, a single (or low number of) long alert period(s) "correctly" capture all outbreaks, resulting in high sensitivity while being minimally penalized as those alerts were long enough to eventually overlap with an outbreak such that they are never false alerts (@fig-multi-outbreak-matching-accuracy).
This can be observed by very large #emph[negative] detection delays (i.e., the alert predated the outbreak's initiation) and the high proportion of the time series being under alert status (@fig-multi-outbreak-matching-delay, and @fig-multi-outbreak-matching-alert-proportion, respectively).

// However, unlike the main results, at high levels of dynamical noise ($gt.eq Lambda (6)$), high testing rates ($gt.eq$ 50%) result in a marked #emph[increase] in outbreak detection accuracy.
// In these scenarios the optimal outbreak alert threshold falls to 0.39 (daily) test positive cases in a 7-day moving average: 3 positive test cases in a week would be sufficient to trigger an outbreak alert (@fig-multi-outbreak-matching-alert-threshold, @fig-multi-outbreak-matching-accuracy).
//

#figure(
  image(
    "supplemental_files/plots/multi-outbreak-matching_optimal-thresholds_alert-threshold-plot.svg",
  ),
  caption: [The optimal alert threshold of outbreak detection systems (that maximizes outbreak detection arithmetic mean) under different testing rates and noise structures, where a single alert can be matched to multiple outbreaks. Each imperfect test uses the same value for both its sensitivity and specificity (either 85% or 90%). Circular markers represent tests with 0-day turnaround times, and triangular markers represent tests with result delays. $Lambda(4)$ indicates the mean noise incidence is 4 times higher than the mean measles incidence, for example.],
)
<fig-multi-outbreak-matching-alert-threshold>

#figure(
  image(
    "supplemental_files/plots/multi-outbreak-matching_optimal-thresholds_accuracy-plot.svg",
  ),
  caption: [The accuracy of outbreak detection systems under different testing rates and noise structures, at their respective (arithmetic mean) optimal alert thresholds, where a single alert can be matched to multiple outbreaks. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Each imperfect test uses the same value for both its sensitivity and specificity (either 85% or 90%). Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays. $Lambda(4)$ indicates the mean noise incidence is 4 times higher than the mean measles incidence, for example.],
)
<fig-multi-outbreak-matching-accuracy>

#figure(
  image(
    "supplemental_files/plots/multi-outbreak-matching_optimal-thresholds_delays-plot.svg",
  ),
  caption: [The detection delay of outbreak detection systems under different testing rates and noise structures, at their respective (arithmetic mean) optimal alert thresholds, where a single alert can be matched to multiple outbreaks. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Each imperfect test uses the same value for both its sensitivity and specificity (either 85% or 90%). Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays. $Lambda(4)$ indicates the mean noise incidence is 4 times higher than the mean measles incidence.],
)
<fig-multi-outbreak-matching-delay>

#figure(
  image(
    "supplemental_files/plots/multi-outbreak-matching_optimal-thresholds_prop-alert-plot.svg",
  ),
  caption: [The proportion of the time series in alert of outbreak detection systems under different testing rates and noise structures, at their respective (arithmetic mean) optimal alert thresholds, where a single alert can be matched to multiple outbreaks. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Each imperfect test uses the same value for both its sensitivity and specificity (either 85% or 90%). Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays. $Lambda(4)$ indicates the mean noise incidence is 4 times higher than the mean measles incidence.],
)
<fig-multi-outbreak-matching-alert-proportion>

