#import "template.typ": article
#import "metadata.typ": meta

#show: article.with(..meta)

#show figure.where(kind: image): set figure(supplement: "Supplemental Figure")
#show figure.where(kind: table): set figure(supplement: "Supplemental Table")
#set math.equation(numbering: "1", supplement: "Supplemental Equation")

== Tables

#let optimal_thresholds = csv(
  "supplemental_files/tables/optimal-thresholds.csv",
)
#figure(
  table(
    columns: 13,
    fill: (x, y) => {
      if y == 0 { gray }
      if y == 1 { gray }
    },
    align: center,
    [], table.cell(
      colspan: 2,
      align: center,
      "Test Characteristic",
    ), table.cell(colspan: 10, align: center, "Testing Rate"),
    ..optimal_thresholds.flatten()
  ),
  caption: [The optimal outbreak alert thresholds for imperfect and perfect diagnostic tests (that maximizes outbreak detection accuracy) under dynamical and static noise structures where the average daily noise incidence is 8 times the average daily measles incidence $Lambda(8)$. The test sensitivity equals the test specificity for all diagnostic tests.],
)
<tbl_od-optimal-thresholds>



// #let accuracy = csv("supplemental_files/tables/optimal-thresholds_accuracy.csv")
// #figure(
//     table(
//     columns: 13,
//     fill: (x, y) => {
//       if y == 0 {gray}
//       if y == 1 {gray}
//     },
//     align: center,
//     [], table.cell(colspan: 2, align: center, "Test Characteristic"), table.cell(colspan: 10, align: center, "Testing Rate"),
//     ..accuracy.flatten()
//   ),
//   caption: [Mean outbreak detection accuracy for imperfect and perfect diagnostic tests, at their respective optimal alert thresholds, under dynamical and Poisson noise structures where the average daily noise incidence is 8 times the average daily measles incidence $Lambda(8)$. The test sensitivity equals the test specificity for all diagnostic tests.]
// )
// <tbl-optimal-thresholds-accuracy>

// #let delays = csv("supplemental_files/tables/optimal-thresholds_detection-delays.csv")
// #figure(
//     table(
//     columns: 13,
//     fill: (x, y) => {
//       if y == 0 {gray}
//       if y == 1 {gray}
//     },
//     align: center,
//     [], table.cell(colspan: 2, align: center, "Test Characteristic"), table.cell(colspan: 10, align: center, "Testing Rate"),
//     ..delays.flatten()
//   ),
//   caption: [Mean outbreak detection delays (days) for imperfect and perfect diagnostic tests, at their respective optimal alert thresholds, under dynamical and Poisson noise structures where the average daily noise incidence is 8 times the average daily measles incidence $Lambda(8)$. The test sensitivity equals the test specificity for all diagnostic tests.]
// )
// <tbl-optimal-thresholds-delays>
//
// #let unavoidable = csv("supplemental_files/tables/optimal-thresholds_unavoidable-cases.csv")
// #figure(
//     table(
//     columns: 9,
//     fill: (x, y) => {
//       if y == 0 {gray}
//       if y == 1 {gray}
//     },
//     align: center,
//     [], table.cell(colspan: 2, align: center, "Test Characteristic"), table.cell(colspan: 6, align: center, "Testing Rate"),
//     ..unavoidable.flatten()
//   ),
//   caption: [Mean unavoidable cases per annum (scaled to Ghana's 2022 population) for imperfect and perfect diagnostic tests, at their respective optimal alert thresholds, under dynamical and Poisson noise structures where the average daily noise incidence is 8 times the average daily measles incidence $Lambda(8)$. The test sensitivity equals the test specificity for all diagnostic tests.]
// )
// <tbl-optimal-thresholds-unavoidable>
//

== Figures

=== Arithmetic Mean Optimizations

#figure(
  image("supplemental_files/plots/optimal-thresholds_n-alerts-plot.svg"),
  caption: [The number of alerts of outbreak detection systems under different testing rates and noise structures, at their respective optimal alert thresholds. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Imperfect tests have the same values for sensitivity and specificity. Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays. $Lambda(4)$ indicates the mean noise incidence is 4 times higher than the mean measles incidence, for example.],
)
<fig-num-alerts>

#figure(
  image("supplemental_files/plots/optimal-thresholds_alert-duration-plot.svg"),
  caption: [The duration of alerts of outbreak detection systems under different testing rates and noise structures, at their respective optimal alert thresholds. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Imperfect tests have the same values for sensitivity and specificity. Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays. $Lambda(4)$ indicates the mean noise incidence is 4 times higher than the mean measles incidence, for example.],
)
<fig-alert-duration>


#figure(
  image("supplemental_files/plots/optimal-thresholds_unavoidable-plot.svg"),
  caption: [The number of unavoidable cases of an outbreak detection systems under different testing rates and noise structures, at their respective optimal alert thresholds. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Imperfect tests have the same values for sensitivity and specificity. Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays. $Lambda(4)$ indicates the mean noise incidence is 4 times higher than the mean measles incidence.],
)
<fig-unavoidable>

=== F-1 Score Optimizations

#figure(
  image(
    "supplemental_files/plots/f1_optimal-thresholds_alert-threshold-plot.svg",
  ),
  caption: [The optimal alert threshold of outbreak detection systems (that maximizes outbreak detection F-1 score) under different testing rates and noise structures. Imperfect tests have the same values for sensitivity and specificity. Circular markers represent tests with 0-day turnaround times, and triangular markers represent tests with result delays. $Lambda(4)$ indicates the mean noise incidence is 4 times higher than the mean measles incidence, for example.],
)
<fig-f1-alert-threshold>

#figure(
  image("supplemental_files/plots/f1_optimal-thresholds_accuracy-plot.svg"),
  caption: [The accuracy of outbreak detection systems under different testing rates and noise structures, at their respective (F-1 score) optimal alert thresholds. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Imperfect tests have the same values for sensitivity and specificity. Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays. $Lambda(4)$ indicates the mean noise incidence is 4 times higher than the mean measles incidence, for example.],
)
<fig-f1-accuracy>

#figure(
  image("supplemental_files/plots/f1_optimal-thresholds_delays-plot.svg"),
  caption: [The detection delay of outbreak detection systems under different testing rates and noise structures, at their respective (F-1 score) optimal alert thresholds. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Imperfect tests have the same values for sensitivity and specificity. Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays. $Lambda(4)$ indicates the mean noise incidence is 4 times higher than the mean measles incidence.],
)
<fig-f1-delay>

#figure(
  image("supplemental_files/plots/f1_optimal-thresholds_prop-alert-plot.svg"),
  caption: [The proportion of the time series in alert of outbreak detection systems under different testing rates and noise structures, at their respective (F-1 score) optimal alert thresholds. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Imperfect tests have the same values for sensitivity and specificity. Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays. $Lambda(4)$ indicates the mean noise incidence is 4 times higher than the mean measles incidence.],
)
<fig-f1-alert-proportion>
