#import "template.typ": article
//
#show: article.with(
        title: "Individual and Population Level Uncertainty Interact to Determine Performance of Outbreak Detection",
        header-title: "true",
        authors: (
	  "Callum R.K. Arnold": (
            affiliation: ("PSU-Bio", "CIDD"),
            corresponding: "true",
            email: "contact\@callumarnold.com",
	  ),
	  "Alex C. Kong": (
            affiliation: ("Hopkins-IH"),
	  ),
	  "Amy K. Winter": (
            affiliation: ("UGA"),
	  ),
	  "William J. Moss": (
            affiliation: ("Hopkins-IH", "Hopkins-Epi"),
	  ),
	 "Bryan N. Patenaude": (
            affiliation: ("Hopkins-IH"),
	  ),
	  "Matthew J. Ferrari": (
            affiliation: ("PSU-Bio", "CIDD"),
	  ),
	),
        affiliations: (
          "PSU-Bio": "Department of Biology, Pennsylvania State University, University Park, PA, USA 16802",
          "CIDD": "Center for Infectious Disease Dynamics, Pennsylvania State University, University Park, PA, USA 16802",
          "Hopkins-IH": "Department of International Health, Johns Hopkins Bloomberg School of Public Health, Baltimore, MD, USA 21205",
          "UGA": "Department of Epidemiology, College of Public Health, University of Georgia, Athens, GA, USA 30602",
          "Hopkins-Epi": "Department of Epidemiology, Johns Hopkins Bloomberg School of Public Health, Baltimore, MD, USA 21205"
	),
)

#show figure.where(kind: image): set figure(supplement: "Supplemental Figure")
#show figure.where(kind: table): set figure(supplement: "Supplemental Table")
#set math.equation(
  numbering: "1",
  supplement: "Supplemental Equation"
)

== Tables

#let accuracy = csv("supplemental_files/tables/optimal-thresholds_accuracy.csv")
#figure(
    table(
    columns: 9,
    fill: (x, y) => {
      if y == 0 {gray}
      if y == 1 {gray}
    },
    align: center,
    [], table.cell(colspan: 2, align: center, "Test Characteristic"), table.cell(colspan: 6, align: center, "Testing Rate"),
    ..accuracy.flatten()
  ),
  caption: [Mean outbreak detection accuracy for imperfect and perfect diagnostic tests, at their respective optimal alert thresholds, under dynamical and Poisson noise structures where the average daily noise incidence is 8 times the average daily measles incidence $Lambda(8)$. The test sensitivity equals the test specificity for all diagnostic tests.]
)
<tbl-optimal-thresholds-accuracy>

#let delays = csv("supplemental_files/tables/optimal-thresholds_detection-delays.csv")
#figure(
    table(
    columns: 9,
    fill: (x, y) => {
      if y == 0 {gray}
      if y == 1 {gray}
    },
    align: center,
    [], table.cell(colspan: 2, align: center, "Test Characteristic"), table.cell(colspan: 6, align: center, "Testing Rate"),
    ..delays.flatten()
  ),
  caption: [Mean outbreak detection delays (days) for imperfect and perfect diagnostic tests, at their respective optimal alert thresholds, under dynamical and Poisson noise structures where the average daily noise incidence is 8 times the average daily measles incidence $Lambda(8)$. The test sensitivity equals the test specificity for all diagnostic tests.]
)
<tbl-optimal-thresholds-delays>

#let unavoidable = csv("supplemental_files/tables/optimal-thresholds_unavoidable-cases.csv")
#figure(
    table(
    columns: 9,
    fill: (x, y) => {
      if y == 0 {gray}
      if y == 1 {gray}
    },
    align: center,
    [], table.cell(colspan: 2, align: center, "Test Characteristic"), table.cell(colspan: 6, align: center, "Testing Rate"),
    ..unavoidable.flatten()
  ),
  caption: [Mean unavoidable cases per annum (scaled to Ghana's 2022 population) for imperfect and perfect diagnostic tests, at their respective optimal alert thresholds, under dynamical and Poisson noise structures where the average daily noise incidence is 8 times the average daily measles incidence $Lambda(8)$. The test sensitivity equals the test specificity for all diagnostic tests.]
)
<tbl-optimal-thresholds-unavoidable>


== Figures

#figure(
  image("supplemental_files/plots/optimal-thresholds_alert-duration-plot.svg"),
  caption: [The duration of alerts of outbreak detection systems under different testing rates and noise structures, at their respective optimal alert thresholds. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Imperfect tests have the same values for sensitivity and specificity. Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays. $Lambda(4)$ indicates the mean noise incidence is 4 times higher than the mean measles incidence, for example.]
)
<fig-alert-duration>

#figure(
  image("supplemental_files/plots/optimal-thresholds_n-alerts-plot.svg"),
  caption: [The number of alerts of outbreak detection systems under different testing rates and noise structures, at their respective optimal alert thresholds. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Imperfect tests have the same values for sensitivity and specificity. Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays. $Lambda(4)$ indicates the mean noise incidence is 4 times higher than the mean measles incidence, for example.]
)
<fig-num-alerts>
