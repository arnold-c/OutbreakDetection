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
        // bib: "/Users/cfa5228/Documents/Repos/OutbreakDetection/manuscript/manuscript_files/OD.bib",
	// bib: "OD.bib",
)


= Results
== Tables

#let accuracy = csv("supplemental-appendix_files/tables/optimal-thresholds_accuracy.csv")
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
  caption: [Mean outbreak detection accuracy of each testing scenario at their specific optimal thresholds, when the average noise incidence is 8 times higher than the average measles incidence. A) the noise structure is dynamical, and the seasonality is in-phase with the measles incidence. B) the noise structure is Poisson only.]
)
<tbl-optimal-thresholds-accuracy>

#let unavoidable = csv("supplemental-appendix_files/tables/optimal-thresholds_unavoidable-cases.csv")
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
  caption: [Mean unavoidable cases per annum of each testing scenario at their specific optimal thresholds, scaled up to Ghana’s 2022 population, when the average noise incidence is 8 times higher than the average measles incidence. A) the noise structure is dynamical, and the seasonality is in-phase with the measles incidence. B) the noise structure is Poisson only.]
)
<tbl-optimal-thresholds-unavoidable>

#let delays = csv("supplemental-appendix_files/tables/optimal-thresholds_detection-delays.csv")
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
  caption: [Mean outbreak alert delay (days) of each testing scenario at their specific optimal thresholds, when the average noise incidence is 8 times higher than the average measles incidence. A) the noise structure is dynamical, and the seasonality is in-phase with the measles incidence. B) the noise structure is Poisson only.]
)
<tbl-optimal-thresholds-delays>


== Figures

#figure(
  image("supplemental-appendix_files/plots/optimal-thresholds_prop-outbreak-plot.svg"),
  caption: [The difference between the proportion of the time series in outbreak for outbreak detection systems under different testing rates and noise structures. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays.]
)
<fig-outbreak-proportion>

#figure(
  image("supplemental-appendix_files/plots/optimal-thresholds_alert-duration-plot.svg"),
  caption: [The difference between the alert durations for outbreak detection systems under different testing rates and noise structures. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays.]
)
<fig-alert-duration>

#figure(
  image("supplemental-appendix_files/plots/optimal-thresholds_n-alerts-plot.svg"),
  caption: [The difference between the number of alerts under different testing rates and noise structures. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate. Solid lines represent tests with 0-day turnaround times, and dashed lines represent tests with result delays.]
)
<fig-num-alerts>


#set bibliography(style: "elsevier-vancouver")

#bibliography("OD.bib")
