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

#show raw.where(block: true): block.with(
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
    new_title_block +
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
                body,
)


= Results
<results>
== Tables
<tables>
#figure([
#let accuracy = csv("optimal-thresholds_accuracy.csv")

#table(
  columns: 9,
    fill: (x, y) => {
      if y == 0 {gray}
      if y == 1 {gray}
    },
  [], table.cell(colspan: 2, align: center, "Test Characteristic"), table.cell(colspan: 6, align: center, "Testing Rate"),
  ..accuracy.flatten()
)
], caption: figure.caption(
position: top, 
[
The outbreak detection accuracy for each testing scenario, at the optimal thresholds.
]), 
kind: "quarto-float-tbl", 
supplement: "Table", 
)
<tbl-optimal-thresholds-accuracy>


#figure([
#let optimal_thresholds = csv("optimal-thresholds_unavoidable-cases.csv")

#table(
  columns: 9,
    fill: (x, y) => {
      if y == 0 {gray}
      if y == 1 {gray}
    },
  [], table.cell(colspan: 2, align: center, "Test Characteristic"), table.cell(colspan: 6, align: center, "Testing Rate"),
  ..optimal_thresholds.flatten()
)
], caption: figure.caption(
position: top, 
[
Unavoidable cases per annum of each testing scenario at their specific optimal thresholds, scaled up to Ghana’s 2022 population. A) the noise structure is dynamical, and the seasonality is in-phase with the measles incidence. B) the noise structure is Poisson only.
]), 
kind: "quarto-float-tbl", 
supplement: "Table", 
)
<tbl-optimal-thresholds-unavoidable>


#figure([
#let optimal_thresholds = csv("optimal-thresholds_detection-delays.csv")

#table(
  columns: 9,
    fill: (x, y) => {
      if y == 0 {gray}
      if y == 1 {gray}
    },
  [], table.cell(colspan: 2, align: center, "Test Characteristic"), table.cell(colspan: 6, align: center, "Testing Rate"),
  ..optimal_thresholds.flatten()
)
], caption: figure.caption(
position: top, 
[
Outbreak alert delay (days) of each testing scenario at their specific optimal thresholds. A) the noise structure is dynamical, and the seasonality is in-phase with the measles incidence. B) the noise structure is Poisson only.
]), 
kind: "quarto-float-tbl", 
supplement: "Table", 
)
<tbl-optimal-thresholds-delays>


== Figures
<figures>
#figure([
#box(image("supplemental-appendix_files/figure-typst/fig-outbreak-proportion-output-1.png"))
], caption: figure.caption(
position: bottom, 
[
The difference between proportion of the time series in outbreak for outbreak detection systems under different testing rates and noise structures. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-outbreak-proportion>


#figure([
#box(image("supplemental-appendix_files/figure-typst/fig-alert-duration-output-1.png"))
], caption: figure.caption(
position: bottom, 
[
The alert durations of outbreak detection systems under different testing rates and noise structures. The shaded bands illustrate the 80% central interval, and the solid/dashed lines represent the mean estimate.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-alert-duration>





#set bibliography(style: "springer-vancouver")

#bibliography("OD.bib")

