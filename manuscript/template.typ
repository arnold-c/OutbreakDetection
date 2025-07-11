#let set_linespacing(
  top-edge: "cap-height",
  bottom-edge: "baseline",
  leading: 0.65em,
  first-line-indent: 0em,
  spacing: 1.2em,
  body,
) = {
  set text(top-edge: top-edge, bottom-edge: bottom-edge)
  set par(
    leading: leading,
    first-line-indent: first-line-indent,
    spacing: spacing,
  )

  body
}

#let set_headings(
  top-edge: "cap-height",
  bottom-edge: "baseline",
  leading: 0.65em,
  first-line-indent: 0em,
  spacing: 1.2em,
  h2-above-adjustment: 1.5em,
  h2-below-adjustment: 0.5em,
  h3-above-adjustment: 1em,
  h3-below-adjustment: 0.5em,
  h4-above-adjustment: 1em,
  h4-below-adjustment: 0.5em,
  body,
) = {
  show heading.where(level: 2): it => block(
    above: spacing + h2-above-adjustment,
    below: spacing + h2-below-adjustment,
  )[
    #set text(1.03em, weight: "black")
    #it.body
  ]

  show heading.where(level: 3): it => block(
    above: spacing + h3-above-adjustment,
    below: spacing + h3-below-adjustment,
  )[
    #set text(1.01em, weight: "black")
    #it.body
  ]

  show heading.where(level: 4): it => block(
    above: spacing + h4-above-adjustment,
    below: spacing + h4-below-adjustment,
  )[
    #set text(1.01em, weight: "bold", style: "italic")
    #it.body
  ]

  body
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
  // The paper's author summary.
  author-summary: [],
  // The paper's keywords.
  keywords: (),
  // The path to a bibliography file if you want to cite some external
  // works.
  // bib: none,
  // bib-title: "References",
  // Word count
  word-count: false,
  // Line numbers
  line-numbers: false,
  line-spacing: 1em,
  // Paper's content
  body,
) = {
  set document(title: title, author: authors.keys())
  set page(numbering: "1", number-align: center)
  // set text(font: ("Linux Libertine", "STIX Two Text", "serif"), lang: "en")
  show footnote.entry: it => [
    #set par(hanging-indent: 0.7em)
    #set align(left)
    #numbering(it.note.numbering, ..counter(footnote).at(it.note.location())) #it.note.body
  ]
  // show figure.caption: emph
  show figure.where(kind: table): set figure.caption(position: top)
  set table(fill: (x, y) => {
    if y == 0 { gray }
  })

  set table(fill: (x, y) => {
    if y == 0 { gray }
  })
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
        // Write author's name
        let au_name = if au_meta.keys().contains("name") { au_meta.name } else {
          au
        }
        text([#au_name])

        // Get labels of author's affiliation
        let au_inst_id = au_meta.affiliation
        let au_inst_primary = ""
        // Test whether the author belongs to multiple affiliations
        if type(au_inst_id) == array {
          // If the author belongs to multiple affiliations,
          // record the first affiliation as the primary affiliation,
          au_inst_primary = affiliations.at(au_inst_id.first())
          // and convert each affiliation's label to index
          let au_inst_index = au_inst_id.map(id => (
            inst_keys.position(key => key == id) + 1
          ))
          // Output affiliation
          super([#(au_inst_index.map(id => [#id]).join([,]))])
        } else if (type(au_inst_id) == str) {
          // If the author belongs to only one affiliation,
          // set this as the primary affiliation
          au_inst_primary = affiliations.at(au_inst_id)
          // convert the affiliation's label to index
          let au_inst_index = inst_keys.position(key => key == au_inst_id) + 1
          // Output affiliation
          super([#au_inst_index])
        }

        // Corresponding author
        if (
          au_meta.keys().contains("corresponding")
            and (
              au_meta.corresponding == "true" or au_meta.corresponding == true
            )
        ) {
          [#super[,]#footnote(numbering: "*")[
              Corresponding author. #au_name. Address:
              #if (
                not au_meta.keys().contains("address") or au_meta.address == ""
              ) {
                [#au_inst_primary.]
              }
              #if au_meta.keys().contains("email") {
                [Email: #link("mailto:" + au_meta.email.replace("\\", "")).]
              }
            ]]
        }

        if (
          au_meta.keys().contains("equal-contributor")
            and (
              au_meta.equal-contributor == "true"
                or au_meta.equal-contributor == true
            )
        ) {
          if ai == 0 {
            [#super[,]#footnote(
                numbering: "*",
              )[Equal contributors.]<ec_footnote>]
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
          text(size: 0.8em, [#super([#(ik + 1)]) #(affiliations.at(key))])
          linebreak()
        }
      ])]
  ]

  if header-title == "true" {
    header-title = title
  }

  set page(header: [
    #set text(8pt)
    #align(right)[#header-title]
  ])
  if word-count {
    import "@preview/wordometer:0.1.4": total-words, word-count
    show: word-count.with(exclude: (heading, table, figure.caption))
  }

  // Abstract and keyword block
  if abstract != [] {
    pagebreak()

    block([
      #if word-count {
        import "@preview/wordometer:0.1.4": (
          total-words,
          word-count,
          word-count-of,
        )
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

    if author-summary != [] {
      block([
        #if word-count {
          import "@preview/wordometer:0.1.4": (
            total-words,
            word-count,
            word-count-of,
          )
          text(weight: "bold", [Word count: ])
          text([#word-count-of(exclude: (heading))[#author-summary].words])
        }
        #heading([Author Summary])
        #author-summary
      ])

      v(1em)
    }
  }

  // Display contents

  pagebreak()

  let line_spacings = if line-spacing == 1 {
    (
      top-edge: "cap-height",
      bottom-edge: "baseline",
      leading: 0.65em,
      first-line-indent: 0em,
      spacing: 1.2em,
    )
  } else if line-spacing == 1.5 {
    (
      top-edge: "cap-height",
      bottom-edge: "baseline",
      leading: 1.0em,
      first-line-indent: 1em,
      spacing: 1.8em,
    )
  } else if line-spacing == 2 {
    (
      top-edge: "cap-height",
      bottom-edge: "baseline",
      leading: 1.3em,
      first-line-indent: 2em,
      spacing: 2.4em,
    )
  }

  show: set_linespacing.with(..line_spacings)
  show: set_headings.with(..line_spacings)


  // show heading.where(level: 1): it => block(above: 1.5em, below: 0.5em)[
  //   #set text(13pt, weight: "black")
  //   #it.body
  // ]
  //
  // show heading.where(level: 2): it => block(above: 1em, below: 0.5em)[
  //   #set text(11pt, weight: "black")
  //   #it.body
  // ]
  //
  //
  // set par(first-line-indent: 0em)

  if word-count {
    import "@preview/wordometer:0.1.4": total-words, word-count, word-count-of
    text(weight: "bold", [Word count: ])
    text([#word-count-of(exclude: (
        heading,
        table,
        figure.caption,
        <additional-info>,
      ))[#body].words])
  }

  if line-numbers {
    set par.line(numbering: "1")
    show figure: set par.line(numbering: none)
    body
  } else {
    body
  }
}

