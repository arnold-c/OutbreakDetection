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
  $if(title)$
    title: "$title$",
  $endif$
  $if(header-title)$
    header-title: "$header-title$",
  $endif$
  $if(by-author)$
    authors: (
    $for(by-author)$
      $if(it.name.literal)$
          "$it.name.literal$": (
            affiliation: ($for(it.affiliations)$"$it.id$"$sep$, $endfor$),
            $if(it.attributes.corresponding)$corresponding: "true",$endif$
            $if(it.attributes.equal-contributor)$equal-contributor: "true",$endif$
            $if(it.email)$email: "$it.email$",$endif$
      ),$endif$
    $endfor$
    ),
  $endif$
  $if(by-affiliation)$
    affiliations: (
    $for(by-affiliation)$
      "$it.id$": "$it.address$"$sep$,
    $endfor$
    ),
  $endif$
  $if(date)$
    date: $date$,
  $endif$
  $if(bibliography)$
    bib: "$bibliography$",
    $if(bibliographytitle)$
      bib-title: "$bibliographytitle$",
    $endif$
  $endif$
  $if(keywords)$
    keywords: ($for(keywords)$"$it$"$sep$,$endfor$),
  $endif$
  $if(abstract)$
    abstract: [$abstract$],
  $endif$
  $if(word-count)$
    word-count: $word-count$,
  $endif$
  $if(line-numbers)$
    line-numbers: $line-numbers$,
  $endif$
  body,
)
