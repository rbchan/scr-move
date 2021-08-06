(TeX-add-style-hook
 "camera_telemetry"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "12pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "vmargin=1in" "hmargin=1in") ("parskip" "parfill")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art12"
    "newtxtext"
    "newtxmath"
    "geometry"
    "amsmath"
    "parskip"
    "hyperref"
    "natbib"
    "bm"
    "amsfonts"
    "graphicx"
    "abstract"
    "lineno"
    "setspace"
    "caption")
   (TeX-add-symbols
    "bs"
    "bsi"
    "bx"
    "bxj"
    "by"
    "bu"
    "bui"
    "but"
    "buit"
    "ed"
    "cS")
   (LaTeX-add-labels
    "eq:ou-ct"
    "eq:ud"
    "eq:p"
    "ssec:tr-deploy"
    "ssec:n0"
    "eq:multinom"
    "eq:q-star"
    "eq:p-tilde"
    "tab:defs"
    "tab:post"
    "fig:ou-concept"
    "fig:study-area"
    "fig:bucks"
    "fig:path-cam"
    "fig:upost"
    "fig:sigma"
    "fig:hr"
    "fig:path-ud")
   (LaTeX-add-bibliographies
    "mybib"))
 :latex)

