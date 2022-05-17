(TeX-add-style-hook
 "output_julia"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "margin=0.7in") ("parskip" "parfill") ("inputenc" "utf8")))
   (TeX-run-style-hooks
    "latex2e"
    "lin_syms"
    "article"
    "art10"
    "geometry"
    "parskip"
    "inputenc"
    "amsmath"
    "amssymb"
    "amsfonts"
    "amsthm"))
 :latex)

