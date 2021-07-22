(TeX-add-style-hook
 "summary_report"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "margin=1in") ("inputenc" "utf8")))
   (TeX-run-style-hooks
    "latex2e"
    "../Output/DBH_model/out"
    "../Output/Lotka_Volterra_realistic/out"
    "../Output/hydons_model/out"
    "../Output/Lotka_Volterra/out"
    "book"
    "bk10"
    "geometry"
    "inputenc"
    "amsmath"
    "amssymb"
    "amsfonts"
    "amsthm"
    "physics"))
 :latex)

