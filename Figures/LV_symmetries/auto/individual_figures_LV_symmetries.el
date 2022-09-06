(TeX-add-style-hook
 "individual_figures_LV_symmetries"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz" "crop" "convert={density=200,outext=.png}" "border=0.4cm")))
   (TeX-run-style-hooks
    "latex2e"
    "./Input/tau_trans"
    "./Input/u_trans"
    "./Input/v_trans"
    "standalone"
    "standalone10"
    "pgfplots"
    "amsmath"
    "physics"
    "xcolor")
   (LaTeX-add-xcolor-definecolors
    "clr_1"
    "clr_2"
    "clr_3"
    "clr_4"))
 :latex)

