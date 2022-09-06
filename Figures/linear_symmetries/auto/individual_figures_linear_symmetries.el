(TeX-add-style-hook
 "individual_figures_linear_symmetries"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz" "crop" "convert={density=200,outext=.png}" "border=0.4cm")))
   (TeX-run-style-hooks
    "latex2e"
    "./Input/stable_r"
    "./Input/stable_theta"
    "./Input/saddle_r"
    "./Input/saddle_theta"
    "standalone"
    "standalone10"
    "pgfplots"
    "amsmath"
    "physics"
    "xcolor")
   (LaTeX-add-xcolor-definecolors
    "stab_r_1"
    "stab_r_2"
    "stab_r_3"
    "stab_theta_1"
    "stab_theta_2"
    "stab_theta_3"
    "saddle_r_1"
    "saddle_r_2"
    "saddle_r_3"
    "saddle_theta_1"
    "saddle_theta_2"
    "saddle_theta_3"))
 :latex)

