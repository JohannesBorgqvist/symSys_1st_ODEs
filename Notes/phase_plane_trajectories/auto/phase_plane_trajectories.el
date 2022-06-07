(TeX-add-style-hook
 "phase_plane_trajectories"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "12pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "margin=1in") ("inputenc" "utf8")))
   (TeX-run-style-hooks
    "latex2e"
    "./Input/Intro"
    "./Input/fibre_preserving_symmetries"
    "./Input/LV"
    "./Input/BZ"
    "./Input/Brusselator"
    "article"
    "art12"
    "geometry"
    "inputenc"
    "amsmath"
    "amssymb"
    "amsfonts"
    "amsthm"
    "physics"))
 :latex)

