# DESCRIPTIONS OF THE INPUT FILES
**Author**: Johannes Borgqvist<br>
**Date**: 2021-12-16<br>
This folder contains all the input files defining the ODE models for which we want to calculate the respective symmetries. These file should be xlsx-files and they must have the same structure as the ones presented here. Currently we have included the models corresponding to Hydon's model, the DBH model, the SIR model and the Lotka-Volterra model. Again, we want to emphasise that we have tested multiple models, and for most of them the symmetry calculations never terminate due to how slow the symbolic matrix calculations in SymPy are (for a description see the main README file with relative path "../README.md"). 
