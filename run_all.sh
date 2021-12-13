#!/bin/bash
#=======================================================================
cd ./Code/
echo "============================================================"
echo "RUNNING SCRIPT 1 OUT OF 6"
echo ""
echo ""
echo "Calculating the generator for Hydon's model."
echo "The results are stored in the folder:"
echo "Output"
echo "============================================================"
echo " "
echo ""
python3 launch_symmetry_calculations.py hydons_model 2
echo "============================================================"
echo "RUNNING SCRIPT 2 OUT OF 6"
echo ""
echo ""
echo "Calculating the generator for the DBH model."
echo "The results are stored in the folder:"
echo "Output"
echo "============================================================"
echo " "
echo ""
python3 launch_symmetry_calculations.py DBH_model 2
echo "============================================================"
echo "RUNNING SCRIPT 3 OUT OF 6"
echo ""
echo ""
echo "Calculating the generator for the Lotka-Volterra (LV) model."
echo "The results are stored in the folder:"
echo "Output"
echo "============================================================"
echo " "
echo ""
python3 launch_symmetry_calculations.py LV 2
echo "============================================================"
echo "RUNNING SCRIPT 4 OUT OF 6"
echo ""
echo ""
echo "Calculating the generator for the SIR model with ansätze of degree 1."
echo "The results are stored in the folder:"
echo "Output"
echo "============================================================"
echo " "
echo ""
python3 launch_symmetry_calculations.py SIR 1
echo "============================================================"
echo "RUNNING SCRIPT 5 OUT OF 6"
echo ""
echo ""
echo "Calculating the generator for the SIR model with ansätze of degree 2."
echo "The results are stored in the folder:"
echo "Output"
echo "============================================================"
echo " "
echo ""
python3 launch_symmetry_calculations.py SIR 2
echo "============================================================"
echo "RUNNING SCRIPT 6 OUT OF 6"
echo ""
echo ""
echo "Calculating the generator for the Brusselator model."
echo "The results are stored in the folder:"
echo "Output"
echo "============================================================"
echo " "
echo ""
python3 launch_symmetry_calculations.py Brusselator 1
echo "============================================================"
echo "GENERATING REPORT"
echo ""
echo ""
echo "A summary report is generated in the Output folder"
echo "============================================================"
echo " "
echo ""
python3 generate_report.py && cd ../Output && pdflatex summary_report.tex && pdflatex summary_report.tex && evince summary_report.pdf &
echo "============================================================"
echo "EVERYTHING WAS EXECUTED SUCCESFULLY!"
echo ""
echo ""
echo "============================================================"
