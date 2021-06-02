#!/bin/bash
#=======================================================================
cd ./Code/
echo "============================================================"
echo "RUNNING SCRIPT 1 OUT OF 1"
echo ""
echo ""
echo "Calculating the generator for Hydons_model."
echo "The results are stored in the folder:"
echo "Output"
echo "============================================================"
echo " "
echo ""
python3 launch_symmetry_calculations.py hydons_model 1
echo "============================================================"
echo "GENERATING REPORT"
echo ""
echo ""
echo "A summary report is generated in the Output folder"
echo "============================================================"
echo " "
echo ""
python3 generate_report.py
cd ../Output
pdflatex summary_report.tex

