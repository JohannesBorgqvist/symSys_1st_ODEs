#=================================================================================
#=================================================================================
# Script:"generate_report"
# Date: 2021-06-02
# Implemented by: Johannes Borgqvist
# Description:
# The script generates a summary report
# from all the calculations that have
# been completed.
#=================================================================================
#=================================================================================
# Import Libraries
#=================================================================================
#=================================================================================
# To create a data folder
import os
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# FUNCTION 1: "create_summary_document"
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# The function creates a LaTeX-document which summarises all the logged reports from all the outputs that exist in the Output folder.
def create_summary_document():
    # Open the report with the "w" option so that
    # we create a new file (we do not want to append
    # this time)
    f = open("../Output/summary_report.tex", "w")
    str_temp = "\\documentclass{book}\n" # The document class
    str_temp += "% General document formatting\n" # Clarifying comment
    str_temp += "\\usepackage[margin=1in]{geometry}\n" # The margins of the document
    str_temp += "\\usepackage[utf8]{inputenc}\n" # Define characters
    str_temp += "% Related to math\n" # Clarifying comment
    str_temp += "\\usepackage{amsmath,amssymb,amsfonts,amsthm}\n" # Math packages
    str_temp += "% For writing derivatives easily\n" # Clarifying comment
    str_temp += "\\usepackage{physics}\n" # The margins of the document    
    str_temp += "\\begin{document}\n" # The document starts
    str_temp += "\\title{\\textbf{Summary of symmetry calculations}}\n" # Add a title
    str_temp += "\\date{\\today}\n" # Add the date
    str_temp += "\\maketitle\n\\tableofcontents\n\\clearpage\n" # Create the title
    # Write the entire document
    f.write("%s"%(str_temp))
    # Find all reports to include in the document
    r = list_subdirectories("../Output/")
    # Loop over all files and include them
    for i in range(len(r)):
        # File name
        file_name = r[i]
        # Add a chapter title
        chapter_name = file_name
        chapter_name = chapter_name.replace("../Output/","")
        chapter_name = chapter_name.replace("/out.tex","")
        chapter_name = chapter_name.replace("_","\_")
        chapter_name = "\\textsc{" + chapter_name + "}"
        # Add the chapter
        str_temp = "\\chapter{" + chapter_name + "}\n"
        # Add the input file
        str_temp += "\\input{" + file_name + "}\n\n"
    # Now we can include all files
    f.write("%s"%(str_temp))
    # The document ends
    str_temp = "\\end{document}" # The document ends
    f.write("%s"%(str_temp))    
    # Close our report
    f.close()
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# FUNCTION 1: "create_summary_document"
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# The function takes a directory as input and then it lists all subdirectories and all files in that folder in that folder. 
def list_subdirectories(dir):                                                            # All the files                                  
    r = []                                                                               # All the subdirectories     
    subdirs = [x[0] for x in os.walk(dir)]                                               # Loop over subdirectories   
    for subdir in subdirs:                                                                   # Loop over files        
        files = os.walk(subdir).__next__()[2]                                                    # If we have any files in the directory?
        if (len(files) > 0):                                                                     # Loop over files
            for file in files:
                # Only save the "out.tex" files
                if file == "out.tex":
                    # Save each directory
                    r.append(os.path.join(subdir, file))                                                                         
    return r


# Call the function
create_summary_document()
