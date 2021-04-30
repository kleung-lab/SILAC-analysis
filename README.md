# SILACanalysis
To process and analyze SILAC data generated from Skyline export. Used in following manuscripts - 

Chen, W., Mou, K. Y., Solomon, P., Aggarwal, R., Leung, K. K., & Wells, J. A. (2021). Large remodeling of the Myc-induced cell surface proteome in B cells and prostate cells creates new opportunities for immunotherapy. Proceedings of the National Academy of Sciences of the United States of America, 118(4), 1–11.  10.1073/pnas.2018861118

Leung, K. K., Wilson, G. M., Kirkemo, L. L., Riley, N. M., Coon, J. J., & Wells, J. A. (2020). Broad and thematic remodeling of the surfaceome and glycoproteome on isogenic cells transformed with driving proliferative oncogenes. Proceedings of the National Academy of Sciences of the United States of America, 117 (14), 7764–7775. https://doi.org/10.1073/pnas.1917947117

Wei, J., Leung, K., Truillet, C., Ruggero, D., Wells, J. A., & Evans, M. J. (2020). Profiling the Surfaceome Identifies Therapeutic Targets for Cells with Hyperactive mTORC1 Signaling. Molecular & Cellular Proteomics : MCP, 19 (2), 294–307. https://doi.org/10.1074/mcp.RA119.001785

Leung KK, Nguyen A, Shi T, Tang L, Ni X, Escoubet L, MacBeth KJ, DiMartino J, Wells JA. Multiomics of azacitidine-treated AML cells reveals variable and convergent targets that remodel the cell-surface proteome.  Proc Natl Acad Sci U S A. 2019 Jan 8;116(2):695-700. doi: 10.1073/pnas.1813666116. Epub 2018 Dec 24. 10.1073/pnas.1813666116

System requirements
1. System capable of running R (tested with v3.5.1) and RStudio (tested with v1.1.456). No non-standard hardware is required.

Installation instructions
1. To install, download "Master" folder. Install time should be minimal.

# General instructions
General steps for processing skyline output data
1. process data using CalRatio function, (need skyline export file, parameter file, SILAC_v2.R file)
2. save data so you are not re-processing data each time
3. analyze data

Detailed Instructions/Example
1. Add Skyline output data and parameter file to "Master" folder. Detailed information on the parameter file can be found below. In this example, the CD4+ hypoxia data in Figure 5 will be analyzed. The Skyline output data are contained in "CD4_Combined_Skyline_Output.csv". Note that data from both PNGase and tryptic fractions has been added. The parameter file in this example is "parameter_CD4.csv".
2. Open "Workflow.R" in RStudio.
2. Set the working directory to the "Master" folder.
4. Set source to "./Master/SILAC_Script/SILAC_v2.R".
5. Set "input" to the name of the Skyline output file (e.g. "CD4_Combined_Skyline_Output.csv").
6. Set "para" to the name of the parameter file (e.g. "parameter_CD4.csv").
7. Input desired name for output .csv file (e.g. "Script_Output_CD4.csv").
8. Set desired volcano plot properties. Default is all protein labels (custom = "", customlab = "normal"). If you would like to label a specific protein, change inputs to include protein names of interest (custom = "custom', customlab = "IL18RA", "CD70")

Detailed info on how calRatio works.

calRatio function in SILAC_v2.R
1.	concatenate skyline files for a given project (put both PNG and tryptic fraction data together, either in skyline or in R using rbind) 
2.	load skyline and parameter file
3.	 calculate l/h ratio for all peptide and all replicates
4.	 make all NaN into all NA’s
5.	 make all 0 to NA
6.	 make inf into NA   
7.	 calculate log2 of the ratio
8.	 for each unique replicate seen in skyline export, do the following for the raw data that are “used” 
    1.	 make filter (fil.dotp) for a given dotp, can filter based on either H or L ratio > dotp or BOTH > dotp
    2.	 make filter (fil.PNG) for non-deamidated peptides for PNG fractions
    3.	 make filter (fil.pro) for quantified peptides < certain values - this filter is based on fil.dotp & fil.PNG
    4.	 merge the three filters as one (fil), add to data table
    5.	 normalize data to center and/or sd (sd should be F since we always want to know absolute fold change)
    6.	 flip l/h ratio to h/l ratio if exp_label is H and write as "enrich" column
    7.	 calculate/aggregate "enrich" by Protein.Accession + Protein.Gene + exp_ref and output median, p value, and count at the replicate level
        1.	p values are expressed as -log10 value (this is a mann whitney one sample test of all peptides for a given protein compared to median of 0)
    8.  concatentate peptide level data into one big skypep file for all replicates
9.	For each unique exp_ref (with multiple replicates), calculate or aggregate "enrich" by Protein.Accession + Protein.Gene + exp_ref and output median, p value, and count
10.	export R objects to global environment
    1.	 skypep
    2.	 skypro
    3.	 skypro_ind 
    3.	 parameter

Parameter file fields
- "raw"- raw file name, has to match exactly to the “file name” from skyline report 
- "exp_name" - a unique name for each experiment, typically contain a identifier, H/L, a number for replicate, PNG/tryptic 
- "used" - T/F for analyzing, can easily use this to exclude datasets from being processed/analyzed without deleting records of the experiments, usually used to exclude replicate datasets that looks funny and have been repeated 
- "fraction" - PNG or tryp, PNG gets rid of peptides without deamidation sites,  
- "exp_label” - H/L for direction of a given replicate for treatment group 
- "exp_ref” - names to define which replicates to merge  
- "cell” - cell name 
- "user” - name 
- "project” - single common name for a project identifier 
- "dotp” - cut off of isotope dot product used 
- "dotp.fil”  - AND/OR determine whether both/either heavy and light idotp has to pass 
- "center” - T/F normalize to mean, typically T 
- “sd” - T/F normalize to sd, typically F         
- "peptides”  - minimum number of well quantified peptides (calculated by dotp) needed to be used to calculate the protein log2 ratio, 2 + for tryptic fraction, 1 + for PNGase 


```

