# Data-of-simulation-study
This folder contains RStudio files used to simulate and analyze data used in the study that was described by the manuscript "Performance Comparison of Unanchored Matching-Adjusted Indirect Comparison and Naïve Indirect Comparison for Estimating Relative Treatment Effect: A Simulation Study". This word document serves as a "read me" file.

Author's Statement: The code for this study has been modified based on the work of Antonio Remiro-Azócar et al., and specific technical details and annotations can be referred to at: Link to the research paper (https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1511),and Link to the GitHub repository (https://github.com/remiroazocar/population_adjustment_simstudy).

This code is encoded using UTF-8 and is divided into 5 main sections. 
The "**setting**" part is used to configure basic parameters for the simulation study. Depending on different simulated relative efficacies, it's separated into three files. The only distinction is in the 22nd line where "b_trt" is set to log(1), log(0.6), and log(0.3) respectively.
The "**generatedata**" section is responsible for creating all simulation study data, including individual patient data (IPD) for populations treated with Drug A and Drug B. It also generates a large sample of patients (P) for calculating marginal relative efficacy.
The "**functions**" part stores the MAIC function equations and performance metric calculation equations used in this study.
The "**indirect comparison**" segment executes MAIC and NIC, calculating and storing the relative efficacy results.
The "**analysis**" portion combines data, computes performance metrics, and generates result tables.
The "**plotting**" part creates nested loop plots.

*Note: Since the three relative efficacies are set separately, when running the calculations later, the three parts of the code need to be executed separately. The performance result data from the three parts should be merged into an Excel file for the final step of generating nested loop plots.
