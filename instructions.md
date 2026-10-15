# Instructions for Running This Study

This document provides detailed instructions for running the different components of this study.

---
### Steps

1)  Download this entire repository (you can download as a zip folder using Code → Download ZIP, or you can use GitHub Desktop).
2)  Open the project `EUPAS1000000443.Rproj` from the repository directory in RStudio (when inside the project, you will see its name on the top-right of your RStudio session)
3)  Open the `codeToRun.R` file - this is the only file you should need to interact with.
4)  Load the `renv` package, run `restore()`, and then load the required study packages.
5)  Populate the required parameters, e.g. add database specific parameters (name of database, schema name with OMOP data, schema name to write results, table name stem for results to be saved in the result schema), minimum cell count... To populate the arguments of `cdmFromCon()` from the loaded `CDMConnector` package, see [Examples](https://darwin-eu.github.io/CDMConnector/articles/a04_DBI_connection_examples.html) for different dbms.
6)  Run the rest of the content. Results will be generated in the `results` folder.
