# Readme
### Repository for codebase to generate results, tables, and figures for the Doctoral Thesis titled: “Improving force values for the unfolding thresholds of talin’s rod subdomains through force-domain simulations designed for the spatial and temporal ranges of talin.''

#### Key Files
* `startup.m` generates any prerequisite data `.mat` files and add relevant folders to Matlab's current path. This needs to be run before any other script.
* `generate_results.m` Runs all simulation required to generate the results' data. Output is stored in the `results\` folder. This is the exact code used to generate the data for the thesis results.
* `generate_tables.m` Generates all main results tables. Output is stored in the `tables\` folder. This is the exact code used to generate the main table in the thesis.
* `generate_figures.m` Generates all the figures that are derived from results data. Output is stored in the `figures\` folder. This is the exact code used to generate the figures shown in the thesis.
* `generate_optim_comparision_figure.m` This figures generation is separated from the main file as it is computationally expensive and can stall on machines with less than 16 GB RAM. Generates optimisation comparison tiled figure. Output is stored in the `figures\` folder. This is the exact code used to generate the figures shown in the thesis.


#### Folders
* `classes/` holds the class definitions and member functions for the key classes that handle all simulation objects and that generate simulation objects from pdb data.
* `funtions_and_scripts/` holds helper functions.
* `pdb_data/` contains the local copy of the pdb files used within the simulation.
 
#### System Requirements
* MATLAB 2020b or newer
    * Bioinformatics Toolbox
    * Optimization Toolbox
    * Global Optimization Toolbox
    * Parallel Computing Toolbox
* Recommendation: System with 16-32 GB of Memory and 6-12 CPU Cores with Hyperthreading.