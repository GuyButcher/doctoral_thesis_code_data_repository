clc
fprintf("---: PDB and Talin Sim :---\n")
fprintf("Executing Startup:...\n");

fprintf("Adding 'classes' path.\n");
addpath("classes/");

fprintf("Adding 'functions_and_scripts' path.\n");
addpath("functions_and_scripts/");

fprintf("Adding 'pdb_data' path.\n");
addpath("pdb_data/");

fprintf("Checking if 'results/' directory exists.\n");
if(not(isfolder('results/')))
    fprintf("Creating 'results' directory.\n");
    mkdir('results');
end

fprintf("Adding 'results' path.\n");
addpath("results/");

fprintf("Checking if 'tables/' directory exists.\n");
if(not(isfolder('tables')))
    fprintf("Creating 'tables/' directory.\n");
    mkdir('tables');
end

fprintf("Adding 'tables' path.\n");
addpath("tables/");

fprintf("Checking if 'figures/' directory exists.\n");
if(not(isfolder('figures')))
    fprintf("Creating 'figures/' directory.\n");
    mkdir('figures');
end

fprintf("Adding 'figures' path.\n");
addpath("figures/");

fprintf("Running 'gen_residue_properties'.\n");
gen_residue_properties;

fprintf("Loading variable 'helix_indeces' to Workspace.\n");
load("helix_indeces.mat");

fprintf("Loading variable 'pdb_filepaths' to Workspace.\n")
pdb_filepaths = Get_Talin_PDB_Filepaths();

warning('off', 'MATLAB:print:ContentTypeImageSuggested');

fprintf("Startup Complete.\n\n\n");