addpath("classes/");

addpath("functions_and_scripts/");

addpath("pdb_data/");

addpath("results/");

addpath("tables/");

addpath("figures/");

gen_residue_properties;

load("helix_indeces.mat");

pdb_filepaths = Get_Talin_PDB_Filepaths();

warning('off', 'MATLAB:print:ContentTypeImageSuggested');