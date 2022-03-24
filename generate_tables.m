%% Loading and Generating Data if not already done.

% Static Force Data
for i = 1:13
    domain_names(i) = sprintf("R%i", i);
end
forces = Simulation.Run_All_Bundle_Static_Forces();

% Optimisation Force Data
if not(isfile('results\force_sim_and_optim_collision_results.mat'))
    generate_results;
    clear;
    clc;
end

table_data = load('results\force_sim_and_optim_collision_results.mat');
table_data = table_data.table_data;

%% Static Force Simulation Results Table

clearvars tab output;

for i = 1:13
    domain_names(i) = sprintf("R%i", i);
end
forces = Simulation.Run_All_Bundle_Static_Forces();

tab.bundle = domain_names';
tab.force = round(forces,3,'significant')' * 10^12;

output = struct2table(tab);
output.Properties.VariableNames = {'Bundle','Force, \si{\pico\newton}'};
save('tables/static_force_table', 'output');
tableLatex(output, true, false, 'tables/static_force_table');

%% Optimisation Results Pre-Collision System



clearvars tab output;

tab.bundle = table_data.names';
tab.orig_force = round(table_data.original_force, 3, 'significant')' * 10^12;
tab.optim_force = round(table_data.optim_force, 3, 'significant')' * 10^12;
tab.force_delta = round(table_data.force_delta, 3, 'significant')' * 10^12;

output = struct2table(tab);
output.Properties.VariableNames = {...
    'Bundle', ...
    'Original Force, \si{\pico\newton}', ...
    'Optimised Force, \si{\pico\newton}', ...
    'Force Delta, \si{\pico\newton}'};

save('tables/optim_pre_collision_table', 'output');
tableLatex(output, true, true, 'tables/optim_pre_collision_table');

%% Optimisation Results After Collision System

clearvars tab output;

tab.bundle = table_data.names';
tab.orig_force = round(table_data.original_force, 3, 'significant')' * 10^12;
tab.optim_col_force = round(table_data.optim_col_force, 3, 'significant')' * 10^12;
tab.force_delta = round(table_data.optim_col_force_delta, 3, 'significant')' * 10^12;

output = struct2table(tab);
output.Properties.VariableNames = {...
    'Bundle', ...
    'Original Force, \si{\pico\newton}', ...
    'Optimised Force, \si{\pico\newton}', ...
    'Force Delta, \si{\pico\newton}'};

save('tables/optim_after_collision_table', 'output');
tableLatex(output, true, true, 'tables/optim_after_collision_table');

%% Table Three Generation

clearvars tab output;

tab.bundle = table_data.names';
tab.position_delta = process_cells_for_table(table_data.optim_col_position);
tab.rotation_delta = process_cells_for_table(table_data.optim_col_rotation);
tab.position_delta_mean = round(table_data.optim_col_position_mean, 3, 'significant')';
tab.rotation_delta_mean = round(table_data.optim_col_rotation_mean, 3, 'significant')';

output = struct2table(tab);
output.Properties.VariableNames = {...
    'Bundle', ...
    'Position Delta, \si{\angstrom}', ...
    'Rotation Delta, Degrees', ...
    'Mean Position Delta, \si{\angstrom}', ...
    'Mean Rotation Delta, Degrees'};

save('tables/validation_after_collision_table', 'output');
tableLatex(output, true, true, 'tables/validation_after_collision_table');

clearvars -except helix_indeces pdb_filepaths residue_properties