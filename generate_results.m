%% Generate all results of Static Force Sim, Pre-Collision and Collision Optimisation Sim
close all;
clear;
startup;

% Original Force
orig_force_results = Simulation.Run_All_Bundle_Static_Forces();

% Domain Names
for i = 1:13
    domain_names(i) = sprintf("R%i", i);
end

% Optimised Force
optim_results = [];
for i = 1:13
    domain_names(i) = sprintf("R%i", i);
    result = do_bundle_optimisation_rot_trans(i, helix_indeces, false);
    optim_results = [optim_results;result];
end
optim_force_results = [optim_results.final_force];

% Position and Rotation Deltas
for i = 1:13
    optim_results_rotation{i} = optim_results(i).final_rotation_configuration;
    optim_results_position{i} = optim_results(i).final_translation_configuration;
end

% Collisions Optim Sims

% Optimised Force
optim_collision_results = [];
for i = 1:13
    domain_names(i) = sprintf("R%i", i);
    result = do_bundle_optim_rot_trans_collision(i, helix_indeces, false);
    optim_collision_results = [optim_collision_results;result];
end
optim_collision_force_results = [optim_collision_results.final_force];

% Position and Rotation Deltas
for i = 1:13
    optim_collision_results_rotation{i} = optim_collision_results(i).final_rotation_configuration;
    optim_collision_results_position{i} = optim_collision_results(i).final_translation_configuration;
end

% Saving Results

table_data.names = domain_names;
table_data.original_force = orig_force_results;
table_data.optim_force = optim_force_results;
table_data.force_delta = table_data.optim_force - table_data.original_force;
table_data.optim_position = optim_results_position;
table_data.optim_rotation = optim_results_rotation;
table_data.optim_col_force = optim_collision_force_results;
table_data.optim_col_force_delta = table_data.optim_col_force - table_data.original_force;
table_data.optim_col_position = optim_collision_results_position;
table_data.optim_col_rotation = optim_collision_results_rotation;
for i = 1:13
    table_data.optim_col_position_mean(i) = mean(table_data.optim_col_position{i});
    table_data.optim_col_rotation_mean(i) = mean(table_data.optim_col_rotation{i});
end
save('results/force_sim_and_optim_collision_results', 'table_data');

%% Hydrophobic Data

close all;
clear;
startup;

for i = 1:13
    sim = Simulation.Load_Sim_Talin_Subdomain(i);
    sim.Calculate_Bundle_Axis();
    data{i} = sim.Calculate_Bundle_Hydropathy_Data();
end

count = 1;
bar_graph_data_all_hydro_groups = zeros(13,2);
bar_graph_data_all_hydro_combined = zeros(13,1);
bar_graph_data_same_hydro_groups = zeros(13,2);
bar_graph_data_same_hydro_combined = zeros(13,1);
bar_graph_n_close_hydrophobic_groups = zeros(13,2);
bar_graph_n_close_hydrophobic_combined = zeros(13,1);
bar_graph_n_close_hydrophobic_vs_hydrophilic_groups = zeros(13,2);
bar_graph_n_close_hydrophobic_vs_hydrophilic_combined = zeros(13,1);
bar_graph_n_far_hydrophilic_groups = zeros(13,2);
bar_graph_n_far_hydrophilic_combined = zeros(13,1);
bar_graph_n_far_hydrophilic_vs_hydrophobic_groups = zeros(13,2);
bar_graph_n_far_hydrophilic_vs_hydrophobic_combined = zeros(13,1);

for i = 1:length(data)
    current_data = data{i};
    for j = 1:length(current_data)
%         fprintf("Bundle: R%d, Group: %d\n", i, j);
        
        all_hydro_ratios = current_data(j).ratio_close_phobic_to_close_philic;
        all_hydro_ratios = all_hydro_ratios + current_data(j).ratio_far_philic_to_far_phobic;
        all_hydro_ratios = all_hydro_ratios + current_data(j).ratio_close_phobic_to_far_phobic;
        all_hydro_ratios = all_hydro_ratios + current_data(j).ratio_far_philic_to_close_philic;
        opposite_hydro_ratios_index(count) = all_hydro_ratios;
        
        same_hydro_ratios = current_data(j).ratio_close_phobic_to_close_philic;
        same_hydro_ratios = all_hydro_ratios + current_data(j).ratio_far_philic_to_far_phobic;
        same_hydro_ratios_index(count) = same_hydro_ratios;
        
        n_close_phobic = current_data(j).n_close_hydrophobic;
        n_close_phobic = n_close_phobic + current_data(j).n_close_hydrophobic;
        
        
        n_close_hydro_vs = (current_data(j).n_close_hydrophobic - current_data(j).n_close_hydrophilic);
        n_close_hydro_vs = n_close_hydro_vs + (current_data(j).n_close_hydrophobic - current_data(j).n_close_hydrophilic);
%         fprintf("Number: %d\n",n_close_phobic);
%         fprintf("Number: %d\n",n_close_hydro_vs);

        n_far_philic = current_data(j).n_far_hydrophilic;
        n_far_philic = n_far_philic + current_data(j).n_far_hydrophilic;
        n_far_hydro_vs = (current_data(j).n_far_hydrophilic - current_data(j).n_close_hydrophilic);
        n_far_hydro_vs = n_far_hydro_vs + (current_data(j).n_far_hydrophilic - current_data(j).n_close_hydrophilic);
%         fprintf("Number: %d\n",n_far_philic);
%         fprintf("Number: %d\n",n_far_hydro_vs);

                
        if length(current_data) == 1
            bar_graph_data_all_hydro_groups(i,2) = all_hydro_ratios;
            bar_graph_data_same_hydro_groups(i,2) = same_hydro_ratios;
            bar_graph_n_close_hydrophobic_groups(i,2) = n_close_phobic;
            bar_graph_n_close_hydrophobic_vs_hydrophilic_groups(i,2) = n_close_hydro_vs;
            bar_graph_n_far_hydrophilic_groups(i,2) = n_far_philic;
            bar_graph_n_far_hydrophilic_vs_hydrophobic_groups(i,2) = n_far_hydro_vs;
        else
            bar_graph_data_all_hydro_groups(i,j) = all_hydro_ratios;
            bar_graph_data_same_hydro_groups(i,j) = same_hydro_ratios;
            bar_graph_n_close_hydrophobic_groups(i,j) = n_close_phobic;
            bar_graph_n_close_hydrophobic_vs_hydrophilic_groups(i,j) = n_close_hydro_vs;
            bar_graph_n_far_hydrophilic_groups(i,j) = n_far_philic;
            bar_graph_n_far_hydrophilic_vs_hydrophobic_groups(i,j) = n_far_hydro_vs;
        end
        bar_graph_data_all_hydro_combined(i) = bar_graph_data_all_hydro_combined(i) + all_hydro_ratios;
        bar_graph_data_same_hydro_combined(i) = bar_graph_data_same_hydro_combined(i) + same_hydro_ratios;
        bar_graph_n_close_hydrophobic_combined(i) = bar_graph_n_close_hydrophobic_combined(i) + n_close_phobic;
        bar_graph_n_close_hydrophobic_vs_hydrophilic_combined(i) = bar_graph_n_close_hydrophobic_vs_hydrophilic_combined(i) + n_close_hydro_vs;
        bar_graph_n_far_hydrophilic_combined(i) = bar_graph_n_far_hydrophilic_combined(i) + n_far_philic;
        bar_graph_n_far_hydrophilic_vs_hydrophobic_combined(i) = bar_graph_n_far_hydrophilic_vs_hydrophobic_combined(i) + n_far_hydro_vs;
        count = count + 1;
%         disp(all_hydro_ratios);
    end
end

hydrophobic_results.raw = data;
hydrophobic_results.bar_graph_data_all_hydro_groups = bar_graph_data_all_hydro_groups;
hydrophobic_results.bar_graph_data_all_hydro_combined = bar_graph_data_all_hydro_combined;
hydrophobic_results.bar_graph_data_same_hydro_groups = bar_graph_data_same_hydro_groups;
hydrophobic_results.bar_graph_data_same_hydro_combined = bar_graph_data_same_hydro_combined;
hydrophobic_results.bar_graph_n_close_hydrophobic_groups = bar_graph_n_close_hydrophobic_groups;
hydrophobic_results.bar_graph_n_close_hydrophobic_combined = bar_graph_n_close_hydrophobic_combined;
hydrophobic_results.bar_graph_n_close_hydrophobic_vs_hydrophilic_groups = bar_graph_n_close_hydrophobic_vs_hydrophilic_groups;
hydrophobic_results.bar_graph_n_close_hydrophobic_vs_hydrophilic_combined = bar_graph_n_close_hydrophobic_vs_hydrophilic_combined;
hydrophobic_results.bar_graph_n_far_hydrophilic_groups = bar_graph_n_far_hydrophilic_groups;
hydrophobic_results.bar_graph_n_far_hydrophilic_combined = bar_graph_n_far_hydrophilic_combined;
hydrophobic_results.bar_graph_n_far_hydrophilic_vs_hydrophobic_groups = bar_graph_n_far_hydrophilic_vs_hydrophobic_groups;
hydrophobic_results.bar_graph_n_far_hydrophilic_vs_hydrophobic_combined = bar_graph_n_far_hydrophilic_vs_hydrophobic_combined;

save('results/hydrophobic_results', 'hydrophobic_results');

%% Backbone Validation Data

indices = [
    1,1; 1,2; 1,3; 1,4; 1,5;
    2,1; 2,2; 2,3; 2,4;
    3,1; 3,2; 3,3; 3,4;
    4,1; 4,2; 4,3; 4,4;
    5,1; 5,2; 5,3; 5,4; 5,5;
    6,1; 6,2; 6,3; 6,4; 6,5;
    7,1; 7,2; 7,3; 7,4; 7,5;
    8,1; 8,2; 8,3; 8,4;
    9,1; 9,2; 9,3; 9,4; 9,5;
    10,1; 10,2; 10,3; 10,4; 10,5;
    11,1; 11,2; 11,3; 11,4; 11,5;
    12,1; 12,2; 12,3; 12,4; 12,5;
    13,1; 13,2; 13,3; 13,4; 13,5;];

gradients = zeros(length(indices(:,1)),1);
parfor index = 1:length(gradients)
    gradients(index) = bp_atoms_dist_grad(indices(index,1),indices(index,2));    
end

save('results/backbone_validation_results', 'gradients');

clearvars -except helix_indeces pdb_filepaths residue_properties