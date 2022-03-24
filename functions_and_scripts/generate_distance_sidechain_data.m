function data_output = generate_distance_sidechain_data(data_structs, input_data, amino_index)
%GENERATE_DISTANCE_SIDECHAIN_FIGURE Summary of this function goes here
%   Detailed explanation goes here

    load_output = load('functions_and_scripts/amino_acid_names.mat', 'names');
    amino_acid_names = load_output.names;

    num_domains = length(data_structs);
    num_amino_acids = length(data_structs{1});
    
    sidechain_atom_positions = input_data.sidechain_atom_positions;
    sidechain_mean_positions = input_data.sidechain_mean_positions;
    sidechain_atom_distances = input_data.sidechain_atom_distances;
    sidechain_mean_distances = input_data.sidechain_mean_distances;

    all_coordinates = [];
    all_distances = [];
    all_mean_distances = [];
    all_mean_positions = [];
    for domain_index = 1:num_domains
        if(isempty(sidechain_atom_distances{amino_index, domain_index})); continue; end
        num_residues = length(sidechain_atom_positions{amino_index,domain_index}(1,1,:));
        for residue_index = 1:num_residues
            all_coordinates = [all_coordinates; sidechain_atom_positions{amino_index, domain_index}(:,:,residue_index)];
        end
        all_distances = [all_distances;sidechain_atom_distances{amino_index,domain_index}(:)];
        all_mean_distances = [all_mean_distances;vertcat(sidechain_mean_distances{amino_index,:})];
        all_mean_positions = [all_mean_positions;vertcat(sidechain_mean_positions{amino_index,:})];
    end
    
    if(isempty(sidechain_atom_distances{amino_index, domain_index}))
        fprintf("test\n");
    end
    
    % Determine per Residue data
    mean_distance = mean(vertcat(sidechain_mean_distances{amino_index,:}));    
    data_output.distances.residue_level_mean = mean_distance;
    data_output.raw.mean_distances = all_mean_distances;
    data_output.raw.mean_positions = all_mean_positions;
    
    % Determine Stat Ranges
    all_coords_mean = mean(all_coordinates);
    all_distances_mean = mean(all_distances);
    all_distances_median = median(all_distances);
    quartiles = quantile(all_distances,[.25 .5 .75]);
    
    data_output.raw.positions = all_coordinates;
    data_output.raw.distances = all_distances;
    
    data_output.positions.mean = all_coords_mean;
    data_output.distances.first_quartile = quartiles(1);
    data_output.distances.mean = all_distances_mean;
    data_output.distances.std = std(all_distances);
    data_output.distances.median = all_distances_median;
    data_output.distances.third_quartile = quartiles(3);
    data_output.distances.min = min(all_distances);
    data_output.distances.max = max(all_distances);
    
    % Determine number of atoms that fall within these ranges.
    data_output.count.first_quartile = sum(all_distances <= quartiles(1));
    data_output.count.mean = sum(all_distances <= all_distances_mean);
    data_output.count.median = sum(all_distances <= all_distances_median);
    data_output.count.third_quartile = sum(all_distances <= quartiles(3));
    data_output.count.total = length(all_distances);
end

