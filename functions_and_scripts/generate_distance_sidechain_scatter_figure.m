function fig_handle = generate_distance_sidechain_scatter_figure(input_data,data_output, amino_index, num_domains, az, el)

    load_output = load('functions_and_scripts/amino_acid_names.mat', 'names');
    amino_acid_names = load_output.names;
    
    sidechain_atom_positions = input_data.sidechain_atom_positions;
    sidechain_atom_distances = input_data.sidechain_atom_distances;

    fig_handle = figure();
    fig_handle.CurrentAxes = axes();
    hold on
    for domain_index = 1:num_domains
        if(isempty(data_output.raw.positions)); continue; end
        positions = sidechain_atom_positions{amino_index,domain_index};
        if(isempty(positions)); continue; end;
        num_residues = length(sidechain_atom_positions{amino_index,domain_index}(1,1,:));
        for residue_index = 1:num_residues
            scatter3(fig_handle.CurrentAxes, ...
                sidechain_atom_positions{amino_index,domain_index}(:,1,residue_index), ...
                sidechain_atom_positions{amino_index,domain_index}(:,2,residue_index), ...
                sidechain_atom_positions{amino_index,domain_index}(:,3,residue_index),'.');
        end
    end
    
    title(fig_handle.CurrentAxes, sprintf("Amino Acid: %s", amino_acid_names(amino_index)));
    xlabel(fig_handle.CurrentAxes,"X Axis, \r{A}", 'Interpreter', 'latex');
    ylabel(fig_handle.CurrentAxes, "Y Axis, \r{A}", 'Interpreter', 'latex');
    zlabel(fig_handle.CurrentAxes, "Z Axis, \r{A}", 'Interpreter', 'latex');
    
    if(~isempty(sidechain_atom_distances{amino_index, domain_index}))
        apply_stats_ranges_overlay(fig_handle.CurrentAxes, data_output, az, el);
    end
    hold off 
end

