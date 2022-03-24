function [output] = process_sidechain_data(data_structs)
%PROCESS_FOR_AMINO_ACID Summary of this function goes here
%   Detailed explanation goes here

    num_domains = length(data_structs);
    num_amino_acids = length(data_structs{1});

    sidechain_atom_positions = cell(num_amino_acids,num_domains);
    sidechain_mean_positions = cell(num_amino_acids,num_domains);
    sidechain_atom_distances = cell(num_amino_acids,num_domains);
    sidechain_mean_distances = cell(num_amino_acids,num_domains);

    for amino_index = 1:num_amino_acids
        for domain_index = 1:num_domains            
            current_residue = data_structs{domain_index}{amino_index};

            if(isempty(current_residue)); continue; end

            num_residues = length(current_residue);
            num_atoms = length(current_residue(1).raw_sidechain_positions(:,1));

            atom_positions = zeros(num_atoms,3,num_residues);
            mean_positions = zeros(num_residues,3);
            atom_distances = zeros(num_atoms,num_residues);
            mean_distances = zeros(num_residues,1);

            for residue_index = 1:num_residues
                ca_position = current_residue(residue_index).ca_position;
                raw_position = current_residue(residue_index).raw_sidechain_positions;
%                 mean_raw_position = mean(raw_position);

                normalised_position = raw_position - ca_position;
                mean_normalised_position = mean(normalised_position);            
                if(length(mean_normalised_position) == 1)
                    mean_normalised_position = normalised_position;
                end

                vector_ca_to_mean_position = (mean_normalised_position) / vecnorm(mean_normalised_position);
                vector_x_axis = [1 0 0];
                rot = VecToVecRotation(vector_ca_to_mean_position, vector_x_axis);
                % Test check
%                 test = normalised_position - mean_normalised_position;
%                 test = test + mean_normalised_position;
%                 test2 = (rot * test')';
                transformed_atom_position = (rot * normalised_position')';

                mean_transformed_atom_position = mean(transformed_atom_position);
                mean_position = mean(transformed_atom_position);
                if(length(mean_position) == 1)
                    mean_position = transformed_atom_position;
                end
                
                
                atom_positions(:,:,residue_index) = transformed_atom_position - mean_transformed_atom_position;
                mean_positions(residue_index,:) = mean_position;
                
                atom_distances(:,residue_index) = vecnorm((transformed_atom_position - mean_transformed_atom_position)')';
                mean_distances(residue_index) = mean(vecnorm((transformed_atom_position - mean_transformed_atom_position)')');
            end
            output.sidechain_atom_positions{amino_index,domain_index} = atom_positions;
            output.sidechain_mean_positions{amino_index,domain_index} = mean_positions;
            output.sidechain_atom_distances{amino_index,domain_index} = atom_distances;
            output.sidechain_mean_distances{amino_index,domain_index} = mean_distances;
        end
    end

end

