function output = sidechain_atom_extraction()
    silent_startup

    load_output = load('functions_and_scripts\amino_acid_names.mat', 'names');
    amino_acid_names = load_output.names;
    num_aa_t = length(amino_acid_names);

    % Preprocess pdb atoms files...
    number_domains = length(pdb_filepaths(:));
    number_atoms = zeros(number_domains,1);
    for domain_index = 1:number_domains
        pdb_filepath = pdb_filepaths(domain_index);
        pdb_file = pdbread(pdb_filepath);
        domain_atoms(domain_index).atoms = pdb_file.Model.Atom;

        for atom_index = 1:length(domain_atoms(domain_index).atoms)
            domain_atoms(domain_index).atoms(atom_index).AtomName = convertCharsToStrings(domain_atoms(domain_index).atoms(atom_index).AtomName);
            domain_atoms(domain_index).atoms(atom_index).resName = convertCharsToStrings(domain_atoms(domain_index).atoms(atom_index).resName);
        end
        number_atoms(domain_index) = length(domain_atoms(domain_index).atoms) - 1; % To remove the OXT terminating atom.
    end

    domain_structs_array = cell(number_domains,1);

    for domain_index = 1:length(pdb_filepaths)
        amino_acid_structs_array = cell(num_aa_t,1);

        for amino_acid_index = 1:length(amino_acid_names)
            % Get all atoms of this type of amino acid
            index_of_atoms = [domain_atoms(domain_index).atoms.resName] == amino_acid_names(amino_acid_index);
            num_of_atoms = sum(index_of_atoms);
            atoms = domain_atoms(domain_index).atoms(index_of_atoms);


            % Get all unique instances of the amino acid
            residues = unique([atoms.resSeq]');
            number_residues = length(residues);
            if(number_residues == 0); continue; end
            corrected_atoms = [];

            % Determine the number of sidechain atoms in each residue
            temp_atoms = atoms([atoms.resSeq] == residues(1));
            temp_index = ~( ...
                    [temp_atoms.AtomName] == "N" | ...
                    [temp_atoms.AtomName] == "CA" | ...
                    [temp_atoms.AtomName] == "C" | ...
                    [temp_atoms.AtomName] == "O");
            temp_atoms = temp_atoms(temp_index);
            number_sidechain_atoms = length(temp_atoms);

            struct_array = [];

            for residue_index = 1:number_residues
                current_atoms = atoms([atoms.resSeq] == residues(residue_index));
                % Pick out the CA atom positions
                ca_index = [current_atoms.AtomName] == "CA";
                ca_position = horzcat( ...
                    current_atoms(ca_index).X, ...
                    current_atoms(ca_index).Y, ...
                    current_atoms(ca_index).Z);
                sidechain_atoms_index = ~( ...
                    [current_atoms.AtomName] == "N" | ...
                    [current_atoms.AtomName] == "CA" | ...
                    [current_atoms.AtomName] == "C" | ...
                    [current_atoms.AtomName] == "O" | ...
                    [current_atoms.AtomName] == "OXT");
                sidechain_atoms = current_atoms(sidechain_atoms_index);
                sidechain_positions = horzcat( ...
                    [sidechain_atoms.X]', ...
                    [sidechain_atoms.Y]', ...
                    [sidechain_atoms.Z]');
                if(isempty(sidechain_positions)); continue; end
                current_struct.raw_sidechain_positions = sidechain_positions;
                current_struct.ca_position = ca_position;
                corrected_sidechain_positions = sidechain_positions - ca_position;
                current_struct.ca_corrected_positions = corrected_sidechain_positions;
                struct_array = [struct_array,current_struct];
            end
            amino_acid_structs_array{amino_acid_index} = [amino_acid_structs_array{amino_acid_index}, struct_array];
        end 
        domain_structs_array{domain_index} = [domain_structs_array{domain_index}, amino_acid_structs_array];
    end
    output = domain_structs_array;
end