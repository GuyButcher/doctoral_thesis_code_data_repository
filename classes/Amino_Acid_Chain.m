classdef Amino_Acid_Chain < handle
    %AMINO_ACID_CHAIN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        names
        backbone_positions
        sidechain_positions
        sidechain_sizes
        sidechain_sizes_min
        sidechain_sizes_max
        selected_data
        current_figure
        pdb_file_path
    end
    
    methods
        function self = Amino_Acid_Chain(filename)
            %AMINO_ACID_CHAIN Construct an instance of this class
            %   Detailed explanation goes here
                self.names = [];
                self.backbone_positions = [];
                self.sidechain_positions = [];
                self.sidechain_sizes = [];
                self.sidechain_sizes_min = [];
                self.sidechain_sizes_max = [];
                self.current_figure = [];
                self.pdb_file_path = "";
            if (filename == "" || filename == '')
                self.names = [];
                self.backbone_positions = [];
                self.sidechain_positions = [];
                self.sidechain_sizes = [];
                self.sidechain_sizes_min = [];
                self.sidechain_sizes_max = [];
                self.current_figure = [];
                self.pdb_file_path = "";
            elseif (ischar(filename))
                self.Import_PDB_with_Sidechains(filename);
            elseif (isstring(filename))
                filename = convertStringsToChars(filename);
                self.Import_PDB_with_Sidechains(filename);
            end
        end
        function Select_Data(self,range)
            if (isvector(range))
                if (size(range,2) == 1)
                    self.selected_data = range;
                elseif (size(range,2) == 2)
                    self.selected_data = range;
                else
                    error('Amino_Acid_Chain:select_data:arguement_error_one',...
                        'Arguement, range, is not a vector of size [1 1] or [1 2]');
                end
            else
                error('Amino_Acid_Chain:select_data:arguement_error_one',...
                        'Arguement, range, is not a vector of size [1 1] or [1 2]');
            end
        end
        function output = Export_Selected(self)
            output.amino_names = self.names(self.selected_data(1):self.selected_data(2));
            output.amino_pos = self.backbone_positions(self.selected_data(1):self.selected_data(2),:);
            output.sidechain_pos = self.sidechain_positions(self.selected_data(1):self.selected_data(2),:);
        end
        function Make_Figure(self)
        %MAKE_FIGURE_RESIDUE_LABELS Displays a figure containing a
        %graphical representation of the Residue structure.
            if (isempty(self.selected_data))
                chain_length = length(self.backbone_positions);
                self.Print_Figure_Range([1 chain_length],false);
            elseif (isvector(self.selected_data) && isequal(size(self.selected_data),[1 1]))
                self.Print_Figure_Range([1 self.selected_data],false);
            elseif (isvector(self.selected_data) && ( isequal(size(self.selected_data),[1 2]) || isequal(size(self.selected_data),[2 1])))
                self.Print_Figure_Range(self.selected_data,false);
            end
        end
        function fig_handle = Make_Figure_Residue_Labels(self)
        %MAKE_FIGURE_RESIDUE_LABELS Displays a figure containing a
        %graphical representation of the Residue structure with Residue
        %Number lables.
            if (isempty(self.selected_data))
                chain_length = length(self.backbone_positions);
                fig_handle = self.Print_Figure_Range([1 chain_length],true);
            elseif (isvector(self.selected_data) && isequal(size(self.selected_data),[1 1]))
                fig_handle = self.Print_Figure_Range([1 self.selected_data],true);
            elseif (isvector(self.selected_data) && ( isequal(size(self.selected_data),[1 2]) || isequal(size(self.selected_data),[2 1])))
                fig_handle = self.Print_Figure_Range(self.selected_data,true);
            end
        end
        function output = Make_Helix_Object_V2(self)
            %Uses new Helix_Object to make the data structure
            selection(1) = self.selected_data(1);
            selection(2) = self.selected_data(2);
            
            backbone = self.backbone_positions(selection(1):selection(2),:);
            sidechain = self.sidechain_positions(selection(1):selection(2),:);
            sizes = self.sidechain_sizes(selection(1):selection(2));
            names = self.names(selection(1):selection(2));
%             sizes = [];
            output = Helix_Object(backbone,sidechain,sizes);
            output.sidechain_names = names;
            
            output.Properties(self.Generate_Sidechain_Properties([self.selected_data(1),self.selected_data(2)]));
        end
    end
    
    methods (Access = protected)
        function Import_PDB(self,filename)
            pdb_file = pdbread(filename);
            atoms = pdb_file.Model.Atom;
            for i = 1:length(atoms)
                atoms(i).AtomName = convertCharsToStrings(atoms(i).AtomName);
            end
            
            residue_count = length(unique([atoms.resSeq]));
            aminoNames = strings(residue_count,1);
            aminos = zeros(residue_count,3);
            sidechains = zeros(residue_count,3);
            
            for i = 1:residue_count
                index = [atoms.resSeq] == i;
                sequence = atoms(index);
                
                aminoNames(i) = convertCharsToStrings(sequence(2).resName);
                aminos(i,1) = sequence(2).X;
                aminos(i,2) = sequence(2).Y;
                aminos(i,3) = sequence(2).Z;
                r = find([sequence.AtomName] == "CB");
                if (r ~= 0)
                    r = sequence(r);
                    sidechains(i,1) = r.X;
                    sidechains(i,2) = r.Y;
                    sidechains(i,3) = r.Z;
                end
            end

            self.names = aminoNames;
            self.backbone_positions = aminos;
            self.sidechain_positions = sidechains;
            
        end
        function Import_PDB_with_Sidechains(self,filename)
            self.pdb_file_path = filename;
            pdb_file = pdbread(filename);
            atoms = pdb_file.Model.Atom;
            for i = 1:length(atoms)
                atoms(i).AtomName = convertCharsToStrings(atoms(i).AtomName);
            end
            
            residue_indeces = unique([atoms.resSeq]);
            residue_count = length(residue_indeces);
            aminoNames = strings(residue_count,1);
            aminos = zeros(residue_count,3);
            sidechains = zeros(residue_count,3);
            sidechains_sizes = zeros(residue_count,1);
            
            for i = 1:residue_count
                index = [atoms.resSeq] == residue_indeces(i);
                sequence = atoms(index);
                
                aminoNames(i) = convertCharsToStrings(sequence(2).resName);
                aminos(i,1) = sequence(2).X;
                aminos(i,2) = sequence(2).Y;
                aminos(i,3) = sequence(2).Z;
                r = find([sequence.AtomName] == "CB");
                if (r ~= 0)
                    
                    sequence([sequence.AtomName] == "N") = [];
                    sequence([sequence.AtomName] == "CA") = [];
                    sequence([sequence.AtomName] == "C") = [];
                    sequence([sequence.AtomName] == "O") = [];
                    num_atoms = length(sequence);
                    
                    sidechains(i,1) = sum([sequence.X])/num_atoms;
                    sidechains(i,2) = sum([sequence.Y])/num_atoms;
                    sidechains(i,3) = sum([sequence.Z])/num_atoms;
                    
                    distances = zeros(num_atoms,1);
                    
                     for j = 1:num_atoms
                        distances(j) = sqrt( ... 
                            (sidechains(i,1) - sequence(j).X)^2 + ...
                            (sidechains(i,2) - sequence(j).Y)^2 + ...
                            (sidechains(i,3) - sequence(j).Z)^2);
                    end
                    sidechains_sizes(i) = sum(distances)/length(distances);
                    sidechain_max(i) = max(distances);
                    sidechain_min(i) = min(distances);                   
                end
            end
 
            self.names = aminoNames;
            self.backbone_positions = aminos;
            self.sidechain_positions = sidechains;
            self.sidechain_sizes = sidechains_sizes;
            self.sidechain_sizes_min = sidechain_min;
            self.sidechain_sizes_max = sidechain_max;
        end
        
        function fig_handle = Print_Figure_Range(self,range,lables)
            backbone = self.backbone_positions(range(1):range(2),:);
            sidechain = self.sidechain_positions(range(1):range(2),:);
            
            self.current_figure = figure();
            fig_handle = self.current_figure;
            hold on
            plot3(backbone(:,1),backbone(:,2),backbone(:,3),'Color','blue');
            
            if (lables == [])
                lables = false;
            end
            
            if (lables)
                distance = 5;
                for i = 1:length(backbone(:,1))
                    if (mod(i,distance) == 0)
                        text(backbone(i,1),backbone(i,2),backbone(i,3),int2str(i));
                    end
                end
            end            
            
            count = 0;
            for i = 1:length(sidechain)
                if (sidechain(i-count,1) == 0 && sidechain(i-count,2) == 0 && sidechain(i-count,3) == 0)
                    backbone(i-count,:) = [];
                    sidechain(i-count,:) = [];
                    count = count + 1;
                end
            end
            
            plot3(sidechain(:,1),sidechain(:,2),sidechain(:,3),'.','Color','red');
            plot3([backbone(:,1)';sidechain(:,1)'],[backbone(:,2)';sidechain(:,2)'],[backbone(:,3)';sidechain(:,3)'],'Color','green');
            
%             if (lables == [])
%                 lables = false;
%             end
%             
%             if (lables)
%                 distance = 5;
%                 for i = 1:length(backbone(:,1))
%                     if (mod(i,distance) == 0)
%                         text(backbone(i,1),backbone(i,2),backbone(i,3),int2str(i));
%                     end
%                 end
%             end
            hold off
            axis equal
        end
        function output = Generate_Sidechain_Properties(self,range)
            gen_residue_properties;
            sidechains = self.names(range(1):range(2));
            properties = zeros(length(sidechains),3);
            for i = 1:length(sidechains)
                name = sidechains(i);
                index = find(upper(residue_properties.name_short) == name);
                properties(i,1) = residue_properties.mass(index);
                properties(i,2) = residue_properties.charge(index);
                properties(i,3) = residue_properties.hydrophobicity_kd(index);
            end
            output = properties;
        end
        function output = Generate_Sidechain_Properties_OLD(self,range)
            residue_properties = load("residueProperties.mat");
            residue_properties = residue_properties.residueProperties;
            sidechains = self.names(range(1):range(2));
            properties = zeros(length(sidechains),3);
            for i = 1:length(sidechains)
                name = sidechains(i);
                index = find(upper(residue_properties.Names) == name);
                properties(i,1) = residue_properties.Masses(index);
                properties(i,2) = residue_properties.Charges(index);
                properties(i,3) = residue_properties.Hydropathy(index);
            end
            output = properties;
        end
    end
    methods(Static)
        function output = Get_Backbone_Atom_Positions(filename,one_index,two_index)              
            pdb_file = pdbread(filename);
            atoms = pdb_file.Model.Atom;
            for i = 1:length(atoms)
                atoms(i).AtomName = convertCharsToStrings(atoms(i).AtomName);
            end
            
            residue_indeces = unique([atoms.resSeq]);
            residue_count = length(residue_indeces);
            
            if (nargin == 1)
                one_index = 1;
                two_index = residue_count;
            end
            
            
            aminoNames = strings(two_index - one_index,1);
            ca = zeros(two_index - one_index,3);
            c = zeros(two_index - one_index,3);
            o = zeros(two_index - one_index,3);
            n = zeros(two_index - one_index,3);
            
            for i = one_index:two_index
                index = [atoms.resSeq] == residue_indeces(i);
                sequence = atoms(index);
                
                aminoNames(i) = convertCharsToStrings(sequence(2).resName);
                cur_bb = sequence([sequence.AtomName] == "CA");
                ca(i,1) = cur_bb.X;
                ca(i,2) = cur_bb.Y;
                ca(i,3) = cur_bb.Z;
                cur_bb = sequence([sequence.AtomName] == "C");
                c(i,1) = cur_bb.X;
                c(i,2) = cur_bb.Y;
                c(i,3) = cur_bb.Z;
                cur_bb = sequence([sequence.AtomName] == "O");
                o(i,1) = cur_bb.X;
                o(i,2) = cur_bb.Y;
                o(i,3) = cur_bb.Z;
                cur_bb = sequence([sequence.AtomName] == "N");
                n(i,1) = cur_bb.X;
                n(i,2) = cur_bb.Y;
                n(i,3) = cur_bb.Z;
            end
            
            output.names = aminoNames(one_index:two_index);
            output.ca_positions = ca(one_index:two_index,:);
            output.c_positions = c(one_index:two_index,:);
            output.o_positions = o(one_index:two_index,:);
            output.n_positions = n(one_index:two_index,:);
        end
    end
end

