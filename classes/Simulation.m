 classdef Simulation < handle
    %SIMULATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        chain
        helices
        baseFrame
        data
        simulation_Constants
    end
    
    methods
        function self = Simulation()
            self.chain = [];
            self.helices = [];
            self.baseFrame = Coordinate_Frame();
            self.simulation_Constants.constant = 1 / (4 * pi * 8.854e-12);
            self.simulation_Constants.charge_unit = 1.6e-19;
            self.simulation_Constants.distance_unit = 1e-9;
        end
        
        function Load_Amino_Acid_PDB_File(self,filepath)
            self.chain = Amino_Acid_Chain(filepath);
        end
        function Load_Helix(self, domain_index, helix_index, helix_indeces)
            self.Create_Helix(...
                helix_indeces(domain_index,helix_index,1),...
                helix_indeces(domain_index,helix_index,2));
        end
        function Load_Helices_Order(self, domain_index, helix_inclusion_indeces, helix_indeces)
            if (length(helix_inclusion_indeces) > 5); return; end
            if (length(helix_inclusion_indeces) < 1); return; end
            if (max(helix_inclusion_indeces) > 5); return; end
            if (min(helix_inclusion_indeces) < 1); return; end
            
            for index = 1:length(helix_inclusion_indeces)
                self.Load_Helix(domain_index,helix_inclusion_indeces(index),helix_indeces);
            end
        end
        function output = Create_Helix(self, range1, range2)
            selectedRange = -1;
            switch nargin
                case 2
                    if (isvector(range1))
                        if (size(range1,2) == 1)
                            selectedRange = range1;
                        elseif (size(range,2) == 2)
                            selectedRange = range1;
                        else
                            error('Amino_Acid_Chain:select_data:arguement_error_one',...
                        'Arguement, range, is not a vector of size [1 1] or [1 2]');
                        end
                    else
                        error('Amino_Acid_Chain:select_data:arguement_error_one',...
                        'Arguement, range, is not a vector of size [1 1] or [1 2]');
                    end
                case 3
                    selectedRange = [range1, range2];
                otherwise
                    error('Need at least one input');
            end
            if (selectedRange == -1)
                error('Need at least one input');
            end
            index = length(self.helices) + 1;
            self.chain.Select_Data(selectedRange); 
            self.helices = [self.helices,self.chain.Make_Helix_Object_V2()];
            output = index;
        end
        
        % Figure Generation
        function fig_handle = Make_Figure(self)
            fig_handle = self.baseFrame.Generate_Complete_Figure();
        end
        function fig_handle = Make_Figure_With_Attributes(self,attributes)
            if(isempty(attributes))
                attributes = self.Generate_Attributes_Struct();
            end
            fig_handle = self.baseFrame.Generate_Complete_Figure_With_Attributes(attributes);
        end
        function fig_handle = Make_Figure_Dots(self)
            fig_handle = self.baseFrame.Generate_Complete_Figure_Dots();
        end
        function fig_handle = Make_Figure_Hydrophobic(self)
            fig_handle = self.baseFrame.Generate_Complete_Figure_Hydrophobic();
        end
        function fig_hanle =  Make_Figure_Coordinate_Frames(self)
            fig_hanle = self.baseFrame.Generate_Coordinate_Frame_Graph();
        end
        
        % Bundle Axis Calculation
        function positions = Calculate_Bundle_Axis(self)
            helix_count = length(self.helices);
            
            he = zeros(helix_count, 3);
            hs = zeros(helix_count, 3);
            
            for i = 1:helix_count
                if (mod(i,2) == 1)
                    hs(i,:) = self.helices(i).parent_frame.Get_Startpoint_World_Frame();
                    he(i,:) = self.helices(i).parent_frame.Get_Endpoint_World_Frame();
                else
                    hs(i,:) = self.helices(i).parent_frame.Get_Endpoint_World_Frame();
                    he(i,:) = self.helices(i).parent_frame.Get_Startpoint_World_Frame();
                end
            end
            
            if (helix_count ~= 5)
                hsmean = mean(hs(1:4,:));
                hemean = mean(he(1:4,:));

                pos_one = hsmean;
                pos_two = hemean;
                
                positions = [pos_one;pos_two];
            else
                midpoints = zeros(5,3);
                for i = 1:5
                    midpoints(i,:) = (hs(i,:) + he(i,:))/2;
                end
                midpoints_distances = zeros(4,1);
                for i = 1:4
                    midpoints_distances(i) = norm( midpoints(1,:) - midpoints(1 + i,:) );
                end
                
                % determine the two closest helices by midpoint positions
                indices = ismember(midpoints_distances,min(midpoints_distances));
                indices = indices | ismember(midpoints_distances,min(midpoints_distances(~indices)));
                helices_index = 2:5;
                helices_index = helices_index(indices);
                
                self.data.three_grouped_helices = [1,helices_index];

                % 3 helix part

                hs3mean = mean(hs([1,helices_index(1),helices_index(2)],:));
                he3mean = mean(he([1,helices_index(1),helices_index(2)],:));

                pos3_one = hs3mean;
                pos3_two = he3mean;

                % 4 helix part

                hs4mean = mean(hs(2:5,:));
                he4mean = mean(he(2:5,:));

                pos4_one = hs4mean;
                pos4_two = he4mean;
                
                pos_one = [pos3_one;pos3_two];
                pos_two = [pos4_one;pos4_two];
                
                positions = zeros(2,3,2);
                positions(:,:,1) = pos_one;
                positions(:,:,2) = pos_two;
            end
            self.data.bundle_axis_positions = positions;            
        end
        function fig_handle = Calculate_Bundle_Axis_Plot(self)
            fig = self.Make_Figure_Dots();
            fig = fig.CurrentAxes;
            
            if (isfield(self.data,'bundle_axis_positions'))
                positions = self.data.bundle_axis_positions;
            else
                positions = self.Calculate_Bundle_Axis();
            end
            
            helix_count = length(self.helices);
            
            if (helix_count ~= 5) % 4 Helix Bundle
                pos_one = positions(1,:);
                pos_two = positions(2,:);
                direction = pos_two - pos_one;
                expos_one = pos_one + (-0.5*direction);
                expos_two = pos_two + (0.5*direction);

                hold(fig,'on')
                plot3(fig, [pos_one(1),pos_two(1)],[pos_one(2),pos_two(2)],[pos_one(3),pos_two(3)],'k-*','LineWidth',2)
                plot3(fig, [expos_one(1),expos_two(1)],[expos_one(2),expos_two(2)],[expos_one(3),expos_two(3)],'k-','LineWidth',2)
                hold(fig,'off')
            else                  % 5 Helix Bundle
                pos3_one = positions(1,:,1);
                pos3_two = positions(2,:,1);
                direction3 = pos3_two - pos3_one;

                expos3_one = pos3_one + (-0.5*direction3);
                expos3_two = pos3_two + (0.5*direction3);

                hold(fig,'on');
                plot3(fig, [pos3_one(1),pos3_two(1)],[pos3_one(2),pos3_two(2)],[pos3_one(3),pos3_two(3)],'k-*','LineWidth',2);
                plot3(fig, [expos3_one(1),expos3_two(1)],[expos3_one(2),expos3_two(2)],[expos3_one(3),expos3_two(3)],'k-','LineWidth',2);
                hold(fig,'off');

                pos4_one = positions(1,:,2);
                pos4_two = positions(2,:,2);

                % Get extended axis positions for plot.

                direction4 = pos4_two - pos4_one;

                expos4_one = pos4_one + (-0.5*direction4);
                expos4_two = pos4_two + (0.5*direction4);

                hold(fig,'on');
                plot3(fig, [pos4_one(1),pos4_two(1)],[pos4_one(2),pos4_two(2)],[pos4_one(3),pos4_two(3)],'k-*','LineWidth',2);
                plot3(fig, [expos4_one(1),expos4_two(1)],[expos4_one(2),expos4_two(2)],[expos4_one(3),expos4_two(3)],'k-','LineWidth',2);
                hold(fig,'off');
            end
            fig_handle = fig;
        end
        
        function stored_distances = Calculate_Axis_Sidechain_Distances(self)
            if(~isfield(self.data,'bundle_axis_positions'))
                error("Bundle axis not calculated")
            end
            
            if (length(self.helices) ~= 5)
                stored_distances = cell(4,1);
                for i = 1:length(stored_distances)
                    stored_distances{i} = zeros(self.helices(i).Get_Number_Sidechains(),1);
                end
                pos_one = self.data.bundle_axis_positions(1,:,1);
                pos_two = self.data.bundle_axis_positions(2,:,1);
                for i = 1:length(self.helices)
                    for j = 1:self.helices(i).Get_Number_Sidechains()
                        sc_pos = self.helices(i).parent_frame.Get_Position_World_Frame(self.helices(i).Get_Sidechain_Position(j)')';
                        i_pos = Point_Line_Intersection(pos_one, pos_two, sc_pos);
                        stored_distances{i}(j,:) = norm(i_pos - sc_pos);
                    end
                end
            else
                % First 3 are the 3 helix axis, the remaining four are of
                % the four helix axis.
                stored_distances = cell(7,1);
                index = self.data.three_grouped_helices;
                potential_index_values = [1,2,3,4,5];
                other_indeces = ismember(potential_index_values, index);
                other_index = potential_index_values(~other_indeces);
                managed_full_index = horzcat(index, other_index, index(2), index(3));
                index = managed_full_index;
                % Calc for 3 helix component
                for i = 1:3
                    stored_distances{i} = zeros(self.helices(index(i)).Get_Number_Sidechains(),1);
                end
                pos_one = self.data.bundle_axis_positions(1,:,1);
                pos_two = self.data.bundle_axis_positions(2,:,1);
                for i = 1:3
                    for j = 1:self.helices(index(i)).Get_Number_Sidechains()
                        sc_pos = self.helices(index(i)).parent_frame.Get_Position_World_Frame(self.helices(index(i)).Get_Sidechain_Position(j)')';
                        i_pos = Point_Line_Intersection(pos_one, pos_two, sc_pos);
                        stored_distances{i}(j,:) = norm(i_pos - sc_pos);
                    end
                end
                % Calc for 4 helix component
                for i = 4:7
                    stored_distances{i} = zeros(self.helices(index(i)).Get_Number_Sidechains(),1);
                end
                pos_one = self.data.bundle_axis_positions(1,:,2);
                pos_two = self.data.bundle_axis_positions(2,:,2);
                for i = 4:7
                    for j = 1:self.helices(index(i)).Get_Number_Sidechains()
                        sc_pos = self.helices(index(i)).parent_frame.Get_Position_World_Frame(self.helices(index(i)).Get_Sidechain_Position(j)')';
                        i_pos = Point_Line_Intersection(pos_one, pos_two, sc_pos);
                        stored_distances{i}(j,:) = norm(i_pos - sc_pos);
                    end
                end
            end
            
            
            
            self.data.sidechain_axis_distances = stored_distances;
            
        end
        function positions = Get_Bundle_Axis(self)
            if(isfield(self.data,'bundle_axis_positions'))
                positions = self.data.bundle_axis_positions;
            else
                positions = [];
                fprintf("Unable to get bundle axis positions.\nNo positions stored.\n");
            end
        end
        
        % Hydropathy Bundle Data
        function data = Calculate_Bin_Ratio_Hydro_Data(~,results)
            dist_min = min(results(:,1));
            dist_max = max(results(:,1));
            dist_mid = dist_min + ((dist_max - dist_min)/2);
            hydropathy_min = min(results(:,2));
            hydropathy_max = max(results(:,2));
            hydropathy_mid = hydropathy_max - abs(hydropathy_min);

            close_index = results(:,1) < dist_mid;
            far_index = results(:,1) >= dist_mid;
            n_close = length(results(close_index,1));
            n_far = length(results(far_index,1));

            bin_dist_close = results(close_index,:);
            bin_dist_far = results(far_index,:);
            
            phobic_index = results(:,2) < hydropathy_mid;
            philic_index = results(:,2) >= hydropathy_mid;
            n_phobic = length(results(phobic_index, 2));
            n_philic = length(results(philic_index, 2));
            
            bin_phobic = results(phobic_index,:);
            bin_philic = results(philic_index,:);

            close_hydrophobic_index = results(:,1) < dist_mid & results(:,2) >= hydropathy_mid;
            close_hydrophilic_index = results(:,1) < dist_mid & results(:,2) < hydropathy_mid;

            bin_close_hydrophobic = results(close_hydrophobic_index,:);
            bin_close_hydrophilic = results(close_hydrophilic_index,:);
            n_close_hydrophobic = length(bin_close_hydrophobic(:,1));
            n_close_hydrophilic = length(bin_close_hydrophilic(:,1));

            far_hydrophobic_index = results(:,1) >= dist_mid & results(:,2) >= hydropathy_mid;
            far_hydrophilic_index = results(:,1) >= dist_mid & results(:,2) < hydropathy_mid;

            bin_far_hydrophobic = results(far_hydrophobic_index,:);
            bin_far_hydrophilic = results(far_hydrophilic_index,:);
            n_far_hydrophobic = length(bin_far_hydrophobic(:,1));
            n_far_hydrophilic = length(bin_far_hydrophilic(:,1));

            ratio_close_phobic_to_close_philic = n_close_hydrophobic / n_close_hydrophilic;
            ratio_far_philic_to_far_phobic = n_far_hydrophilic / n_far_hydrophobic;
            ratio_close_phobic_to_far_phobic = n_close_hydrophobic / n_far_hydrophobic;
            ratio_far_philic_to_close_philic = n_far_hydrophilic / n_close_hydrophilic;
            
            ratio_n_close_phobic_to_n_phobic = n_close_hydrophobic / n_phobic;
            ratio_n_close_phobic_to_total = n_close_hydrophobic / (n_phobic + n_philic);
            
            data.is_three_part = 0;
            data.n_close_hydrophobic = n_close_hydrophobic;
            data.n_close_hydrophilic = n_close_hydrophilic;
            data.n_far_hydrophobic = n_far_hydrophobic;
            data.n_far_hydrophilic= n_far_hydrophilic;
            data.ratio_close_phobic_to_close_philic = ratio_close_phobic_to_close_philic;
            data.ratio_far_philic_to_far_phobic = ratio_far_philic_to_far_phobic;
            data.ratio_close_phobic_to_far_phobic = ratio_close_phobic_to_far_phobic;
            data.ratio_far_philic_to_close_philic = ratio_far_philic_to_close_philic;
            data.ratio_n_close_phobic_to_n_phobic = ratio_n_close_phobic_to_n_phobic;
            data.ratio_n_close_phobic_to_total = ratio_n_close_phobic_to_total;
            data.limits.dist_min = dist_min;
            data.limits.dist_mid = dist_mid;
            data.limits.dist_max = dist_max;
            data.limits.hydropathy_min = hydropathy_min;
            data.limits.hydropathy_mid = hydropathy_mid;
            data.limits.hydropathy_max = hydropathy_max;
        end
        function data = Calculate_Helix_Hydropathy_Data_For_Scatter(self, helix_index)
            axis_distances = self.Calculate_Axis_Sidechain_Distances();
            if (helix_index > 5 || helix_index < 1) 
                print("ERROR: helix_index must be a valid helix number");
                return
            end
            
            if length(axis_distances) ~= 7
            
                hydropathy = cell(4,1);
                for i = 1:4
                    hydropathy{i} = self.helices(i).sidechain_properties(:,3);
                end
                comparison = cell(4,1);
                for i = 1:4
                    comparison{i} = horzcat(axis_distances{i},hydropathy{i});
                end

%                 results = vertcat(comparison{1}(:,1), comparison{2}(:,1), comparison{3}(:,1), comparison{4}(:,1));
%                 results = horzcat(results, vertcat(comparison{1}(:,2), comparison{2}(:,2), comparison{3}(:,2), comparison{4}(:,2)));
                
                results = horzcat(comparison{i}(:,1), comparison{i}(:,2));
                
                data = self.Calculate_Bin_Ratio_Hydro_Data(results);

               
            else
                index = self.data.three_grouped_helices;
                potential_index_values = [1,2,3,4,5];
                other_indeces = ismember(potential_index_values, index);
                other_index = potential_index_values(~other_indeces);
                managed_full_index = horzcat(index, other_index, index(2), index(3));
                index = managed_full_index;
                
                hydropathy = cell(7,1);
                for i = 1:7
                    hydropathy{i} = self.helices(index(i)).sidechain_properties(:,3);
                end
                comparison = cell(7,1);
                for i = 1:7
                    comparison{i} = horzcat(axis_distances{i},hydropathy{i});
                end
                
%                 results_3 = [];
%                 for i = 1:3
%                     results_3 = vertcat(results_3, comparison{i}(:,:));
%                 end
%                 
%                 results_4 = [];
%                 for i = 4:7
%                     results_4 = vertcat(results_4, comparison{i}(:,:));
%                 end
                
                index_3 =  managed_full_index == potential_index_values(helix_index);
                index_3(4:7) = 0;
                no_3 = all(index_3 == 0);
                if (not(no_3))
                    results_3 = horzcat(comparison{index_3}(:,1),comparison{index_3}(:,2));
                    data_3 = self.Calculate_Bin_Ratio_Hydro_Data(results_3);
                    data_3.is_three_part = 1;
                else
                    data_3 = [];
                end
                
                index_4 = managed_full_index == potential_index_values(helix_index);
                index_4(1:3) = 0;
                no_4 = all(index_4 == 0);
                
                if (not(no_4))
                    results_4 = horzcat(comparison{index_4}(:,1),comparison{index_4}(:,2));                 
                    data_4 = self.Calculate_Bin_Ratio_Hydro_Data(results_4);
                    data_4.is_three_part = 0;
                else
                    data_4 = [];
                end
                
                data = [data_3, data_4];
            end
        end
        function data = Calculate_Bundle_Hydropathy_Data(self)
            axis_distances = self.Calculate_Axis_Sidechain_Distances();
            
            if length(axis_distances) ~= 7
            
                hydropathy = cell(4,1);
                for i = 1:4
                    hydropathy{i} = self.helices(i).sidechain_properties(:,3);
                end
                comparison = cell(4,1);
                for i = 1:4
                    comparison{i} = horzcat(axis_distances{i},hydropathy{i});
                end

                results = vertcat(comparison{1}(:,1), comparison{2}(:,1), comparison{3}(:,1), comparison{4}(:,1));
                results = horzcat(results, vertcat(comparison{1}(:,2), comparison{2}(:,2), comparison{3}(:,2), comparison{4}(:,2)));
                
                data = self.Calculate_Bin_Ratio_Hydro_Data(results);

               
            else
                index = self.data.three_grouped_helices;
                potential_index_values = [1,2,3,4,5];
                other_indeces = ismember(potential_index_values, index);
                other_index = potential_index_values(~other_indeces);
                managed_full_index = horzcat(index, other_index, index(2), index(3));
                index = managed_full_index;
                
                hydropathy = cell(7,1);
                for i = 1:7
                    hydropathy{i} = self.helices(index(i)).sidechain_properties(:,3);
                end
                comparison = cell(7,1);
                for i = 1:7
                    comparison{i} = horzcat(axis_distances{i},hydropathy{i});
                end
                
                results_3 = [];
                for i = 1:3
                    results_3 = vertcat(results_3, comparison{i}(:,:));
                end
                
                results_4 = [];
                for i = 4:7
                    results_4 = vertcat(results_4, comparison{i}(:,:));
                end
                
                data_3 = self.Calculate_Bin_Ratio_Hydro_Data(results_3);
                data_4 = self.Calculate_Bin_Ratio_Hydro_Data(results_4);
                data = [data_3, data_4];
            end
        end
        function data = Calculate_Bundle_Hydropathy_Data_For_Helix(self, helix_index)
            axis_distances = self.Calculate_Axis_Sidechain_Distances();
            if (helix_index > 5 || helix_index < 1) 
                print("ERROR: helix_index must be a valid helix number");
                return
            end
            
            if length(axis_distances) ~= 7
            
                hydropathy = cell(4,1);
                for i = 1:4
                    hydropathy{i} = self.helices(i).sidechain_properties(:,3);
                end
                comparison = cell(4,1);
                for i = 1:4
                    comparison{i} = horzcat(axis_distances{i},hydropathy{i});
                end
                
                results = horzcat(comparison{i}(:,1), comparison{i}(:,2));
                
                data = self.Calculate_Bin_Ratio_Hydro_Data(results);
                data.results = results;
                data.helix_indeces = helix_index;

               
            else
                index = self.data.three_grouped_helices;
                potential_index_values = [1,2,3,4,5];
                other_indeces = ismember(potential_index_values, index);
                other_index = potential_index_values(~other_indeces);
                managed_full_index = horzcat(index, other_index, index(2), index(3));
                index = managed_full_index;
                
                hydropathy = cell(7,1);
                for i = 1:7
                    hydropathy{i} = self.helices(index(i)).sidechain_properties(:,3);
                end
                comparison = cell(7,1);
                for i = 1:7
                    comparison{i} = horzcat(axis_distances{i},hydropathy{i});
                end
                
                index_3 =  managed_full_index == potential_index_values(helix_index);
                index_3(4:7) = 0;
                no_3 = all(index_3 == 0);
                if (not(no_3))
                    results_3 = horzcat(comparison{index_3}(:,1),comparison{index_3}(:,2));
                    data_3 = self.Calculate_Bin_Ratio_Hydro_Data(results_3);
                    data_3.results = results_3;
                    data_3.is_three_part = 1;
                    data_3.helix_indeces = helix_index;
                else
                    data_3 = [];
                end
                
                index_4 = managed_full_index == potential_index_values(helix_index);
                index_4(1:3) = 0;
                no_4 = all(index_4 == 0);
                
                if (not(no_4))
                    results_4 = horzcat(comparison{index_4}(:,1),comparison{index_4}(:,2));                 
                    data_4 = self.Calculate_Bin_Ratio_Hydro_Data(results_4);
                    data_4.results = results_4;
                    data_4.is_three_part = 0;
                    data_4.helix_indeces = helix_index;
                else
                    data_4 = [];
                end
                
                data = [data_3, data_4];
            end
        end
        
        % Hydropathy Sidechain Scatter Plots
        function fig_handle = Make_Hydropathy_Sidechain_Scatter_Plot_Helix(~, data)
            fig_handle = figure();
            hold on
            scatter(data(:,1), data(:,2));
            xlabel("Distance from Bundle Axis, angstroms");
            ylabel("Hydropathy");
            hold off
        end
        function fig_handle = Make_Hydropathy_Sidechain_Scatter_Plot_3_Part(self)
            self.Calculate_Bundle_Axis();
            self.Get_Bundle_Axis();
            axis_distances = self.Calculate_Axis_Sidechain_Distances();
            
            three_part_indeces = self.data.three_grouped_helices;
            
            hydropathy = cell(3,1);
            for i = 1:3
                hydropathy{i} = self.helices(three_part_indeces(i)).sidechain_properties(:,3);
            end
            comparison = cell(3,1);
            for i = 1:3
                comparison{i} = horzcat(axis_distances{i},hydropathy{i});
            end
            
            fig_handle = figure();
            hold on
            scatter(comparison{1}(:,1),comparison{1}(:,2));
            scatter(comparison{2}(:,1),comparison{2}(:,2));
            scatter(comparison{3}(:,1),comparison{3}(:,2));
            xlabel("Distance from Bundle Axis, angstroms");
            ylabel("Hydropathy");
            hold off            
        end
        function fig_handle = Make_Hydropathy_Sidechain_Scatter_Plot_4_Part(self)
            self.Calculate_Bundle_Axis();
            self.Get_Bundle_Axis();
            axis_distances = self.Calculate_Axis_Sidechain_Distances();
            
            if length(axis_distances) == 4
                hydropathy = cell(4,1);
                for i = 1:4
                    hydropathy{i} = self.helices(i).sidechain_properties(:,3);
                end
                comparison = cell(4,1);
                for i = 1:4
                    comparison{i} = horzcat(axis_distances{i},hydropathy{i});
                end
                fig_handle = figure();
                hold on
                scatter(comparison{1}(:,1),comparison{1}(:,2));
                scatter(comparison{2}(:,1),comparison{2}(:,2));
                scatter(comparison{3}(:,1),comparison{3}(:,2));
                scatter(comparison{4}(:,1),comparison{4}(:,2));
                xlabel("Distance from Bundle Axis, angstroms");
                ylabel("Hydropathy");
                hold off                
            else
                index = self.data.three_grouped_helices;
                potential_index_values = [1,2,3,4,5];
                other_indeces = ismember(potential_index_values, index);
                other_index = potential_index_values(~other_indeces);
                managed_full_index = horzcat(index, other_index, index(2), index(3));
                index = managed_full_index;
                
                hydropathy = cell(7,1);
                for i = 1:7
                    hydropathy{i} = self.helices(index(i)).sidechain_properties(:,3);
                end
                comparison = cell(7,1);
                for i = 1:7
                    comparison{i} = horzcat(axis_distances{i},hydropathy{i});
                end
                
                fig_handle = figure();
                hold on
                scatter(comparison{4}(:,1),comparison{4}(:,2));
                scatter(comparison{5}(:,1),comparison{5}(:,2));
                scatter(comparison{6}(:,1),comparison{6}(:,2));
                scatter(comparison{7}(:,1),comparison{7}(:,2));
                xlabel("Distance from Bundle Axis, angstroms");
                ylabel("Hydropathy");
                hold off  
                
            end
            
            
        end
        function fig_1 = Apply_Labels_Scatter_Plot(~,fig_1, data_raw)
            data = data_raw.limits;
            temp = ylim(fig_1.CurrentAxes);
            data.hydropathy_max = temp(1);
            data.hydropathy_min = temp(2);
            temp = xlim(fig_1.CurrentAxes);
            data.dist_min = temp(1);
            data.dist_max = temp(2);
            vert_line_x = linspace(data.dist_mid, data.dist_mid, 2);
            vert_line_y = linspace(data.hydropathy_min, data.hydropathy_max, 2);
            horz_line_x = linspace(data.dist_min, data.dist_max, 2);
            horz_line_y = linspace(data.hydropathy_mid, data.hydropathy_mid, 2);

            data = data_raw.limits;
            far_phobic(1) = data.dist_mid + ( ( data.dist_max - data.dist_mid ) / 2 );
            far_philic(1) = data.dist_mid + ( ( data.dist_max - data.dist_mid ) / 2 );
            close_phobic(1) = data.dist_mid + ( ( data.dist_min - data.dist_mid ) / 2 );
            close_philic(1) = data.dist_mid + ( ( data.dist_min - data.dist_mid ) / 2 );
            far_phobic(2) = data.hydropathy_mid + ( ( data.hydropathy_max - data.hydropathy_mid ) / 2 );
            far_philic(2) = data.hydropathy_mid + ( ( data.hydropathy_min - data.hydropathy_mid ) / 2 );
            close_phobic(2) = data.hydropathy_mid + ( ( data.hydropathy_max - data.hydropathy_mid ) / 2 );
            close_philic(2) = data.hydropathy_mid + ( ( data.hydropathy_min - data.hydropathy_mid ) / 2 );

            hold(fig_1.CurrentAxes, 'on');
            plot(fig_1.CurrentAxes,vert_line_x, vert_line_y,'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');
            plot(fig_1.CurrentAxes,horz_line_x, horz_line_y,'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');

            text(fig_1.CurrentAxes, far_phobic(1), far_phobic(2), num2str(data_raw.n_far_hydrophobic), 'FontSize', 16);
            text(fig_1.CurrentAxes, far_philic(1), far_philic(2), num2str(data_raw.n_far_hydrophilic), 'FontSize', 16);
            text(fig_1.CurrentAxes, close_phobic(1), close_phobic(2), num2str(data_raw.n_close_hydrophobic), 'FontSize', 16);
            text(fig_1.CurrentAxes, close_philic(1), close_philic(2), num2str(data_raw.n_close_hydrophilic), 'FontSize', 16);

            hold(fig_1.CurrentAxes, 'off');
        end
        function fig_handle = Apply_Fitted_Curve_Scatter_Plot(~, fig_handle, data)
            [curve,~] = fit(data(:,1), data(:,2), 'poly2');
            hold(fig_handle.CurrentAxes, 'on');
            plot(curve);
            hold(fig_handle.CurrentAxes, 'off');
            drawnow;
        end
        function fig_handle = Make_Hydropathy_Scatter_Plot_Helix_Fitted(~, helixData, withLabels)
            if (nargin < 3)
                withLabels = true;
            end
            % Make Figure
            fig_handle = figure();
            fig_handle.CurrentAxes = axes();
            hold(fig_handle.CurrentAxes, 'on');
            
            % Apply Scatter Plot of hydopathy vs distance from core
            scatter(fig_handle.CurrentAxes, helixData.results(:,1), helixData.results(:,2));
            if (withLabels)
                xlabel(fig_handle.CurrentAxes, "Distance from Bundle Axis, angstroms");
                ylabel(fig_handle.CurrentAxes, "Hydropathy");
            end

            % Workout the visual division
            limitsData = helixData.limits;            
            vert_line_x = linspace(limitsData.dist_mid, limitsData.dist_mid, 2);
            vert_line_y = linspace(limitsData.hydropathy_min, limitsData.hydropathy_max, 2);
            horz_line_x = linspace(limitsData.dist_min, limitsData.dist_max, 2);
            horz_line_y = linspace(limitsData.hydropathy_mid, limitsData.hydropathy_mid, 2);
            
            vert_line_y(1) = vert_line_y(1) - (0.05 * abs(vert_line_y(2) - vert_line_y(1)));
            vert_line_y(2) = vert_line_y(2) + (0.05 * abs(vert_line_y(2) - vert_line_y(1)));
            
            horz_line_x(1) = horz_line_x(1) - (0.05 * abs(horz_line_x(2) - horz_line_x(1)));
            horz_line_x(2) = horz_line_x(2) + (0.05 * abs(horz_line_x(2) - horz_line_x(1)));
            
            % Workout placement of count values
            limitsData = helixData.limits;
            far_phobic_pos(1) = limitsData.dist_mid + ( ( limitsData.dist_max - limitsData.dist_mid ) / 2 );
            far_philic_pos(1) = limitsData.dist_mid + ( ( limitsData.dist_max - limitsData.dist_mid ) / 2 );
            close_phobic_pos(1) = limitsData.dist_mid + ( ( limitsData.dist_min - limitsData.dist_mid ) / 2 );
            close_philic_pos(1) = limitsData.dist_mid + ( ( limitsData.dist_min - limitsData.dist_mid ) / 2 );
            far_phobic_pos(2) = limitsData.hydropathy_mid + ( ( limitsData.hydropathy_max - limitsData.hydropathy_mid ) / 2 );
            far_philic_pos(2) = limitsData.hydropathy_mid + ( ( limitsData.hydropathy_min - limitsData.hydropathy_mid ) / 2 );
            close_phobic_pos(2) = limitsData.hydropathy_mid + ( ( limitsData.hydropathy_max - limitsData.hydropathy_mid ) / 2 );
            close_philic_pos(2) = limitsData.hydropathy_mid + ( ( limitsData.hydropathy_min - limitsData.hydropathy_mid ) / 2 );

            % Plot Division lines
            plot(fig_handle.CurrentAxes,vert_line_x, vert_line_y,'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');
            plot(fig_handle.CurrentAxes,horz_line_x, horz_line_y,'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');

            % Plot count values
            text(fig_handle.CurrentAxes, far_phobic_pos(1), far_phobic_pos(2), num2str(helixData.n_far_hydrophobic), 'FontSize', 16);
            text(fig_handle.CurrentAxes, far_philic_pos(1), far_philic_pos(2), num2str(helixData.n_far_hydrophilic), 'FontSize', 16);
            text(fig_handle.CurrentAxes, close_phobic_pos(1), close_phobic_pos(2), num2str(helixData.n_close_hydrophobic), 'FontSize', 16);
            text(fig_handle.CurrentAxes, close_philic_pos(1), close_philic_pos(2), num2str(helixData.n_close_hydrophilic), 'FontSize', 16);
            
            % Calculate Fit Curve and Plot
            [curve,~] = fit(helixData.results(:,1), helixData.results(:,2), 'poly2');
            plot(curve);
            
            % Fix axes limits
            xlim(fig_handle.CurrentAxes, [horz_line_x(1), horz_line_x(2)]);
            ylim(fig_handle.CurrentAxes, [vert_line_y(1), vert_line_y(2)]);
            
            hold(fig_handle.CurrentAxes, 'off');
        end
        
        % Simulation Set-ups
        function Load_PDB_Parallel_Chain(self)
            self.baseFrame.addChild(Joint_Object());
            for i = 1:length(self.helices)
                self.baseFrame.children.addChild(self.helices(i));
                self.helices(i).parent_frame.Set_Frame_to_Helix_Original();
            end
            translation_factor = self.helices(1).original_pos;
            for i = 1:length(self.helices)
                self.helices(i).parent_frame.Translate_Frame(-translation_factor);
                self.helices(i).parent_frame.insertJointBefore();
            end
        end
        
        % Sidechain Distance Calculations
        function output = Get_Distance_Between_Two_Sidechains(self, helix_one_index, helix_two_index, helix_one_sidechain_index, helix_two_sidechain_index)
            distance_unit = 1e-9;
            num_sc_one = self.helices(helix_one_index).Get_Number_Sidechains();
            num_sc_two = self.helices(helix_two_index).Get_Number_Sidechains();
            if helix_one_sidechain_index > num_sc_one
                error("Helix One sidechain index is greater than the number of sidechains");
            end
            if helix_two_sidechain_index > num_sc_two
                error("Helix One sidechain index is greater than the number of sidechains");
            end
            % pos_one = zeros(3,1);
            % pos_two = zeros(3,1);
            
            pos_one = self.helices(helix_one_index).parent_frame.Get_Position_World_Frame(...
                self.helices(helix_one_index).sidechain_positions(helix_one_sidechain_index,:)');
            pos_two = self.helices(helix_two_index).parent_frame.Get_Position_World_Frame(...
                self.helices(helix_two_index).sidechain_positions(helix_two_sidechain_index,:)');
            distance = vecnorm(pos_two - pos_one);
            output = distance * distance_unit;
        end
        function output = Get_Distance_Between_Sidechain_Vector(self, helix_one_index, helix_two_index, helix_one_sidechain_index)
            distance_unit = 1e-9;
            num_sc_one = self.helices(helix_one_index).Get_Number_Sidechains();
            num_sc_two = self.helices(helix_two_index).Get_Number_Sidechains();
            if helix_one_sidechain_index > num_sc_one
                error("Helix One sidechain index is greater than the number of sidechains");
            end
            % pos_one = zeros(3,1);
            % pos_two = zeros(3,num_sc_two);
            distances = zeros(1,num_sc_two);
            
            pos_one = self.helices(helix_one_index).Get_Position_World_Frame(...
                self.helices(helix_one_index).sidechain_positions(helix_one_sidechain_index,:)');
            pos_two = self.helices(helix_two_index).Get_Position_World_Frame(...
                self.helices(helix_two_index).sidechain_positions(:,:)');
            distances = vecnorm(minus(pos_two, pos_one));
            output = distances;% * distance_unit;
        end
        function output = Get_All_Sidechain_Distaces_Two_Helices(self, helix_one_index, helix_two_index)
            num_sc_one = self.helices(helix_one_index).Get_Number_Sidechains();
            num_sc_two = self.helices(helix_two_index).Get_Number_Sidechains();
            distances = zeros(num_sc_one,num_sc_two);
            for i = 1:num_sc_one
                distances(i,:) = self.Get_Distance_Between_Sidechain_Vector(helix_one_index,helix_two_index,i);
            end
            output = distances;
        end
        
        % Sidechain Radii Calculations
        function output = Get_Sidechain_Radius(self, helix_index, sidechain_index)
            output = self.helices(helix_index).sidechain_sizes(sidechain_index);
        end
        function output = Get_All_Sidechain_Radii_Helix(self, helix_index)
            output = self.helices(helix_index).sidechain_sizes;
        end
        
        % Sidechain Intersection Calculations
        function output = Run_Sidechain_Intersection_Two_Helices(self, helix_one_index, helix_two_index)
            sidechain_distances = self.Get_All_Sidechain_Distaces_Two_Helices(helix_one_index, helix_two_index);
            h1_sidechain_radii = self.Get_All_Sidechain_Radii_Helix(helix_one_index);
            h2_sidechain_radii = self.Get_All_Sidechain_Radii_Helix(helix_two_index);
            
            num_sc_1 = length(h1_sidechain_radii);
            num_sc_2 = length(h2_sidechain_radii);
            sidechain_radii_combined = zeros(num_sc_1 ,num_sc_2);
            for i = 1:num_sc_1
                sidechain_radii_combined(i,:) = h1_sidechain_radii(i) + h2_sidechain_radii';
            end
            
            sidechain_distance_radii_comparison = sidechain_distances - sidechain_radii_combined;
            output.comparison = sidechain_distance_radii_comparison;
            output.distances = sidechain_distances;
            output.radii_combined = sidechain_radii_combined;
        end
        function output = Run_Sidechain_Intersection(self)
            combinations = nchoosek(1:length(self.helices),2);
            for i = 1:length(combinations)
                results(i) = self.Run_Sidechain_Intersection_Two_Helices(combinations(i,1),combinations(i,2));
                sidechain_radii_comparison{i} = results(i).comparison;
                sidechain_distances{i} = results(i).distances;
                sidechain_radii_combined{i} = results(i).distances;
            end
            
            coords = cell(length(combinations),1);
            for i = 1:length(results)
                linear_indeces = find(sidechain_radii_comparison{i} < 0);
                coords{i} = zeros(length(linear_indeces),2);
                [rows,cols] = ind2sub(size(sidechain_radii_comparison{i}),find(sidechain_radii_comparison{i} < 0));
                coords{i} = horzcat(rows,cols);
            end
            
            collisions = [];
            for i = 1:length(sidechain_radii_comparison)
                if ~isempty(coords{i})
                    collision.helices = combinations(i,:);
                    collision.sidechains = coords{i};
                    collision.intersection = sidechain_radii_comparison{i}(coords{i}(1),coords{i}(2));
                    collision.distance = sidechain_distances{i}(coords{i}(1),coords{i}(2));
                    collision.radii_combined = sidechain_radii_combined{i}(coords{i}(1),coords{i}(2));
                    collisions = vertcat(collisions,collision);
                end
            end
            output = collisions;
        end
        
        % Interaction Force Simulations       
        function output = Run_Static_Force(self, helix_one_index, helix_two_index)
            force = 0;
            constant = 1 / (4 * pi * 8.854e-12);
            charge_unit = 1.6e-19;
            distance_unit = 1e-9;
            
            num_sc_one = self.helices(helix_one_index).Get_Number_Sidechains();
            num_sc_two = self.helices(helix_two_index).Get_Number_Sidechains();
            charge_one = 0;
            charge_two = zeros(num_sc_two, 1);
            pos_one = zeros(3,1);
            pos_two = zeros(3,num_sc_two);
            
            for i = 1:num_sc_one
                pos_one = self.helices(helix_one_index).parent_frame.Get_Position_World_Frame(...
                    self.helices(helix_one_index).sidechain_positions(i,:)');
                charge_one = self.helices(helix_one_index).sidechain_properties(i,2);
                % Vectorised helix two calcs
                pos_two = self.helices(helix_two_index).parent_frame.Get_Position_World_Frame(...
                    self.helices(helix_two_index).sidechain_positions(:,:)');
                charge_two = self.helices(helix_two_index).sidechain_properties(:,2);
                distance = vecnorm(minus(pos_two, pos_one));
%                 calculated_force = ( constant * ( charge_one*charge_unit  * charge_two*charge_unit ) ./ ( distance*distance_unit )' .^ 2);
                calculated_force = (constant * charge_one *charge_unit * charge_two * charge_unit) ./ ((distance * distance_unit)' .^ 2);
                force = force + sum(calculated_force);
            end
            self.data.force = force;
            output = force;
        end
        function output = Static_Force_Along_Length_Helices(self,helix_one_index, helix_two_index)
            constant = 1 / (4 * pi * 8.854e-12);
            charge_unit = 1.6e-19;
            distance_unit = 1e-9;
            
            num_sc_one = self.helices(helix_one_index).Get_Number_Sidechains();
            num_sc_two = self.helices(helix_two_index).Get_Number_Sidechains();
            
            forces_one = zeros(num_sc_one,1);
            forces_two = zeros(num_sc_two,1);
            
            charge_one = 0;
            charge_two = zeros(num_sc_two, 1);
            pos_one = zeros(3,1);
            pos_two = zeros(3,num_sc_two);
            
            for i = 1:num_sc_one
                pos_one = self.helices(helix_one_index).parent_frame.Get_Position_World_Frame(...
                    self.helices(helix_one_index).sidechain_positions(i,:)');
                charge_one = self.helices(helix_one_index).sidechain_properties(i,2);
                for j = 1:num_sc_two
                    pos_two = self.helices(helix_two_index).parent_frame.Get_Position_World_Frame(...
                        self.helices(helix_two_index).sidechain_positions(j,:)');
                    charge_two = self.helices(helix_two_index).sidechain_properties(j,2);
                    
                    distance = vecnorm(minus(pos_two,pos_one));
                    force = (constant * charge_one *charge_unit * charge_two * charge_unit) / ((distance * distance_unit)' ^2);
                    forces_one(i) = forces_one(i) + force;
                    forces_two(j) = forces_two(j) + force;
                end
            end
            output.forces_one = forces_one;
            output.forces_two = forces_two;
        end
        function [forceSum, forces, combinations] = Run_Static_Force_On_Specified_Bundle(self, bundle_pdb_path, bundle_helix_indeces)
            self.Load_Amino_Acid_PDB_File(bundle_pdb_path);
            num_helices = length(bundle_helix_indeces(:,:,1));
            if (bundle_helix_indeces(:,5,1) == 0)
                num_helices = num_helices - 1;
            end
            for i = 1:num_helices
                self.Create_Helix(bundle_helix_indeces(1,i,1),bundle_helix_indeces(1,i,2));
            end
            self.Load_PDB_Parallel_Chain();
            combinations = nchoosek(1:num_helices,2);
            forces = zeros(length(num_helices), 1);
            for i = 1:length(combinations)
                forces(i) = self.Run_Static_Force(combinations(i,1), combinations(i,2));
            end
            forceSum = sum(forces);
        end
        function [forceSum, forces, combinations] = Run_Bundle_Static_Force(self)
            num_helices = length(self.helices);
            combinations = nchoosek(1:num_helices,2);
            forces = zeros(num_helices,1);
            for i = 1:length(combinations)
                forces(i) = self.Run_Static_Force(combinations(i,1), combinations(i,2));
            end
            forceSum = sum(forces);
        end
        function output = Run_Interaction_Mapping(self, helix_index_one, helix_index_two, angles)
            if (nargin < 4)
                angles = 1:10:360;
            end
            num_angles = length(angles);
            force = zeros(num_angles, num_angles);
            constant = 1 / (4 * pi * 8.854e-12);
            charge_unit = 1.6e-19;
            distance_unit = 1e-9;
            
            num_sc_one = self.helices(helix_index_one).Get_Number_Sidechains();
            num_sc_two = self.helices(helix_index_two).Get_Number_Sidechains();
            charge_one = 0;
            charge_two = zeros(num_sc_two, 1);
            pos_one = zeros(3,1);
            pos_two = zeros(3,num_sc_two);
            
            for a = 1:num_angles
                self.helices(helix_index_one).parent_frame.Set_Rotation_AA([1,0,0], angles(a));
                for b = 1:num_angles
                    self.helices(helix_index_two).parent_frame.Set_Rotation_AA([1,0,0], angles(b));
                    for i = 1:num_sc_one
                        pos_one = self.helices(helix_index_one).parent_frame.Get_Position_World_Frame(...
                            self.helices(helix_index_one).sidechain_positions(i,:)');
                        charge_one = self.helices(helix_index_one).sidechain_properties(i,2);
                        pos_two = self.helices(helix_index_two).parent_frame.Get_Position_World_Frame(...
                            self.helices(helix_index_two).sidechain_positions(:,:)');
                        charge_two = self.helices(helix_index_two).sidechain_properties(:,2);
                        distance = vecnorm(minus(pos_two, pos_one));
                        calculated_force = ( constant * ( charge_one*charge_unit * charge_two*charge_unit ) ./ ( distance*distance_unit )'.^2);
                        force(a,b) = force(a,b) + sum(calculated_force);
                    end
                end
            end
            self.data.force = force;
            output = force;
            
        end
        function fig_handle = Generate_Mapping_Figure(self,forces)
            switch (nargin)
                case 1
                    % No forces included, therefore used stored values
                    fig_handle = figure();
                    surf(self.data, 'EdgeAlpha', 0.25);
                    colormap(jet);
                case 2
                    fig_handle = figure();
                    surf(forces, 'EdgeAlpha', 0.25);
                    colormap(jet);
                otherwise
                    error("Invalid Number of Input Arguements");
            end
        end
        
        % --Dev-- Unbounded Optimsation Functions
        function output = Run_Wiggle_Translation(self, moving_helix_index, static_helix_index)
            length_helix = ceil(self.helices(moving_helix_index).position_Two(1));
            start_distance = length_helix / 4;
            distance_length = (length_helix / 2) + 1;
            saved_position = self.helices(moving_helix_index).parent_frame.translation_from_parent;
            forces = zeros(distance_length, 1);
            
            self.helices(moving_helix_index).parent_frame.Translate_Frame([start_distance, 0, 0]);
            forces(1) = self.Run_Static_Force(moving_helix_index, static_helix_index);
            for i = 2:distance_length
                self.helices(moving_helix_index).parent_frame.Translate_Frame([-1,0,0]);
                forces(i) = self.Run_Static_Force(moving_helix_index, static_helix_index);
            end
            
            self.helices(moving_helix_index).parent_frame.Set_Position_Frame(saved_position);
            
            output.axis = start_distance:-1:-start_distance;
            output.forces = forces;
        end
        function output = Run_Wiggle_Rotation(self, moving_helix_index, static_helix_index)
            start_angle = 30.0; % Degrees
            angle_length = (2 * start_angle) + 1.0;
            saved_rotation = self.helices(moving_helix_index).parent_frame.rotation_from_parent;
            forces = zeros(angle_length, 1);
            
            self.helices(moving_helix_index).parent_frame.Rotate_Frame_AA([1,0,0], start_angle);
            forces(1) = self.Run_Static_Force(moving_helix_index, static_helix_index);
            for i = 2:angle_length
                self.helices(moving_helix_index).parent_frame.Rotate_Frame_AA([1,0,0], -1);
                forces(i) = self.Run_Static_Force(moving_helix_index, static_helix_index);
            end
            
            self.helices(moving_helix_index).parent_frame.Set_Rotation_RM(saved_rotation);
            
            output.axis = start_angle:-1:-start_angle;
            output.forces = forces;
        end
        function output = Run_Test_Optimisation_Rotation(self, angle)
            self.helices(1).parent_frame.Set_Rotation_AA([1,0,0],angle);
            output = self.Run_Static_Force(1,2);
        end
        function output = Run_Test_Optimisation_Translation(self, position)
            self.helices(1).parent_frame.Set_Position_Frame([position;0;0]);
            output = self.Run_Static_Force(1,2);
        end
        function output = Run_Test_Translation_Optimisation(self,offset)
            self.helices(2).parent_frame.Set_Position_Axial_By_Offset(self,offset);
            output = self.Run_Static_Force(1,2);
        end
        function output = Run_Bundle_Optimisation_Rotation(self, angles)
            self.helices(1).parent_frame.Set_Rotation_AA([1,0,0],angles(1));
            self.helices(2).parent_frame.Set_Rotation_AA([1,0,0],angles(2));
            self.helices(3).parent_frame.Set_Rotation_AA([1,0,0],angles(3));
            self.helices(4).parent_frame.Set_Rotation_AA([1,0,0],angles(4));
            combinations = nchoosek(1:4,2);
            if (length(angles) > 4)
                self.helices(5).parent_frame.Set_Rotation_AA([1,0,0],angles(5));
                combinations = nchoosek(1:5,2);
            end
            forces = zeros(length(combinations),1);
            for i = 1:length(combinations)
                forces(i) = self.Run_Static_Force(combinations(i,1), combinations(i,2));
            end
            sum_forces = sum(forces);
            output = sum_forces;
        end
        
        % Bounded Optimisation Functions
        function output = Bounded_Rotation_Optimisation(self, angles)
            num_helices = length(angles);
            for i = 1:num_helices
                self.helices(i).parent_frame.Set_Rotation_AA([1,0,0],angles(i));
            end
            if num_helices > 2
                combinations = nchoosek(1:num_helices,2);
                forces = zeros(length(combinations(:,1)),1);
                for i = 1:length(forces)
                    forces(i) = self.Run_Static_Force(combinations(i,1),combinations(i,2));
                end
                sum_forces = sum(forces);
                output = sum_forces;
            else
                output = self.Run_Static_Force(1,2);
            end
            
        end
        function output = Bounded_Translation_Optimisation(self, positions)
            num_helices = length(positions);
            for i = 1:num_helices
                self.helices(i).parent_frame.Set_Position_Axial_By_Offset(positions(i));
            end
            if num_helices > 2
                combinations = nchoosek(1:num_helices,2);
                forces = zeros(length(combinations(:,1)),1);
                for i = 1:length(forces)
                    forces(i) = self.Run_Static_Force(combinations(i,1),combinations(i,2));
                end
                sum_forces = sum(forces);
                output = sum_forces;
            else
                output = self.Run_Static_Force(1,2);
            end
        end
        function output = Bounded_RotationTranslation_Optimisation(self, angles_and_positions)
        % BOUNDED_ROTATIONTRANSLATION_OPTIMISATION  Optimises the
        % alpha-helix bundle via rotaion and translation of alpha-helices
        %
        %   OUTPUT =
        %   Bounded_RotationTranslation_Optimisation(ANGLES_AND_POSITIONS)
        %   the angles_and_positions vector arguement variabe contains the
        %   specified angles and positions to alter the alpha-helices by as
        %   dictated by the GlobalSearch or MultiStart optimisation
        %   functions. The first half the vector contains the angles
        %   information and the second half of the variabel conatins the
        %   position data.
        
            num_helices = length(angles_and_positions)/2;
            angles = angles_and_positions(1:(num_helices));
            positions = angles_and_positions(num_helices+1:2*num_helices);
            for i = 1:num_helices
                self.helices(i).parent_frame.Set_Rotation_AA([1,0,0],angles(i));
                self.helices(i).parent_frame.Set_Position_Axial_By_Offset(positions(i));
            end
            if num_helices > 2
                combinations = nchoosek(1:num_helices,2);
                forces = zeros(length(combinations(:,1)),1);
                for i = 1:length(forces)
                    forces(i) = self.Run_Static_Force(combinations(i,1),combinations(i,2));
                end
                sum_forces = sum(forces);
                output = sum_forces;
            else
                output = self.Run_Static_Force(1,2);
            end
        end
        function output = Bounded_RotationTranslation_Collision(self, angles_and_positions)
        % BOUNDED_ROTATIONTRANSLATION_COLLISION Optimises the
        % alpha-helix bundle via rotaion and translation of alpha-helices
        %
        %   OUTPUT =
        %   Bounded_RotationTranslation_Optimisation(ANGLES_AND_POSITIONS)
        %   the angles_and_positions vector arguement variabe contains the
        %   specified angles and positions to alter the alpha-helices by as
        %   dictated by the GlobalSearch or MultiStart optimisation
        %   functions. The first half the vector contains the angles
        %   information and the second half of the variabel conatins the
        %   position data.
        
            num_helices = length(angles_and_positions)/2;
            angles = angles_and_positions(1:(num_helices));
            positions = angles_and_positions(num_helices+1:2*num_helices);
            for i = 1:num_helices
                self.helices(i).parent_frame.Set_Rotation_AA([1,0,0],angles(i));
                self.helices(i).parent_frame.Set_Position_Axial_By_Offset(positions(i));
            end
            if num_helices > 2
                combinations = nchoosek(1:num_helices,2);
                forces = zeros(length(combinations(:,1)),1);
                for i = 1:length(forces)
                    forces(i) = self.Run_Static_Force(combinations(i,1),combinations(i,2));
                end
                sum_forces = sum(forces);
                output = sum_forces;
            else
                output = self.Run_Static_Force(1,2);
            end
            
            intersection_force_coefficient = 1 * 1e-6;
            collisions = self.Run_Sidechain_Intersection();
            if(~isempty(collisions))
                intersections = [collisions.intersection]';
                intersections = abs(intersections * intersection_force_coefficient);
                output = output + sum(intersections);
            end
            
        end
        
        % nonlcon Bounds Functions
        function [c,ceq] = Nonlcon_Maximum_Output_Limit(self, angles_and_positions)
            ceq = [];
            output = self.Bounded_RotationTranslation_Optimisation(angles_and_positions);
            if (output < (-45 * 1e-12))
                c = 1;
            else
                c = [];
            end            
        end
        
        % GlobalSearchh Wrapper Function
        function [x,fval] = GS_Rotation_Optimisation(self, num_helices, x0_coefficient, bounds_coefficient)
            % Setting num_helices to a lower number than the number of
            % helices will result in only the first number of helices
            % denoted by 'num_helices' to be rotated as part of 
            rng default;
            fun = @self.Bounded_Rotation_Optimisation;
            x0 = ones(1,num_helices) * x0_coefficient;
            bounds = ones(1,num_helices) * bounds_coefficient;
            problem = createOptimProblem('fmincon','objective',fun,...
                'x0',       x0,...
                'lb',       -bounds,...
                'ub',       bounds,...
                'options',  optimoptions(@fmincon));
            gs = GlobalSearch('Display','iter','PlotFcn',@gsplotbestf);
            [x,fval] = run(gs,problem);
        end
        function [x,fval] = GS_Translation_OPtimisation(self, num_helices, x0_coefficient, bounds_coefficient)
            % Setting num_helices to a lower number than the number of
            % helices will result in only the first number of helices
            % denoted by 'num_helices' to be rotated as part of 
            rng default;
            fun = @self.Bounded_Translation_Optimisation;
            x0 = ones(1,num_helices) * x0_coefficient;
            bounds = ones(1,num_helices) * bounds_coefficient;
            problem = createOptimProblem('fmincon','objective',fun,...
                'x0',       x0,...
                'lb',       -bounds,...
                'ub',       bounds,...
                'options',  optimoptions(@fmincon));
            ms = MultiStart('UseParallel',true,'Display','iter');
            myPool = parpool();
            [x,fval,eflag,output,manymins] = run(ms,problem,num_starting_points);
            delete(myPool);
        end
        % MultiStart Wrapper Function
        function [x,fval,eflag,output,manymins] = MS_Rotation_Optimisation(self, num_helices, x0_coefficient, bounds_coefficient, num_starting_points)
            % Setting num_helices to a lower number than the number of
            % helices will result in only the first number of helices
            % denoted by 'num_helices' to be rotated as part of 
            rng default;
            fun = @self.Bounded_Rotation_Optimisation;
            x0 = ones(1,num_helices) * x0_coefficient;
            bounds = ones(1,num_helices) * bounds_coefficient;
            problem = createOptimProblem('fmincon','objective',fun,...
                'x0',       x0,...
                'lb',       -bounds,...
                'ub',       bounds,...
                'options',  optimoptions(@fmincon));
            ms = MultiStart('UseParallel',true,'Display','iter');
            myPool = parpool();
            [x,fval,eflag,output,manymins] = run(ms,problem,num_starting_points);
            delete(myPool);
        end
        function [x,fval,eflag,output,manymins] = MS_Rotation_Optimisation_Selective(self, num_helices, x0_coefficient, bounds_coefficient, num_starting_points, selection)
            % Setting num_helices to a lower number than the number of
            % helices will result in only the first number of helices
            % denoted by 'num_helices' to be rotated as part of 
            rng default;
            fun = @self.Bounded_Rotation_Optimisation;
            x0 = ones(1,num_helices) * x0_coefficient;
            bounds = ones(1,num_helices) * bounds_coefficient;
            bounds = bounds.*selection;
            problem = createOptimProblem('fmincon','objective',fun,...
                'x0',       x0,...
                'lb',       -bounds,...
                'ub',       bounds,...
                'options',  optimoptions(@fmincon));
            ms = MultiStart('UseParallel',true,'Display','off');
            [x,fval,eflag,output,manymins] = run(ms,problem,num_starting_points);
        end
        function [x,fval,eflag,output,manymins] = MS_Translation_Optimisation(self, num_helices, x0_coefficient, bounds_coefficient, num_starting_points)
            % Setting num_helices to a lower number than the number of
            % helices will result in only the first number of helices
            % denoted by 'num_helices' to be rotated as part of 
            rng default;
            fun = @self.Bounded_Translation_Optimisation;
            x0 = ones(1,num_helices) * x0_coefficient;
            bounds = ones(1,num_helices) * bounds_coefficient;
            problem = createOptimProblem('fmincon','objective',fun,...
                'x0',       x0,...
                'lb',       -bounds,...
                'ub',       bounds,...
                'options',  optimoptions(@fmincon));
            ms = MultiStart('UseParallel',true,'Display','iter');
            myPool = parpool();
            [x,fval,eflag,output,manymins] = run(ms,problem,num_starting_points);
            delete(myPool);
        end
        function [x,fval,eflag,output,manymins] = MS_Translation_Optimisation_Selective(self, num_helices, x0_coefficient, bounds_coefficient, num_starting_points, selection)
            % Setting num_helices to a lower number than the number of
            % helices will result in only the first number of helices
            % denoted by 'num_helices' to be rotated as part of 
            rng default;
            fun = @self.Bounded_Translation_Optimisation;
            x0 = ones(1,num_helices) * x0_coefficient;
            bounds = ones(1,num_helices) * bounds_coefficient;
            bounds = bounds.*selection;
            problem = createOptimProblem('fmincon','objective',fun,...
                'x0',       x0,...
                'lb',       -bounds,...
                'ub',       bounds,...
                'options',  optimoptions(@fmincon));
            ms = MultiStart('UseParallel',true,'Display','off');
            [x,fval,eflag,output,manymins] = run(ms,problem,num_starting_points);
        end
        function [x,fval,eflag,output,manymins] = MS_RotationTranslation_Optimisation_Selective(self, num_helices, x0_rotation_coefficient, x0_translation_coefficient, bounds_rotation_coefficient, bounds_traslation_coefficient, num_starting_points, selection)
            % Setting num_helices to a lower number than the number of
            % helices will result in only the first number of helices
            % denoted by 'num_helices' to be rotated as part of 
            rng default;
            fun = @self.Bounded_RotationTranslation_Optimisation;
            x0 = [ones(1,num_helices) * x0_rotation_coefficient,ones(1,num_helices) * x0_translation_coefficient];
            bounds = [ones(1,num_helices) * bounds_rotation_coefficient,ones(1,num_helices) * bounds_traslation_coefficient];
            bounds = bounds.*[selection,selection];
            problem = createOptimProblem('fmincon','objective',fun,...
                'x0',       x0,...
                'lb',       -bounds,...
                'ub',       bounds,...
                'options',  optimoptions(@fmincon));
            ms = MultiStart('UseParallel',true,'Display','off');
            [x,fval,eflag,output,manymins] = run(ms,problem,num_starting_points);
        end
        function [x,fval,eflag,output,manymins] = MS_RotationTranslation_Collision(self, num_helices, x0_rotation_coefficient, x0_translation_coefficient, bounds_rotation_coefficient, bounds_traslation_coefficient, num_starting_points, selection)
            % Setting num_helices to a lower number than the number of
            % helices will result in only the first number of helices
            % denoted by 'num_helices' to be rotated as part of 
            rng default;
            fun = @self.Bounded_RotationTranslation_Collision;
            nonlcon = @self.Nonlcon_Maximum_Output_Limit;
            x0 = [ones(1,num_helices) * x0_rotation_coefficient,ones(1,num_helices) * x0_translation_coefficient];
            bounds = [ones(1,num_helices) * bounds_rotation_coefficient,ones(1,num_helices) * bounds_traslation_coefficient];
            bounds = bounds.*[selection,selection];
            problem = createOptimProblem('fmincon','objective',fun,...
                'x0',       x0,...
                'lb',       -bounds,...
                'ub',       bounds,...
                'options',  optimoptions(@fmincon),...
                'nonlcon',  nonlcon);
            ms = MultiStart('UseParallel',true,'Display','off');
            [x,fval,eflag,output,manymins] = run(ms,problem,num_starting_points);
        end
        
        % Force Optimisation Simulations
        function myPool = Initialise_Optim_Rot_Sim(self,domain_index,helix_indeces,helices_to_include,helices_to_rotate)
            self.data.optim.helices_to_include = helices_to_include;
            self.data.optim.helices_to_rotate = helices_to_rotate;
            self.data.optim.domain_index = domain_index;
            
            filepaths = Get_Talin_PDB_Filepaths();
            domain_filepath = filepaths(domain_index);
            
            self.Load_Amino_Acid_PDB_File(domain_filepath);
            self.Load_Helices_Order(domain_index,helices_to_include, helix_indeces);
            self.Load_PDB_Parallel_Chain();
            
            % Check if pool exists, if not create new one.
            exists = gcp('nocreate');
            if(isempty(exists))
                myPool = parpool();
            end
            return;
        end
        function results = Run_Optim_Rot_Sim(self, inital_conditions_coefficient, bounds_coefficient, num_starting_points)
           if (~isfield(self.data,'optim')); return; end
           
           self.data.optim.map_logical = ismember(self.data.optim.helices_to_include,self.data.optim.helices_to_rotate);
           self.data.optim.map_indeces = find(self.data.optim.map_logical);
           num_helices = length(self.data.optim.helices_to_include);
                      
           [x,fval,eflag,output,manymins] = self.MS_Rotation_Optimisation_Selective(num_helices, inital_conditions_coefficient, bounds_coefficient,num_starting_points,self.data.optim.map_logical);
           % Update Simulation helices with final rotations from optim sim
           for index = 1:length(x)
               self.helices(index).parent_frame.Set_Rotation_AA([1,0,0],x(index));
           end
           
            results.domain_index = self.data.optim.domain_index;
            results.helices_included = self.data.optim.helices_to_include;
            results.helices_rotated = self.data.optim.helices_to_rotate;
            results.inital_conditions_coefficient = inital_conditions_coefficient;
            results.bounds_coefficient = bounds_coefficient;
            results.final_force = fval;
            results.final_configuration = x;
        end
        
        function myPool = Initialise_Optim_Trans_Sim(self,domain_index,helix_indeces,helices_to_include,helices_to_rotate)
            self.data.optim.helices_to_include = helices_to_include;
            self.data.optim.helices_to_rotate = helices_to_rotate;
            self.data.optim.domain_index = domain_index;
            
            filepaths = Get_Talin_PDB_Filepaths();
            domain_filepath = filepaths(domain_index);
            
            self.Load_Amino_Acid_PDB_File(domain_filepath);
            self.Load_Helices_Order(domain_index,helices_to_include, helix_indeces);
            self.Load_PDB_Parallel_Chain();
            
            % Check if pool exists, if not create new one.
            exists = gcp('nocreate');
            if(isempty(exists))
                myPool = parpool();
            end
            return;
        end
        function results = Run_Optim_Trans_Sim(self, inital_conditions_coefficient, bounds_coefficient, num_starting_points)
           if (~isfield(self.data,'optim')); return; end
           
           self.data.optim.map_logical = ismember(self.data.optim.helices_to_include,self.data.optim.helices_to_rotate);
           self.data.optim.map_indeces = find(self.data.optim.map_logical);
           num_helices = length(self.data.optim.helices_to_include);
                      
           [x,fval,eflag,output,manymins] = self.MS_Translation_Optimisation_Selective(num_helices, inital_conditions_coefficient, bounds_coefficient,num_starting_points,self.data.optim.map_logical);
           % Update Simulation helices with final rotations from optim sim
           for index = 1:length(x)
%                self.helices(index).parent_frame.Set_Position_Axial_By_Offset([1,0,0],x(index));
                self.helices(index).parent_frame.Set_Position_Axial_By_Offset(x(index));
           end
           
            results.domain_index = self.data.optim.domain_index;
            results.helices_included = self.data.optim.helices_to_include;
            results.helices_rotated = self.data.optim.helices_to_rotate;
            results.inital_conditions_coefficient = inital_conditions_coefficient;
            results.bounds_coefficient = bounds_coefficient;
            results.final_force = fval;
            results.final_configuration = x;
        end
        
        function myPool = Initialise_Optim_RotTrans_Sim(self,domain_index,helix_indeces,helices_to_include,helices_to_rotate)
            self.data.optim.helices_to_include = helices_to_include;
            self.data.optim.helices_to_rotate = helices_to_rotate;
            self.data.optim.domain_index = domain_index;
            
            filepaths = Get_Talin_PDB_Filepaths();
            domain_filepath = filepaths(domain_index);
            
            self.Load_Amino_Acid_PDB_File(domain_filepath);
            self.Load_Helices_Order(domain_index,helices_to_include, helix_indeces);
            self.Load_PDB_Parallel_Chain();
            
            % Check if pool exists, if not create new one.
            exists = gcp('nocreate');
            if(isempty(exists))
                myPool = parpool();
            end
            return;
        end
        function results = Run_Optim_RotTrans_Sim(self, inital_conditions_rotation_coefficient, inital_conditions_translation_coefficient, bounds_rotation_coefficient, bounds_translation_coefficient, num_starting_points)
           if (~isfield(self.data,'optim')); return; end
           
           self.data.optim.map_logical = ismember(self.data.optim.helices_to_include,self.data.optim.helices_to_rotate);
           self.data.optim.map_indeces = find(self.data.optim.map_logical);
           num_helices = length(self.data.optim.helices_to_include);
                      
%            [x,fval,eflag,output,manymins] = self.MS_Translation_Optimisation_Selective(num_helices, inital_conditions_coefficient, bounds_coefficient,num_starting_points,self.data.optim.map_logical);
           [x,fval,~,~,~] = self.MS_RotationTranslation_Optimisation_Selective(num_helices, inital_conditions_rotation_coefficient, inital_conditions_translation_coefficient, bounds_rotation_coefficient, bounds_translation_coefficient,num_starting_points,self.data.optim.map_logical);
           % Update Simulation helices with final rotations from optim sim
           for index = 1:length(x)/2
               self.helices(index).parent_frame.Set_Rotation_AA([1,0,0],x(index));
           end
           for index = (length(x)/2)+1:length(x)
               self.helices(index-(length(x)/2)).parent_frame.Set_Position_Axial_By_Offset(x(index));
           end
           
            results.domain_index = self.data.optim.domain_index;
            results.helices_included = self.data.optim.helices_to_include;
            results.helices_rotated = self.data.optim.helices_to_rotate;
            results.inital_conditions_coefficient = [inital_conditions_rotation_coefficient,inital_conditions_translation_coefficient];
            results.bounds_coefficient = [bounds_rotation_coefficient,bounds_translation_coefficient];
            results.final_force = fval;
            results.final_rotation_configuration = x(1:length(x)/2);
            results.final_translation_configuration = x(length(x)/2+1:length(x));
        end
        
        function myPool = Initialise_Optim_RotTrans_Collision_Sim(self,domain_index,helix_indeces,helices_to_include,helices_to_rotate)
            self.data.optim.helices_to_include = helices_to_include;
            self.data.optim.helices_to_rotate = helices_to_rotate;
            self.data.optim.domain_index = domain_index;
            
            filepaths = Get_Talin_PDB_Filepaths();
            domain_filepath = filepaths(domain_index);
            
            self.Load_Amino_Acid_PDB_File(domain_filepath);
            self.Load_Helices_Order(domain_index,helices_to_include, helix_indeces);
            self.Load_PDB_Parallel_Chain();
            
            % Check if pool exists, if not create new one.
            exists = gcp('nocreate');
            if(isempty(exists))
                myPool = parpool();
            end
            return;
        end
        function myPool = Initialise_Custom_Optim_RotTrans_Collision_Sim(self,domain_index,helices_to_include,helices_to_adjust)
            self.data.optim.helices_to_include = helices_to_include;
            self.data.optim.helices_to_rotate = helices_to_adjust;
            self.data.optim.domain_index = domain_index;
            
            exists = gcp('nocreate');
            if(isempty(exists))
                myPool = parpool();
            end
            return;
        end
        function results = Run_Optim_RotTrans_Collision_Sim(self, inital_conditions_rotation_coefficient, inital_conditions_translation_coefficient, bounds_rotation_coefficient, bounds_translation_coefficient, num_starting_points)
           if (~isfield(self.data,'optim')); return; end
           
           self.data.optim.map_logical = ismember(self.data.optim.helices_to_include,self.data.optim.helices_to_rotate);
           self.data.optim.map_indeces = find(self.data.optim.map_logical);
           num_helices = length(self.data.optim.helices_to_include);
                      
%            [x,fval,eflag,output,manymins] = self.MS_Translation_Optimisation_Selective(num_helices, inital_conditions_coefficient, bounds_coefficient,num_starting_points,self.data.optim.map_logical);
           [x,fval,~,~,~] = self.MS_RotationTranslation_Collision(num_helices, inital_conditions_rotation_coefficient, inital_conditions_translation_coefficient, bounds_rotation_coefficient, bounds_translation_coefficient,num_starting_points,self.data.optim.map_logical);
           % Update Simulation helices with final rotations from optim sim
           for index = 1:length(x)/2
               self.helices(index).parent_frame.Set_Rotation_AA([1,0,0],x(index));
           end
           for index = (length(x)/2)+1:length(x)
               self.helices(index-(length(x)/2)).parent_frame.Set_Position_Axial_By_Offset(x(index));
           end
           
            results.domain_index = self.data.optim.domain_index;
            results.helices_included = self.data.optim.helices_to_include;
            results.helices_rotated = self.data.optim.helices_to_rotate;
            results.inital_conditions_coefficient = [inital_conditions_rotation_coefficient,inital_conditions_translation_coefficient];
            results.bounds_coefficient = [bounds_rotation_coefficient,bounds_translation_coefficient];
            results.final_force = fval;
            results.final_rotation_configuration = x(1:length(x)/2);
            results.final_translation_configuration = x(length(x)/2+1:length(x));
        end
        
    end
    methods (Static)
        % Helper Functions
        function attributes = Generate_Attributes_Struct()
            attributes.mode = "";
            attributes.helix_colours = [[-1,-1,-1];[-1,-1,-1];[-1,-1,-1];[-1,-1,-1];[-1,-1,-1]];
            attributes.primary_sidechain_stick_colour = '';
        end
        function pdb_filepaths = Load_PDB_Filepaths()
            pdb_filepaths = [ ...
                "pdb_data/R1.pdb",...
                "pdb_data/R2.pdb",...
                "pdb_data/R3.pdb",...
                "pdb_data/R4.pdb",...
                "pdb_data/R5.pdb",...
                "pdb_data/R6.pdb",...
                "pdb_data/R7R8.pdb",...
                "pdb_data/R7R8.pdb",...
                "pdb_data/R9.pdb",...
                "pdb_data/R10.pdb",...
                "pdb_data/R11.pdb",...
                "pdb_data/R12.pdb",...
                "pdb_data/R13.pdb",...
                ];
        end
        function sim = Load_Sim_Helix(domain_index, helix_index)
            pdb_filepaths = Simulation.Load_PDB_Filepaths();
            helix_indeces = load("helix_indeces.mat");
            helix_indeces = helix_indeces.helix_indeces;
            
            sim = Simulation();
            sim.Load_Amino_Acid_PDB_File(pdb_filepaths(domain_index));
            sim.Create_Helix(helix_indeces(domain_index, helix_index, 1),helix_indeces(domain_index, helix_index, 2));
        end
        function sim = Load_Sim_Talin_Subdomain(index)
            pdb_filepaths = Simulation.Load_PDB_Filepaths();
            helix_indeces = load("helix_indeces.mat");
            helix_indeces = helix_indeces.helix_indeces;
            
            sim = Simulation();
            sim.Load_Amino_Acid_PDB_File(pdb_filepaths(index));
            
            helix_count = length(helix_indeces(index,:,1));
            if (helix_indeces(index,5,1) == 0)
                helix_count = helix_count - 1;
            end
            
            for h = 1:helix_count
                sim.Create_Helix(helix_indeces(index,h,1),helix_indeces(index, h, 2));
            end
            sim.Load_PDB_Parallel_Chain();
        end
        
        % Run Pre-Configured Simulations
        function forces = Run_All_Bundle_Static_Forces()
            pdb_filepaths = Simulation.Load_PDB_Filepaths();
            helix_indeces = load("helix_indeces.mat");
            helix_indeces = helix_indeces.helix_indeces;
            
            for i = 1:13
               sim = Simulation();
               [sum, ~, ~] = sim.Run_Static_Force_On_Specified_Bundle(pdb_filepaths(i), helix_indeces(i,:,:));
               forces(i) = sum;
            end
        end
        % All-in-One Force Optimisation Sim
%         function results = Optim_Simulation
        function results = Optim_Simulation_Bundle_Bounds(domain_index,Rot_Bounds_Max, Trans_Bounds_Max, Rot_Bounds_Delta, Trans_Bounds_Delta)
            pdb_filepaths = Simulation.Load_PDB_Filepaths();
            helix_indeces = load("helix_indeces.mat");
            helix_indeces = helix_indeces.helix_indeces;

            rot_index_vector = round(linspace(1,Rot_Bounds_Max, 20));
            trans_index_vector = round(linspace(1,Trans_Bounds_Max, 10));

            rot_index_vector_size = length(rot_index_vector);
            trans_index_vector_size = length(trans_index_vector);

            num_total_results = rot_index_vector_size*trans_index_vector_size*13;
            
            domain_results = [];
            for rot_index = 1:length(rot_index_vector)
                fprintf("Domain: %d, Round: (%d/%d)\n", domain_index, rot_index,length(rot_index_vector));
                for trans_index = 1:length(trans_index_vector)
                    % Determine Number of Helices in domain.
                    num_helices = 5;
                    if (helix_indeces(domain_index,5,1) == 0)
                        num_helices = 4;
                    end

                    helices_to_include = 1:num_helices;
                    helices_to_rotate = helices_to_include;
                    inital_conditions_rotaion_coefficient = 0;
                    inital_conditions_translation_coefficient = 0;

                    bounds_rotation_coefficient = rot_index_vector(rot_index);
                    bounds_translation_coefficient = trans_index_vector(trans_index);

                    sim = Simulation();
                    sim.Initialise_Optim_RotTrans_Sim(domain_index,helix_indeces,helices_to_include,helices_to_rotate);

                    result = sim.Run_Optim_RotTrans_Sim(...
                        inital_conditions_rotaion_coefficient,...
                        inital_conditions_translation_coefficient,...
                        bounds_rotation_coefficient,...
                        bounds_translation_coefficient,...
                        100);
                    domain_results = [domain_results,result];
                end
            end
%                 filename_str = sprintf('optim_paper/domain_results_r%d', domain_index);
%                 save(filename_str,'domain_results');
%                 fprintf("Domain R%d results saved.\n",domain_index);
%                 results = [results,domain_results];
            results = domain_results;
        end
    end
end