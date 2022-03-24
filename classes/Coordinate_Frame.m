classdef Coordinate_Frame < handle
    %COORDINATE_FRAME Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        base
        children
        parent
        object
        
        position_start
        position_end
        
        rotation_from_parent
        translation_from_parent
        
        locked_dof
    end
    
    methods
        function self = Coordinate_Frame(object,parent,rotation,translation)
            
            % Init Properties
            self.base = false;
            self.children = [];
            self.parent = [];
            self.object= [];
            
            self.position_start = zeros(3,1);
            self.position_end = zeros(3,1);
            
            self.rotation_from_parent = eye(3);
            self.translation_from_parent = zeros(3,1);
            
            self.locked_dof = [0;0;0];
                        
            
            % Arguement to Property assignment
            switch nargin
                case 0 % No Argmuents. Assumes the creation of Base Frame
                    self.base = true;
                case 1
                    self.object = object;
                case 2
                    self.object = object;
                    self.parent = parent;
                case 3
                    self.object = object;
                    self.parent = parent;
                    self.rotation_from_parent = rotation;
                case 4
                    self.object = object;
                    self.parent = parent;
                    self.rotation_from_parent = rotation;
                    self.translation_from_parent = Vector_Helper.Make_Column(translation);
                otherwise
                    fprintf("ERROR: unexpected number of input arguements.\n");
            end
            
            % Object Handling
            switch (class(self.object))
                case 'double'
                    %fprintf("Initialising Base Frame. \n");
                case 'Helix_Object'
                    self.object.parent_frame = self;
                    self.Process_Helix();
                    self.position_end = [norm(self.object.position_Two - self.object.position_One);0;0];
                case 'Joint_Object'
                    self.object.parent_frame = self;
                    self.position_end = self.object.position_Two';
                otherwise
                    fprintf("ERROR: Unexpected object type.\n");
            end
        end
        function addChild(self,object,rotation,translation)
            switch nargin
                case 2
                    self.children = [self.children,Coordinate_Frame(object,self)];
                case 3
                    self.children = [self.children,Coordinate_Frame(object,self,rotation)];
                case 4
                    self.children = [self.children,Coordinate_Frame(object,self,rotation,Vector_Helper.Make_Column(translation))];
            end
        end
        function insertJointBefore(self)
            % BIG WARNING: ASSUME LINEAR CHAIN!!!
            joint = Coordinate_Frame(Joint_Object());
            joint.rotation_from_parent = self.rotation_from_parent;
            joint.translation_from_parent = self.translation_from_parent;
            self.translation_from_parent = zeros(3,1);
            self.rotation_from_parent = eye(3);
            
            joint.parent = self.parent;
            joint.children = self;
            index = 1;
            for i = 1:length(self.parent.children)
                if (isequal(self.parent.children(i),self))
                    index = i;
                end
            end
            self.parent.children(index) = joint;
            self.parent = joint;
        end
        function ProcessHelixRotation_OLD(self)
            % Update original positions and orientiations before
            % translations
            % get relative position of first sidechain that is not Gly
            first_sidechain_pos = [0,0,0];
            for i = 1:length(self.object.sidechain_positions)
                first_sidechain_pos = self.object.sidechain_positions(i,:);
                if ~all(first_sidechain_pos == 0)
                    break
                end
            end
            % get direction_vector from helix_axis to first sidechain. This
            % is considered the UPWARDS direction for the helix, relative
            % to the pdb axes
            intersection_point = Point_Line_Intersect(self.object.position_One, self.object.position_Two, first_sidechain_pos);
            direction_vector =  first_sidechain_pos - intersection_point
            direction_unit_vector = direction_vector / norm(direction_vector);
            self.object.original_up = direction_unit_vector;
            self.object.original_pos = self.object.position_Two;
            
            % Workaround for Gly not being removed.
            indexs = find(self.object.sidechain_positions(:,1) == 0);  
            
            % Get translation vector from p1 to origin and translate p1, p2
            % and sidechains
            toOriginTranslate = self.object.position_One;
            self.object.position_One = self.object.position_One - toOriginTranslate;
            self.object.position_Two = self.object.position_Two - toOriginTranslate;
            self.object.sidechain_positions = self.object.sidechain_positions - toOriginTranslate;
            
                      
            % Workaround for Gly not being removed.
            for i = 1:length(indexs)
                self.object.sidechain_positions(indexs(i),:) = [0,0,0];
            end
            
            % Get Rotation matrix from current helix vector to [1;0;0] and
            % apply to p1, p2 and sidechains
            R = self.Get_Object_Conversion_Rotation_Matrix();
            self.object.position_One = (R * self.object.position_One')'; 
            self.object.position_Two = (R * self.object.position_Two')';
            for i = 1:length(self.object.sidechain_positions(:,1))
                self.object.sidechain_positions(i,:) = (R * self.object.sidechain_positions(i,:)')';
            end
            
            % Get new first the y-axis sidechain direction and rotate helix around
            % x-axis such that it is aligned with 
            % Check whether first sidechain or last sidechain is closest to
            % local origin
            first_sidechain_pos = self.object.sidechain_positions(1,:);
            last_sidechain_pos = self.object.sidechain_positions(length(self.object.sidechain_positions(:,1)),:);
            dist_to_first = norm(first_sidechain_pos - 0);
            dist_to_last = norm(last_sidechain_pos - 0);
            % If Euclidian distance to last sidechain is shorter than to
            % first, helix to be flipped along X-axis length.
            if (dist_to_last < dist_to_first)
                % Translate pos two to origin
                self.object.position_One = self.object.position_One - self.object.position_Two;
                self.object.sidechain_positions = self.object.sidechain_positions - self.object.position_Two;
                self.object.position_Two = self.object.position_Two - self.object.position_Two;
                % Rotate everything around the Y axis by the difference in
                % angle between Z axis and sidechain
                % Workaround for Gly not being removed.
                for i = 1:length(indexs)
                    self.object.sidechain_positions(indexs(i),:) = [0,0,0];
                end
                R = self.rotationMatrix([0;0;1],180);
                self.object.position_One = (R * self.object.position_One')';
                self.object.position_Two = (R * self.object.position_Two')';
                for i = 1:length(self.object.sidechain_positions(:,1))
                    self.object.sidechain_positions(i,:) = (R * self.object.sidechain_positions(i,:)')';
                end
            end
            % 
            % get direction of sidechain
                first_sidechain_pos = self.object.sidechain_positions(1,:);
                sidechain_direction = [0,first_sidechain_pos(2),first_sidechain_pos(3)];
                sidechain_direction = sidechain_direction / norm(sidechain_direction);
                angleFromZ = acosd(dot([0,0,1],sidechain_direction));
                angleFromY = acosd(dot([0,1,0],sidechain_direction));
                if (angleFromY <= 90.0)
                    R = self.rotationMatrix([1;0;0],angleFromZ);
                else
                    R = self.rotationMatrix([1;0;0],-angleFromZ);
                end
                % Workaround for Gly not being removed.
                indexs = find(self.object.sidechain_positions(:,1) == 0);
                % Workaround for Gly not being removed.
                for i = 1:length(indexs)
                    self.object.sidechain_positions(indexs(i),:) = [0,0,0];
                end
                self.object.position_One = (R * self.object.position_One')';
                self.object.position_Two = (R * self.object.position_Two')';
                for i = 1:length(self.object.sidechain_positions(:,1))
                    self.object.sidechain_positions(i,:) = (R * self.object.sidechain_positions(i,:)')';
                end
        end
        % New function created on 22 sep 2019
        function Process_Helix(self)
            % Tansform helix to allign to local coordinate frame
            
            % Pre-process sidechain positions because of Glycine
            indexs = find(self.object.sidechain_positions(:,1) == 0);
            % Translate positions and sidechains by position one
            translation_vec = self.object.position_One;
            self.object.position_One = self.object.position_One - translation_vec;
            self.object.position_Two = self.object.position_Two - translation_vec;
            self.object.sidechain_positions = self.object.sidechain_positions - translation_vec;
            % Post-process sidechain positions because of Glycine
            for i = 1:length(indexs)
                self.object.sidechain_positions(indexs(i),:) = [0,0,0];
            end
            
            % Apply rotation to helix to align central axis with X axis
            R = self.Get_Object_Conversion_Rotation_Matrix();
            self.object.position_One = (R * self.object.position_One')';
            self.object.position_Two = (R * self.object.position_Two')';
            for i = 1:length(self.object.sidechain_positions)
                self.object.sidechain_positions(i,:) = (R * self.object.sidechain_positions(i,:)');
            end
            
            % Apply rotation to helix to align first non-glycine sidechain
            % with the Z axis 
            first_sidechain_pos = [0,0,0];
            for i = 1:length(self.object.sidechain_positions)
                first_sidechain_pos = self.object.sidechain_positions(i,:);
                if (~all(first_sidechain_pos == 0))
                    break
                end
            end
            sidechain_dir_vec = [0,first_sidechain_pos(2),first_sidechain_pos(3)];
            sidechain_dir_vec = sidechain_dir_vec / norm(sidechain_dir_vec);
            angle_from_y = acosd(dot([0,1,0],sidechain_dir_vec));
            angle_from_z = acosd(dot([0,0,1],sidechain_dir_vec));
            if (angle_from_y <= 90.0)
                R = self.rotationMatrix([1;0;0],angle_from_z);
            else
                R = self.rotationMatrix([1;0;0],-angle_from_z);
            end
            % Pre-process sidechain positions because of Glycine
            indexs = find(self.object.sidechain_positions(:,1) == 0);
            % Apply Rotation
            self.object.position_One = (R * self.object.position_One')';
            self.object.position_Two = (R * self.object.position_Two')';
            for i = 1:length(self.object.sidechain_positions)
                self.object.sidechain_positions(i,:) = (R * self.object.sidechain_positions(i,:)')';
            end
            % Post-process sidechain positions because of Glycine
            for i = 1:length(indexs)
                self.object.sidechain_positions(indexs(i),:) = [0,0,0];
            end
        end
        
        function output = Get_Object_Conversion_Rotation_Matrix(self)
            vec_one = self.object.position_Two - self.object.position_One;
            vec_one = vec_one / norm(vec_one);
            vec_two = [1,0,0];
            vec_cross = cross(vec_one,vec_two);
            vec_dot = dot(vec_one,vec_two);
            %             vec_skew = [0,-vec_cross(3),vec_cross(2);vec_cross(3),0,-vec_cross(1);-vec_cross(2),vec_cross(1),0];
            %             output = eye(3) + vec_skew + (vec_skew^2 * ((1 - vec_dot)/(norm(vec_cross)^2)));
            output = self.rotationMatrix(vec_cross,acosd(vec_dot));
            
        end
        
        function output = Get_Figure_Data_World_Frame(self)
            figure_data = struct("type",[],"data",[]);
            switch (class(self.object))
                case 'Helix_Object'
                    figure_data.type = "Helix_Object";
                    % Get this object helix positions for figure generation in world frame
                    helixData = self.object.Get_Figure_Data();
                    for i = 1:length(helixData.positions(:,1))
                        % Deal with first helix position
                        helixPos = helixData.positions(1,:)';
                        helixPos = self.Get_Position_World_Frame(helixPos);
                        helixData.positions(1,:) = helixPos';
                        
                        % Deal with second helix position
                        helixPos = helixData.positions(2,:)';
                        helixPos = self.Get_Position_World_Frame(helixPos);
                        helixData.positions(2,:) = helixPos';
                    end
                    % Get world frame data for cylinder data
                    for i = 1:length(helixData.cylinder.X(1,:))
                        % Deal with first cylinder head
                        cylinderPos = [helixData.cylinder.X(1,i);helixData.cylinder.Y(1,i);helixData.cylinder.Z(1,i)];
                        cylinderPos = self.Get_Position_World_Frame(cylinderPos);
                        helixData.cylinder.X(1,i) = cylinderPos(1);
                        helixData.cylinder.Y(1,i) = cylinderPos(2);
                        helixData.cylinder.Z(1,i) = cylinderPos(3);
                        
                        % Deal with second cylinder head
                        
                        cylinderPos = [helixData.cylinder.X(2,i);helixData.cylinder.Y(2,i);helixData.cylinder.Z(2,i)];
                        cylinderPos = self.Get_Position_World_Frame(cylinderPos);
                        helixData.cylinder.X(2,i) = cylinderPos(1);
                        helixData.cylinder.Y(2,i) = cylinderPos(2);
                        helixData.cylinder.Z(2,i) = cylinderPos(3);
                    end
                    for i = 1:length(helixData.sidechains(:,1))
                        % Deal with sidechain positions
                        sidechainPos = helixData.sidechains(i,:)';
                        sidechainPos = self.Get_Position_World_Frame(sidechainPos);
                        helixData.sidechains(i,:) = sidechainPos';
                    end
                    for i = 1:length(helixData.sidechain_mount.X)
                        mountPos = [helixData.sidechain_mount.X(i);helixData.sidechain_mount.Y(i);helixData.sidechain_mount.Z(i)];
                        mountPos = self.Get_Position_World_Frame(mountPos);
                        helixData.sidechain_mount.X(i) = mountPos(1);
                        helixData.sidechain_mount.Y(i) = mountPos(2);
                        helixData.sidechain_mount.Z(i) = mountPos(3);
                    end
                    figure_data.data = helixData;
                case 'Joint_Object'
                    figure_data.type = "Joint_Object";
                    jointData.positions = [self.object.position_One;self.object.position_Two];
                    figure_data.data = jointData;
                otherwise
                    fprintf("ERROR: Unexpected object type.\n");
            end
            output = figure_data;
        end
        function h = Generate_Figure_From_Data(~,figure_data)
            switch (figure_data.type)
                case "Helix_Object"
                    X = figure_data.data.cylinder.X;
                    Y = figure_data.data.cylinder.Y;
                    Z = figure_data.data.cylinder.Z;
                    sidechains = figure_data.data.sidechains;
                    sidechain_mounts(:,1) = figure_data.data.sidechain_mount.X;
                    sidechain_mounts(:,2) = figure_data.data.sidechain_mount.Y;
                    sidechain_mounts(:,3) = figure_data.data.sidechain_mount.Z;
                    sX = [sidechain_mounts(:,1)';sidechains(:,1)'];
                    sY = [sidechain_mounts(:,2)';sidechains(:,2)'];
                    sZ = [sidechain_mounts(:,3)';sidechains(:,3)'];
                    h = figure();
                    axis equal
                    hold on
                    surf(X,Y,Z,'FaceColor','r','EdgeAlpha',0,'FaceAlpha',1.0); % Cylinder Surface
                    patch(X(1,:),Y(1,:),Z(1,:),[1 0 0],'FaceAlpha',1.0); % Cylinder End Face
                    patch(X(2,:),Y(2,:),Z(2,:),[1 0 0],'FaceAlpha',1.0); % Cylinder End Face
                    plot3(sidechains(:,1),sidechains(:,2),sidechains(:,3),'.','Color','blue');
                    plot3(sX,sY,sZ,'Color','green');
                    hold off
                case "Joint_Object"
                    h = figure();
                    axis equal
                    hold on
                    plot3(figure_data.data.positions(:,1),figure_data.data.positions(:,2),figure_data.data.positions(:,3),'Color','blue');
                    hold off
                otherwise
                    fprintf("ERROR: Unexpected figure_data type.\n");
            end
        end
        function h = Generate_Figure_From_Multiple_Data(~,figure_datas)
            h = figure();
            residue_properties = load("residueProperties.mat");
            residue_properties = residue_properties.residueProperties;
            hold on
            for i = 1:length(figure_datas)
                figure_data = figure_datas(i);
                switch (figure_data.type)
                    case "Helix_Object"
                        X = figure_data.data.cylinder.X;
                        Y = figure_data.data.cylinder.Y;
                        Z = figure_data.data.cylinder.Z;
                        sidechains = figure_data.data.sidechains;
                        sidechain_sizes = figure_data.data.sidechain_sizes;
                        sidechain_mounts = [];
                        sidechain_mounts(:,1) = figure_data.data.sidechain_mount.X;
                        sidechain_mounts(:,2) = figure_data.data.sidechain_mount.Y;
                        sidechain_mounts(:,3) = figure_data.data.sidechain_mount.Z;
                        sX = [sidechain_mounts(:,1)';sidechains(:,1)'];
                        sY = [sidechain_mounts(:,2)';sidechains(:,2)'];
                        sZ = [sidechain_mounts(:,3)';sidechains(:,3)'];
                        surf(X,Y,Z,'FaceColor','r','EdgeAlpha',0,'FaceAlpha',1.0); % Cylinder Surface
                        patch(X(1,:),Y(1,:),Z(1,:),[1 0 0],'FaceAlpha',1.0); % Cylinder End Face
                        patch(X(2,:),Y(2,:),Z(2,:),[1 0 0],'FaceAlpha',1.0); % Cylinder End Face
                        plot3(sidechains(:,1),sidechains(:,2),sidechains(:,3),'.','Color','blue');
                        for j = 1:length(sidechains(:,1)) % Sidechain Spheres
                            PlotSphere(sidechains(j,1),sidechains(j,2),sidechains(j,3),sidechain_sizes(j));
                        end
                        plot3(sX,sY,sZ,'Color','green');
                        plot3(sX(:,1),sY(:,1),sZ(:,1),'Color','m');
                        if (1)
                            xlabel('X');
                            ylabel('Y');
                            zlabel('Z');
                        end
                    case "Joint_Object"
                        plot3(figure_data.data.positions(:,1),figure_data.data.positions(:,2),figure_data.data.positions(:,3),'Color','#7E2F8E','LineWidth',2);
                    otherwise
                        fprintf("ERROR: Unexpected figure_data type.\n");
                end
            end
            hold off
            axis equal
        end
        function h = Generate_Figure_From_Multiple_Data_With_Attributes(~,figure_datas,attributes)
            h = figure();
            residue_properties = load("residueProperties.mat");
            residue_properties = residue_properties.residueProperties;
            
            % Check Attributes arguement is valid
            % Do Here...
                       
            hold on
            for i = 1:length(figure_datas)
                figure_data = figure_datas(i);
                switch (figure_data.type)
                    case "Helix_Object"
                        X = figure_data.data.cylinder.X;
                        Y = figure_data.data.cylinder.Y;
                        Z = figure_data.data.cylinder.Z;
                        sidechains = figure_data.data.sidechains;
                        sidechain_sizes = figure_data.data.sidechain_sizes;
                        sidechain_mounts = [];
                        sidechain_mounts(:,1) = figure_data.data.sidechain_mount.X;
                        sidechain_mounts(:,2) = figure_data.data.sidechain_mount.Y;
                        sidechain_mounts(:,3) = figure_data.data.sidechain_mount.Z;
                        sX = [sidechain_mounts(:,1)';sidechains(:,1)'];
                        sY = [sidechain_mounts(:,2)';sidechains(:,2)'];
                        sZ = [sidechain_mounts(:,3)';sidechains(:,3)'];
                        helix_colour = [1,0,0];
                        if (figure_data.data.helix_index <= length(attributes.helix_colours(:,1)))
                            if (~all(attributes.helix_colours(figure_data.data.helix_index,:) == -1))
                                helix_colour = attributes.helix_colours(figure_data.data.helix_index,:);
                            end
                        end
                        % Helix Cylinder
                        surf(X,Y,Z,'FaceColor',helix_colour,'EdgeAlpha',0,'FaceAlpha',1.0); % Cylinder Surface
                        patch(X(1,:),Y(1,:),Z(1,:),helix_colour,'FaceAlpha',1.0); % Cylinder End Face
                        patch(X(2,:),Y(2,:),Z(2,:),helix_colour,'FaceAlpha',1.0); % Cylinder End Face
                        
                        sidechainLineWidth = 0.5;
                        sidechainMarkerSize = 6;
                        sidechainStickLineStyle = '-';
                        primarySidechainLineWidth = 0.5;
                        primarySidechainMarkerSize = 6;
                        if (ismember("primary_sidechain_emphesis",attributes.mode))
                            sidechainLineWidth = 0.5;
                            sidechainStickLineStyle = ':';
                            sidechainMarkerSize = 0.5;
                            primarySidechainLineWidth = 2;
                            primarySidechainMarkerSize = 12;
                        end
                        if (ismember("simple_primary_sidechain_emphesis",attributes.mode))
                            primarySidechainLineWidth = 3;
%                             primarySidechainMarkerSize = 12;
                        end
                        primarySidechainStickColour = 'm';
                        if (~isempty(attributes.primary_sidechain_stick_colour))
                            primarySidechainStickColour = attributes.primary_sidechain_stick_colour;
                        end
                        if (~ismember("hide_normal_sidechains",attributes.mode))
                            % Normal Sidechain Position
                            plot3(sidechains(2:end,1),sidechains(2:end,2),sidechains(2:end,3),...
                                'Marker','.',...
                                'LineStyle','none',...
                                'MarkerSize',sidechainMarkerSize,...
                                'Color','blue');
                            if (~ismember("dots",attributes.mode))
                                for j = 1:length(sidechains(:,1)) % Sidechain Spheres
                                    if 0
                                        PlotSphere(sidechains(j,1),sidechains(j,2),sidechains(j,3),sidechain_sizes(j),'g');
                                    else
                                        PlotSphere(sidechains(j,1),sidechains(j,2),sidechains(j,3),sidechain_sizes(j));
                                    end
%                                     PlotSphere(sidechains(j,1),sidechains(j,2),sidechains(j,3),sidechain_sizes(j));
                                end
                            end
                            % Normal Sidechain Sticks
                            plot3(sX(:,2:end),sY(:,2:end),sZ(:,2:end),...
                                'LineStyle',sidechainStickLineStyle,...
                                'Color',[0,1,0],...
                                'LineWidth',sidechainLineWidth);
                        end
                        % Primary Sidechain Stick
                        plot3(sX(:,1),sY(:,1),sZ(:,1),...
                            'LineStyle','-',...
                            'Color',primarySidechainStickColour,...
                            'LineWidth',primarySidechainLineWidth);
                        % Primary Sidechain Position
                        plot3(sidechains(1,1),sidechains(1,2),sidechains(1,3),...
                            'Marker','.',...
                            'LineStyle','none',...
                            'Color','blue',...
                            'MarkerSize',primarySidechainMarkerSize);
                        if (1)
                            xlabel('X');
                            ylabel('Y');
                            zlabel('Z');
                        end
                    case "Joint_Object"
                        plot3(figure_data.data.positions(:,1),figure_data.data.positions(:,2),figure_data.data.positions(:,3),'Color','#7E2F8E','LineWidth',2);
                    otherwise
                        fprintf("ERROR: Unexpected figure_data type.\n");
                end
            end
            hold off
            axis equal
        end
        function h = Generate_Figure_From_Multiple_Data_Dots(~,figure_datas)
            h = figure();
            residue_properties = load("residueProperties.mat");
            residue_properties = residue_properties.residueProperties;
            hold on
            for i = 1:length(figure_datas)
                figure_data = figure_datas(i);
                switch (figure_data.type)
                    case "Helix_Object"
                        X = figure_data.data.cylinder.X;
                        Y = figure_data.data.cylinder.Y;
                        Z = figure_data.data.cylinder.Z;
                        sidechains = figure_data.data.sidechains;
                        sidechain_sizes = figure_data.data.sidechain_sizes;
                        sidechain_mounts = [];
                        sidechain_mounts(:,1) = figure_data.data.sidechain_mount.X;
                        sidechain_mounts(:,2) = figure_data.data.sidechain_mount.Y;
                        sidechain_mounts(:,3) = figure_data.data.sidechain_mount.Z;
                        sX = [sidechain_mounts(:,1)';sidechains(:,1)'];
                        sY = [sidechain_mounts(:,2)';sidechains(:,2)'];
                        sZ = [sidechain_mounts(:,3)';sidechains(:,3)'];
                        surf(X,Y,Z,'FaceColor','r','EdgeAlpha',0,'FaceAlpha',1.0); % Cylinder Surface
                        patch(X(1,:),Y(1,:),Z(1,:),[1 0 0],'FaceAlpha',1.0); % Cylinder End Face
                        patch(X(2,:),Y(2,:),Z(2,:),[1 0 0],'FaceAlpha',1.0); % Cylinder End Face
                        plot3(sidechains(:,1),sidechains(:,2),sidechains(:,3),'.','Color','blue');
                        plot3(sX,sY,sZ,'Color','green');
                        plot3(sX(:,1),sY(:,1),sZ(:,1),'Color','m');
                        if (1)
                            xlabel('X');
                            ylabel('Y');
                            zlabel('Z');
                        end
                    case "Joint_Object"
                        plot3(figure_data.data.positions(:,1),figure_data.data.positions(:,2),figure_data.data.positions(:,3),'Color','#7E2F8E','LineWidth',2);
                    otherwise
                        fprintf("ERROR: Unexpected figure_data type.\n");
                end
            end
            hold off
            axis equal
        end
        function h = Generate_Figure_From_Multiple_Data_Thick_Lines(~,figure_datas)
            h = figure();
            axis equal
            hold on
            for i = 1:length(figure_datas)
                figure_data = figure_datas(i);
                switch (figure_data.type)
                    case "Helix_Object"
                        X = figure_data.data.cylinder.X;
                        Y = figure_data.data.cylinder.Y;
                        Z = figure_data.data.cylinder.Z;
                        sidechains = figure_data.data.sidechains;
                        sidechain_mounts = [];
                        sidechain_mounts(:,1) = figure_data.data.sidechain_mount.X;
                        sidechain_mounts(:,2) = figure_data.data.sidechain_mount.Y;
                        sidechain_mounts(:,3) = figure_data.data.sidechain_mount.Z;
                        sX = [sidechain_mounts(:,1)';sidechains(:,1)'];
                        sY = [sidechain_mounts(:,2)';sidechains(:,2)'];
                        sZ = [sidechain_mounts(:,3)';sidechains(:,3)'];
                        surf(X,Y,Z,'FaceColor','r','EdgeAlpha',0,'FaceAlpha',1.0); % Cylinder Surface
                        patch(X(1,:),Y(1,:),Z(1,:),[1 0 0],'FaceAlpha',1.0); % Cylinder End Face
                        patch(X(2,:),Y(2,:),Z(2,:),[1 0 0],'FaceAlpha',1.0); % Cylinder End Face
                        plot3(sidechains(:,1),sidechains(:,2),sidechains(:,3),'.','Color','blue');
                        plot3(sX,sY,sZ,'Color','green','LineWidth',1);
                        plot3(sX(:,1),sY(:,1),sZ(:,1),'Color','m','LineWidth',2);
                        if (1)
                            xlabel('X');
                            ylabel('Y');
                            zlabel('Z');
                        end
                    case "Joint_Object"
                        plot3(figure_data.data.positions(:,1),figure_data.data.positions(:,2),figure_data.data.positions(:,3),'Color','blue');
                    otherwise
                        fprintf("ERROR: Unexpected figure_data type.\n");
                end
            end
            hold off
        end
        function h = Generate_Figure_From_Multiple_Data_Hydrophobic(~,figure_datas)
            h = figure();
            gen_residue_properties;
            hydro_scale = sort(residue_properties.hydrophobicity_kd);
            hold on
            for i = 1:length(figure_datas)
                figure_data = figure_datas(i);
                switch (figure_data.type)
                    case "Helix_Object"
                        X = figure_data.data.cylinder.X;
                        Y = figure_data.data.cylinder.Y;
                        Z = figure_data.data.cylinder.Z;
                        sidechains = figure_data.data.sidechains;
                        sidechain_sizes = figure_data.data.sidechain_sizes;
                        sidechain_mounts = [];
                        sidechain_mounts(:,1) = figure_data.data.sidechain_mount.X;
                        sidechain_mounts(:,2) = figure_data.data.sidechain_mount.Y;
                        sidechain_mounts(:,3) = figure_data.data.sidechain_mount.Z;
                        sX = [sidechain_mounts(:,1)';sidechains(:,1)'];
                        sY = [sidechain_mounts(:,2)';sidechains(:,2)'];
                        sZ = [sidechain_mounts(:,3)';sidechains(:,3)'];
                        surf(X,Y,Z,'FaceColor','r','EdgeAlpha',0,'FaceAlpha',1.0); % Cylinder Surface
                        patch(X(1,:),Y(1,:),Z(1,:),[1 0 0],'FaceAlpha',1.0); % Cylinder End Face
                        patch(X(2,:),Y(2,:),Z(2,:),[1 0 0],'FaceAlpha',1.0); % Cylinder End Face
                        plot3(sidechains(:,1),sidechains(:,2),sidechains(:,3),'.','Color','blue');
                        
                        hydro_scale(hydro_scale==0) = [];
                        hydro_scale = unique(hydro_scale);
                        
                        hydroMax = max(hydro_scale);
                        hydroMin = min(hydro_scale);
                        sidechainHydropathy = figure_data.data.sidechain_hydropathy;
                        colourRange = linspace(180, 360, length(hydro_scale)); 
                        colours = zeros(length(sidechainHydropathy), 3);
                        
                        % as this is currently setup, red is hydrophobic
                        % and blue is hydrophilic
                        
                        for j = 1:length(sidechains(:,1)) % Sidechain Spheres
                            colour = hsl2rgb(colourRange(hydro_scale==sidechainHydropathy(j)),1,0.5);
                            colour = colour ./ 255;
                            colours(j,:) = colour;
                        end
                        for j = 1:length(sidechains(:,1)) % Sidechain Spheres
                            colour = hsl2rgb(colourRange(hydro_scale==sidechainHydropathy(6)),1,0.5);
                            colour = colour ./ 255;
                            PlotSphere(sidechains(j,1),sidechains(j,2),sidechains(j,3),sidechain_sizes(j),colours(j,:));
                        end
                        plot3(sX,sY,sZ,'Color','green');
                        plot3(sX(:,1),sY(:,1),sZ(:,1),'Color','m');
                        if (1)
                            xlabel('X');
                            ylabel('Y');
                            zlabel('Z');
                        end
                    case "Joint_Object"
                        plot3(figure_data.data.positions(:,1),figure_data.data.positions(:,2),figure_data.data.positions(:,3),'Color','#7E2F8E','LineWidth',2);
                    otherwise
                        fprintf("ERROR: Unexpected figure_data type.\n");
                end
            end
            hold off
            axis equal
        end
        function h = Generate_Complete_Figure(self)
            assert(self.base,"ERROR: Generate_Complete_Figure: is not base...\n")
            h = self.Generate_Figure_From_Multiple_Data(self.CS_Generate_Complete_Figure());
        end
        function h = Generate_Complete_Figure_With_Attributes(self,attributes)
            assert(self.base,"ERROR: Generate_Complete_Figure_With_Attributes: is not base...\n");
            h = self.Generate_Figure_From_Multiple_Data_With_Attributes(self.CS_Generate_Complete_Figure(),attributes);
        end
        function h = Generate_Complete_Figure_Dots(self)
            assert(self.base,"ERROR: Generate_Complete_Figure: is not base...\n");
            h = self.Generate_Figure_From_Multiple_Data_Dots(self.CS_Generate_Complete_Figure());
        end
        function h = Generate_Complete_Figure_Thick_Lines(self)
            assert(self.base,"ERROR: Generate_Complete_Figure: is not base...\n")
            h = self.Generate_Figure_From_Multiple_Data_Thick_Lines(self.CS_Generate_Complete_Figure());
        end
        function h = Generate_Complete_Figure_Hydrophobic(self)
            assert(self.base,"ERROR: Generate_Complete_Figure: is not base...\n")
            h = self.Generate_Figure_From_Multiple_Data_Hydrophobic(self.CS_Generate_Complete_Figure());
        end
        function h = Generate_Coordinate_Frame_Graph(self)
            h = figure();
            xlabel("x");
            ylabel("y");
            zlabel("z");
            axis equal
            hold on
            datas = self.Get_Coordinate_Frame_Data();
            for i = 1:length(datas)
                data = datas(i);
                Draw_Cube(data.pos,0.5);
                for j = 1:size(data.pos_children,2)
                    plot3([data.pos(1),data.pos_children(1,j)],[data.pos(2),data.pos_children(2,j)],[data.pos(3),data.pos_children(3,j)]);
                end
            end
            
            hold off;
            
            function Draw_Cube(centre_pos,radius)
                x = [
                    centre_pos(1)-radius,centre_pos(1)+radius,centre_pos(1)+radius,centre_pos(1)-radius;
                    centre_pos(1)-radius,centre_pos(1)+radius,centre_pos(1)+radius,centre_pos(1)-radius;
                    centre_pos(1)-radius,centre_pos(1)+radius,centre_pos(1)+radius,centre_pos(1)-radius;
                    centre_pos(1)-radius,centre_pos(1)-radius,centre_pos(1)-radius,centre_pos(1)-radius;
                    centre_pos(1)-radius,centre_pos(1)+radius,centre_pos(1)+radius,centre_pos(1)-radius;
                    centre_pos(1)+radius,centre_pos(1)+radius,centre_pos(1)+radius,centre_pos(1)+radius;];
                y = [
                    centre_pos(2)-radius,centre_pos(2)-radius,centre_pos(2)+radius,centre_pos(2)+radius;
                    centre_pos(2)-radius,centre_pos(2)-radius,centre_pos(2)+radius,centre_pos(2)+radius;
                    centre_pos(2)-radius,centre_pos(2)-radius,centre_pos(2)-radius,centre_pos(2)-radius;
                    centre_pos(2)-radius,centre_pos(2)-radius,centre_pos(2)+radius,centre_pos(2)+radius;
                    centre_pos(2)+radius,centre_pos(2)+radius,centre_pos(2)+radius,centre_pos(2)+radius;
                    centre_pos(2)-radius,centre_pos(2)-radius,centre_pos(2)+radius,centre_pos(2)+radius;];
                z = [
                    centre_pos(3)-radius,centre_pos(3)-radius,centre_pos(3)-radius,centre_pos(3)-radius;
                    centre_pos(3)+radius,centre_pos(3)+radius,centre_pos(3)+radius,centre_pos(3)+radius;
                    centre_pos(3)-radius,centre_pos(3)-radius,centre_pos(3)+radius,centre_pos(3)+radius;
                    centre_pos(3)-radius,centre_pos(3)+radius,centre_pos(3)+radius,centre_pos(3)-radius;
                    centre_pos(3)-radius,centre_pos(3)-radius,centre_pos(3)+radius,centre_pos(3)+radius;
                    centre_pos(3)-radius,centre_pos(3)+radius,centre_pos(3)+radius,centre_pos(3)-radius;];
                patch(x',y',z','y');
            end
            
            
        end
        function output = Get_Coordinate_Frame_Data(self)
                output = [];
                if (self.base)
                    data.pos = self.Get_Position_World_Frame(zeros(3,1));
                else
                    data.pos = self.parent.Get_Position_World_Frame(self.translation_from_parent);
                end
                    if (~isempty(self.children))
                        for m = 1:length(self.children)
                            data.pos_children(:,m) = self.Get_Position_World_Frame(self.children(m).translation_from_parent);
                        end
                    else
                        data.pos_children = [];
                    end    
                    output = [output,data];
                for nChildren = 1:length(self.children)
                    output = l_Recursive_Get_Next_Coordinate_Frame_Data(self.children(nChildren),output);
                end
                
                function output = l_Recursive_Get_Next_Coordinate_Frame_Data(a_frame,a_output)
                    coordinate_frame_data.pos = a_frame.parent.Get_Position_World_Frame(a_frame.translation_from_parent);
                    if (~isempty(a_frame.children))
                        for n = 1:length(a_frame.children)
                            coordinate_frame_data.pos_children(:,n) = a_frame.Get_Position_World_Frame(a_frame.children(n).translation_from_parent);
                        end
                    else
                        coordinate_frame_data.pos_children = [];
                    end
                    output = [a_output,coordinate_frame_data];
                    if (~isempty(a_frame.children))
                        for i = 1:length(a_frame.children)
                            output = l_Recursive_Get_Next_Coordinate_Frame_Data(a_frame.children(i),output);
                        end
                    end
                end
            end
        
        function Set_Position(self,positionVector)
            self.translation_from_parent = Vector_Helper.Make_Column(positionVector);
        end
        % TODO: Implement
        function Set_Position_Axial_By_Offset(self,offset)
            %SETS THE TRANSLATION OF THE FRAME ALONG ITS AXIS BASED ON AN OFFSET FROM ORIGINAL POSITION.
            self.translation_from_parent = Vector_Helper.Make_Column([offset,0,0]);
        end
        function Translate_Frame(self,translationVector)
            self.translation_from_parent = self.translation_from_parent + Vector_Helper.Make_Column(translationVector);
        end
        
        function Set_Rotation_RM(self,rotationMatrix)
            self.rotation_from_parent = rotationMatrix;
        end
        function Set_Rotation_AA(self,axis,angle)
            R = Coordinate_Frame.rotationMatrix(axis,angle);
            self.Set_Rotation_RM(R);
        end
        function Set_Orientation(self,axis)
            R = self.Generate_Rotation_from_VecOne_to_VecTwo([1,0,0] * self.rotation_from_parent,axis);
            self.Set_Rotation_RM(R);
        end
        function Rotate_Frame_RM(self,rotationMatrix)
            self.rotation_from_parent = rotationMatrix * self.rotation_from_parent;
        end
        function Rotate_Frame_AA(self,axis,angle)
            R = Coordinate_Frame.rotationMatrix(axis,angle);
            self.Rotate_Frame_RM(R);
        end
        
        function Set_Frame_to_Helix_Original(self)
            self.Set_Frame_to_Helix_Original_Position();
            self.Set_Frame_to_Helix_Original_Orientation();
            self.Set_Frame_to_Helix_Original_Rotation();
        end
        function Set_Frame_to_Helix_Original_Position(self)
            self.Set_Position([0;0;0]);
            self.Set_Position(self.object.original_pos');
        end
        function Set_Frame_to_Helix_Original_Orientation(self)
            switch (class(self.object))
                case ('Helix_Object')
                    self.Set_Orientation(self.object.original_direction);
                otherwise
                    fprintf("Object of class %s, is not a 'Helix_Object'",class(self.object));
            end
        end
        function Set_Frame_to_Helix_Original_Rotation(self)
            switch (class(self.object))
                case ('Helix_Object')
                    first_sc = [0,0,0];
                    for i = 1:length(self.object.sidechain_positions)
                        first_sc = self.object.sidechain_positions(i,:);
                        if (~all(first_sc == 0))
                            break
                        end
                    end
%                     cur_dir = first_sc - Point_Line_Intersect(self.object.position_One, self.object.position_Two,first_sc);
                    cur_dir = self.Get_Position_World_Frame(first_sc') - Point_Line_Intersect(self.Get_Position_World_Frame(self.object.position_One'),self.Get_Position_World_Frame(self.object.position_Two'),self.Get_Position_World_Frame(first_sc'));
                    cur_dir = cur_dir/norm(cur_dir);
                    orig_dir = self.object.original_up;
                    angle = acosd(dot(cur_dir,orig_dir));
                    axis = cross(cur_dir,orig_dir);
                    self.Rotate_Frame_AA(axis,angle);
                otherwise
                    fprintf("Object of class %s, is not a 'Helix_Object'",class(self.object));
            end
                
        end
        
        function output = Get_Translation_To_World_Frame(self)
            output = l_Recursive_Get_Translation_To_World_Frame(self, [0,0,0]');
            function output = l_Recursive_Get_Translation_To_World_Frame(a_frame, a_output)
                if (a_frame.base)
                    output = a_frame.translation_from_parent;
                else
                    output = a_output + l_Recursive_Get_Translation_To_World_Frame(a_frame.parent, a_output);
                end
            end
        end
        
        function output = Get_Rotaion_To_World_Frame(self)
            output = l_Recursive_Get_Rotation_To_World_Frame(self,eye(3));
            function output = l_Recursive_Get_Rotation_To_World_Frame(a_frame, a_output)
                if (a_frame.base)
                    output = a_frame.rotation_from_parent;
                else
                    output = a_output * l_Recursive_Get_Rotation_To_World_Frame(a_frame.parent, a_output);
                end
            end
        end
                
        function output = Get_Position_World_Frame(self,position)
            output = l_Recursive_Get_Endpoint_World_Frame(self,position);
            function output = l_Recursive_Get_Endpoint_World_Frame(a_frame,a_output)
                if a_frame.base
                    output = a_frame.position_start + a_output;
                else
                    a_output = l_Recursive_Get_Endpoint_Parent_Frame(a_frame,a_output);
                    output = l_Recursive_Get_Endpoint_World_Frame(a_frame.parent,a_output);
                end
            end
            function output = l_Recursive_Get_Endpoint_Parent_Frame(a_frame,a_output)
                output = a_frame.parent.position_end + a_frame.translation_from_parent + (a_frame.rotation_from_parent * a_output);
            end
        end
        function output = Get_World_Position_Current_Frame(self,position)
            function output = l_Recursive_Get_Position_Child_Frame(a_frame, a_output)
                if (isempty(a_frame.children))
                    output = (a_output - a_frame.translation_from_parent)/(a_frame.rotation_from_parent);
                end
            end
        end
        
        function output = Get_Endpoint_Current_Frame(self)
            output = self.position_end;
        end
        function output = Get_Endpoint_Parent_Frame(self)
            output = [0;0;0];
            if (~self.base)             
                output = self.parent.position_end + self.translation_from_parent + (self.rotation_from_parent * self.position_end);
            end
            %output = self.position_end + self.parent_frame_position + (self.parent_frame_rotation * self.position_end);
        end
        function output = Get_Endpoint_World_Frame(self)
            output = l_Recursive_Get_Endpoint_World_Frame(self,self.position_end);
            function output = l_Recursive_Get_Endpoint_World_Frame(a_frame,a_output)
                if a_frame.base
                    output = a_frame.position_start + a_output;
                else
                    a_output = l_Recursive_Get_Endpoint_Parent_Frame(a_frame,a_output);
                    output = l_Recursive_Get_Endpoint_World_Frame(a_frame.parent,a_output);
                end
            end
            function output = l_Recursive_Get_Endpoint_Parent_Frame(a_frame,a_output)
                output = a_frame.parent.position_end + a_frame.translation_from_parent + (a_frame.rotation_from_parent * a_output);
                %output = a_frame.parent_frame_position + (a_frame.parent_frame_rotation * a_output);
            end
        end
        function output = Get_Startpoint_World_Frame(self)
            output = l_Recursive_Get_Endpoint_World_Frame(self,self.position_start);
            function output = l_Recursive_Get_Endpoint_World_Frame(a_frame,a_output)
                if a_frame.base
                    output = a_frame.position_start + a_output;
                else
                    a_output = l_Recursive_Get_Endpoint_Parent_Frame(a_frame,a_output);
                    output = l_Recursive_Get_Endpoint_World_Frame(a_frame.parent,a_output);
                end
            end
            function output = l_Recursive_Get_Endpoint_Parent_Frame(a_frame,a_output)
                output = a_frame.parent.position_end + a_frame.translation_from_parent + (a_frame.rotation_from_parent * a_output);
            end
        end
        
        % Function to be used as base frame only. Defined by "CS_" =
        % "Coordinate System", at the beginning of the function.
        function output = CS_Generate_Complete_Figure(self)
            assert(self.base);
            output = [];
            for nChildren = 1:length(self.children)
                output = l_Recursive_Generate_Complete_Figure(self.children(nChildren),output);
            end
            
            function output = l_Recursive_Generate_Complete_Figure(a_frame,a_output)
                figure_data = struct("type",[],"data",[]);
                switch (class(a_frame.object))
                    case "Helix_Object"
                        % Get this object helix positions for figure generation in
                        % world frame
                        figure_data.type = "Helix_Object";
                        helixData = a_frame.object.Get_Figure_Data();
                        % Calculate Helix Metadata
                        if (isempty(a_frame.parent.parent)) % Workaround for old sims without joint objects
                            helix_index_logical = a_frame.parent.children == a_frame.parent;
                        else
                            helix_index_logical = a_frame.parent.parent.children == a_frame.parent;
                        end                        
                        index_sequence = 1:length(helix_index_logical);
                        helixData.helix_index = index_sequence(helix_index_logical);
                        for i = 1:length(helixData.positions(:,1))
                            % Deal with first helix position
                            helixPos = helixData.positions(1,:)';
                            helixPos = a_frame.Get_Position_World_Frame(helixPos);
                            helixData.positions(1,:) = helixPos;
                            
                            % Deal with second helix position
                            helixPos = helixData.positions(2,:)';
                            helixPos = a_frame.Get_Position_World_Frame(helixPos);
                            helixData.positions(2,:) = helixPos;
                        end
                        % Get world frame data for cylinder data
                        for i = 1:length(helixData.cylinder.X(1,:))
                            % Deal with first cylinder head
                            cylinderPos = [helixData.cylinder.X(1,i);helixData.cylinder.Y(1,i);helixData.cylinder.Z(1,i)];
                            cylinderPos = a_frame.Get_Position_World_Frame(cylinderPos);
                            helixData.cylinder.X(1,i) = cylinderPos(1);
                            helixData.cylinder.Y(1,i) = cylinderPos(2);
                            helixData.cylinder.Z(1,i) = cylinderPos(3);
                            
                            % Deal with second cylinder head
                            
                            cylinderPos = [helixData.cylinder.X(2,i);helixData.cylinder.Y(2,i);helixData.cylinder.Z(2,i)];
                            cylinderPos = a_frame.Get_Position_World_Frame(cylinderPos);
                            helixData.cylinder.X(2,i) = cylinderPos(1);
                            helixData.cylinder.Y(2,i) = cylinderPos(2);
                            helixData.cylinder.Z(2,i) = cylinderPos(3);
                        end
                        for i = 1:length(helixData.sidechains(:,1))
                            % Deal with sidechain positions
                            sidechainPos = helixData.sidechains(i,:)';
                            sidechainPos = a_frame.Get_Position_World_Frame(sidechainPos);
                            helixData.sidechains(i,:) = sidechainPos';
                        end
                        for i = 1:length(helixData.sidechain_mount.X)
                            mountPos = [helixData.sidechain_mount.X(i);helixData.sidechain_mount.Y(i);helixData.sidechain_mount.Z(i)];
                            mountPos = a_frame.Get_Position_World_Frame(mountPos);
                            helixData.sidechain_mount.X(i) = mountPos(1);
                            helixData.sidechain_mount.Y(i) = mountPos(2);
                            helixData.sidechain_mount.Z(i) = mountPos(3);
                        end 
                        figure_data.data = helixData;
                        output = [a_output,figure_data];
                    case "Joint_Object"
                        figure_data.type = "Joint_Object";
                        jointData.positions = [a_frame.object.position_One;a_frame.object.position_Two];
                        jointPos = jointData.positions(1,:)';
                        jointPos = a_frame.Get_Position_World_Frame(jointPos);
                        jointData.positions(1,:) = jointPos';
                        jointPos = jointData.positions(2,:)';
                        jointPos = a_frame.Get_Position_World_Frame(jointPos);
                        jointData.positions(2,:) = jointPos';
                        figure_data.data = jointData;
                        output = [a_output,figure_data];
                    otherwise
                        fprintf("ERROR: Unexpected object type.\n");
                        output = a_output;
                end
                % Recursive Search for Helices
                if ~isempty(a_frame.children) % If there are child frames. Continue Recursive Search...
                    for i = 1:length(a_frame.children)
                        output = l_Recursive_Generate_Complete_Figure(a_frame.children(i),output);
                    end
                end
            end
        end
    end
    methods (Static)
        function output = rotationMatrix(a_axis,a_angle)
            %ROTATIONMATRIX Generates a 3x3 rotaion matrix which represents a rotaion
            %of an angle 'a_angle' in degress around an axis a_axis as a unit vecotr
            %   X Axis = [1;0;0]
            %   Y Axis = [0;1;0]
            %   Z Axis = [0;0;1]
            a_axis = Vector_Helper.Make_Column(a_axis);
            
            if norm(a_axis) ~= 1
                a_axis = a_axis/norm(a_axis);
            end
            
            c = cosd(a_angle);
            s = sind(a_angle);
            t = 1-c;
            
            x = a_axis(1);
            y = a_axis(2);
            z = a_axis(3);
            
            output = [  (t*x*x)+c,      (t*x*y)-(z*s),  (t*x*z)+(y*s);
                (t*x*y)+(z*s),  (t*y*y)+c,      (t*y*z)-(x*s);
                (t*x*z)-(y*s),  (t*y*z)+(x*s),  (t*z*z)+c;];
            
        end
        function output = Generate_Rotation_from_VecOne_to_VecTwo(vecOne,VecTwo)
            vec_one = vecOne;
            vec_one = vec_one / norm(vec_one);
            vec_two = VecTwo;
            vec_cross = cross(vec_one,vec_two);
            vec_dot = dot(vec_one,vec_two);
            %             vec_skew = [0,-vec_cross(3),vec_cross(2);vec_cross(3),0,-vec_cross(1);-vec_cross(2),vec_cross(1),0];
            %             output = eye(3) + vec_skew + (vec_skew^2 * ((1 - vec_dot)/(norm(vec_cross)^2)));
            output = Coordinate_Frame.rotationMatrix(vec_cross,acosd(vec_dot));
        end
    end
end

function output = Point_Line_Intersect(p1,p2,p0)
%Calculate position along line defined by points p1 and p2 which is closest
%to point p0
    t = -(dot(p1-p0,p2-p1))/(norm(p2-p1)^2);
    output = p1 + t*(p2-p1);
end

function output = hsl2rgb(h,s,l)
    C = (1 - abs(2*l - 1)) * s;
    X = C * (1 - abs( mod(h/60, 2) - 1 ));
    M = l - C/2;
    inter_output = [0,0,0];
    if ((h >= 0 && h < 60) || h == 360)
        inter_output = [C,X,0];
    elseif (h >= 60 && h < 120)
        inter_output = [X,C,0];
    elseif (h >= 120 && h < 180)
        inter_output = [0,C,X];
    elseif (h >= 180 && h < 240)
        inter_output = [0,X,C];
    elseif (h >= 240 && h < 300)
        inter_output = [X,0,C];
    elseif (h >= 300 && h < 360)
        inter_output = [C,0,X];
    end    
    output = (inter_output + M) .* 255;
end


function PlotSphere(X,Y,Z,r,c)
    if (nargin <= 4)
        c = 'b';
    end
    [x,y,z] = sphere();
    surf(X + (x*r),Y + (y*r),Z + (z*r), 'EdgeAlpha',0,'FaceColor',c);
end