classdef Helix_Object < handle
    %HELIX_OBJECT Summary of this class goes here
    %   Detailed explanation goes here
  
    properties (Access = public)
        parent_frame
        position_One            % [x,y,z]
        position_Two            % [x,y,z]
        radius                  % r
        up_direction            % 
        
        original_pos            % [x,y,z]
        original_ori            % [x,y,z] unit vector
        original_up             % []
        original_direction      % [x,y,z]
        
        sidechain_names
        sidechain_positions     % [x,y,z];[x,y,z];[x,y,z];...
        sidechain_direction     % [x,y,z] unit vector
        sidechain_sizes         % radius (approximated spherical volume)
        sidechain_properties    % Mass,Charge,Hydropathy
    end
    
    methods
        % New function created on 22 sep 2019
        function self = Helix_Object(backbone_positions, sidechain_positions, sidechain_sizes)
        
            % backbone_positions: N-by-3 [[x,y,z];[...];...]
            % sidechain_Positions: N-by-3 [[x,y,z];[...];...]
            
            % Find central axis through helix.
            midpoint = mean(backbone_positions,1);
            residuals = bsxfun(@minus,backbone_positions,midpoint);
            C = (residuals'*residuals)/(length(residuals(:,1))-1);
            [R,~] = svd(C,0);
            
            x = residuals*R(:,1);
            xMin = min(x);
            xMax = max(x);
            Xa = xMin*R(:,1)' + midpoint;
            Xb = xMax*R(:,1)' + midpoint;
            helixVec = Xb - Xa;
            helixUnitVec = helixVec / norm(helixVec);
            
            % Find correct orientaion of helixUnitVec
            % Get pos of first and last sidechain
            first_sidechain_pos = [0,0,0];
            for i = 1:length(sidechain_positions)
                first_sidechain_pos = sidechain_positions(i,:);
                if (~all(first_sidechain_pos == 0))
                    break
                end
            end
            last_sidechain_pos = [0,0,0];
            for i = length(sidechain_positions):-1:1
                last_sidechain_pos = sidechain_positions(i,:);
                if (~all(last_sidechain_pos == 0))
                    break
                end
            end
            % create point along helixDir in positive direction
            direction_test_vec = helixUnitVec*(5) + midpoint;
            % calc distances from point to first and last sidechain pos
            length_to_first = norm(first_sidechain_pos - direction_test_vec);
            length_to_last = norm(last_sidechain_pos - direction_test_vec);
            % Shortest distance is the current positive direction
            % Want +ve dir to be from first towards last. Therefore, if
            % shortest distance is towards first sidechain. Invert
            % helixUnitVec
            if (length_to_first < length_to_last)
                helixUnitVec = -helixUnitVec;
            end
            
            % set pos one and two of helix cylinder.
            % Pos One is the start of the cylinder.
            % Pos Two is the end of the cylinder. 
            % Thus direction is One -> Two
            helixLength = norm(backbone_positions(1,:) - backbone_positions(end,:));
            helixLength = norm((helixUnitVec*(helixLength/2) + midpoint) - midpoint);
            for i = 1:length(sidechain_positions(:,1))
                if(~all(sidechain_positions(i,:) == 0))
                    curLength = norm(sidechain_positions(i,:) - midpoint);
                    if (curLength > helixLength)
                        helixLength = curLength;
                    end
                end
            end
            self.position_One = helixUnitVec*(-helixLength) + midpoint;
            self.position_Two = helixUnitVec*(helixLength) + midpoint;
            
            % Find the radius for the cylinder
            radius = zeros(1,length(backbone_positions(:,1)));
            for i = 1:length(backbone_positions(:,1))
                a = self.position_One - self.position_Two;
                b = backbone_positions(i,:) - self.position_Two;
                radius(i) = norm(cross(a,b))/norm(a);
            end
            self.radius = sum(radius)/length(radius);
            
            % Store sidechain positions
            self.sidechain_positions = sidechain_positions;
            self.sidechain_sizes = sidechain_sizes;
            
            % Save the original starting position of the helix and its
            % direction
            self.original_direction = helixUnitVec;
            self.original_pos = self.position_One;
            
            % Calculate and save helix up direction
            % Defined as perpendicular direction away from the central
            % helix axis that passes through first non-glycine sidechain
            int_point = Point_Line_Intersect(self.position_One,self.position_Two,first_sidechain_pos);
            dir_vec = first_sidechain_pos - int_point;
            dir_unit_vec = dir_vec / norm(dir_vec);
            self.original_up = dir_unit_vec;
        end
        function Properties(self,properties)
            self.sidechain_properties = properties;
        end
        function Save_Original_Pos_Ori(self,pos,ori)
            self.original_pos = pos;
            self.original_ori = ori;
        end
        
        function Generate_Figure(self)
            N = 100;
            theta = linspace(0,2*pi,N);

            helixUnitVec = self.position_Two - self.position_One;
            helixOrthoganol = rand(1,3);
            
            VecA = helixUnitVec; % Unit Vector of Helix Axis
            VecB = cross(VecA,helixOrthoganol);
            VecB = VecB/norm(VecB); % First Orthogonal Vector to Helix Axis
            VecC = cross(VecA,VecB);
            VecC = VecC/norm(VecC); % Second Orthogonal Vector, perpendicular to the axis and VecB
            
            X = zeros(2,N); % initialise vars
            Y = zeros(2,N);
            Z = zeros(2,N);
            
            for i = 0:1 % Create Surface Data for Cylinder
                X(i+1,:) = self.position_One(1) + (self.position_Two(1) - self.position_One(1))*i + self.radius*cos(theta)*VecB(1) + self.radius*sin(theta)*VecC(1);
                Y(i+1,:) = self.position_One(2) + (self.position_Two(2) - self.position_One(2))*i + self.radius*cos(theta)*VecB(2) + self.radius*sin(theta)*VecC(2);
                Z(i+1,:) = self.position_One(3) + (self.position_Two(3) - self.position_One(3))*i + self.radius*cos(theta)*VecB(3) + self.radius*sin(theta)*VecC(3);
            end
            
            % Code to handle Sidechains
            sidechains = self.sidechain_positions;
            count = 0;
            for i = 1:length(sidechains(:,1))
                if all(sidechains(i-count,:) == [0,0,0])
                    sidechains(i-count,:) = [];
                    count = count + 1;
                end
            end
            
            sidechain_mounts = zeros(length(sidechains(:,1)),3);
            for i = 1:length(sidechains(:,1))
                sidechain_mounts(i,:) = Point_Line_Intersect(self.position_One,self.position_Two,sidechains(i,:));
            end
            
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
            for i = 1:length(sidechains(:,1)) % Sidechain Spheres
                PlotSphere(sidechains(i,1),sidechains(i,2),sidechains(i,3),self.sidechain_sizes(i));
            end
            plot3(sX,sY,sZ,'Color','green');
            plot3(sX(:,1),sY(:,1),sZ(:,1),'Color', 'ma');
            if (1)
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
            end
            hold off
        end
        function Generate_Figure_Thick_Lines(self)
            N = 100;
            theta = linspace(0,2*pi,N);

            helixUnitVec = self.position_Two - self.position_One;
            helixOrthoganol = rand(1,3);
            
            VecA = helixUnitVec; % Unit Vector of Helix Axis
            VecB = cross(VecA,helixOrthoganol);
            VecB = VecB/norm(VecB); % First Orthogonal Vector to Helix Axis
            VecC = cross(VecA,VecB);
            VecC = VecC/norm(VecC); % Second Orthogonal Vector, perpendicular to the axis and VecB
            
            X = zeros(2,N); % initialise vars
            Y = zeros(2,N);
            Z = zeros(2,N);
            
            for i = 0:1 % Create Surface Data for Cylinder
                X(i+1,:) = self.position_One(1) + (self.position_Two(1) - self.position_One(1))*i + self.radius*cos(theta)*VecB(1) + self.radius*sin(theta)*VecC(1);
                Y(i+1,:) = self.position_One(2) + (self.position_Two(2) - self.position_One(2))*i + self.radius*cos(theta)*VecB(2) + self.radius*sin(theta)*VecC(2);
                Z(i+1,:) = self.position_One(3) + (self.position_Two(3) - self.position_One(3))*i + self.radius*cos(theta)*VecB(3) + self.radius*sin(theta)*VecC(3);
            end
            
            % Code to handle Sidechains
            sidechains = self.sidechain_positions;
            count = 0;
            for i = 1:length(sidechains(:,1))
                if all(sidechains(i-count,:) == [0,0,0])
                    sidechains(i-count,:) = [];
                    count = count + 1;
                end
            end
            
            sidechain_mounts = zeros(length(sidechains(:,1)),3);
            for i = 1:length(sidechains(:,1))
                sidechain_mounts(i,:) = Point_Line_Intersect(self.position_One,self.position_Two,sidechains(i,:));
            end
            
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
            plot3(sX,sY,sZ,'Color','green','LineWidth',1);
            plot3(sX(:,1),sY(:,1),sZ(:,1),'Color', 'ma','LineWidth',2);
            if (1)
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
            end
            hold off
        end
        function output = Get_Figure_Data(self)
            N = 100;
            theta = linspace(0,2*pi,N);

            helixUnitVec = self.position_Two - self.position_One;
            helixOrthoganol = rand(1,3);
            
            VecA = helixUnitVec; % Unit Vector of Helix Axis
            VecB = cross(VecA,helixOrthoganol);
            VecB = VecB/norm(VecB); % First Orthogonal Vector to Helix Axis
            VecC = cross(VecA,VecB);
            VecC = VecC/norm(VecC);
            
            X = zeros(2,N); % initialise vars
            Y = zeros(2,N);
            Z = zeros(2,N);
            
            for i = 0:1 % Create Surface Data for Cylinder
                X(i+1,:) = self.position_One(1) + (self.position_Two(1) - self.position_One(1))*i + self.radius*cos(theta)*VecB(1) + self.radius*sin(theta)*VecC(1);
                Y(i+1,:) = self.position_One(2) + (self.position_Two(2) - self.position_One(2))*i + self.radius*cos(theta)*VecB(2) + self.radius*sin(theta)*VecC(2);
                Z(i+1,:) = self.position_One(3) + (self.position_Two(3) - self.position_One(3))*i + self.radius*cos(theta)*VecB(3) + self.radius*sin(theta)*VecC(3);
            end
            
            % Code to handle Sidechains
            sidechains = self.sidechain_positions;
            count = 0;
            for i = 1:length(sidechains(:,1))
                if all(sidechains(i-count,:) == [0,0,0])
                    sidechains(i-count,:) = [];
                    count = count + 1;
                end
            end
            
            sidechain_mounts = zeros(length(sidechains(:,1)),3);
            for i = 1:length(sidechains(:,1))
                sidechain_mounts(i,:) = Point_Line_Intersect(self.position_One,self.position_Two,sidechains(i,:));
            end
            
            sX = [sidechain_mounts(:,1)';sidechains(:,1)'];
            sY = [sidechain_mounts(:,2)';sidechains(:,2)'];
            sZ = [sidechain_mounts(:,3)';sidechains(:,3)'];
            
            output.positions = [self.position_One;self.position_Two];
            output.cylinder.X = X;
            output.cylinder.Y = Y;
            output.cylinder.Z = Z;
            output.sidechains = sidechains;
            output.sidechain_sizes = self.sidechain_sizes;
            output.sidechain_hydropathy = self.sidechain_properties(:,3);
            output.sidechain_mount.X = sidechain_mounts(:,1);
            output.sidechain_mount.Y = sidechain_mounts(:,2);
            output.sidechain_mount.Z = sidechain_mounts(:,3);
        end
        
        function output = Get_Position(self)
%             output.x = [self.position_One(1) self.position_Two(1)];
%             output.y = [self.position_One(2) self.position_Two(2)];
%             output.z = [self.position_One(3) self.position_Two(3)];
            output = [self.position_One;self.position_Two];
        end
        function output = Get_Rotation(self)
        end
        
        function output = Get_Position_World_Frame(self, position)
            output = self.parent_frame.Get_Position_World_Frame(position);
        end
        function output = Get_Sidechain_Position_Helix_Frame(self,index)
            if all(size(index) == [1,1]) % If index is single value
                sidechain = self.sidechain_positions(index,:);
                output = sidechain - self.position_One;
            elseif all(size(index) == [1,2]) || all(size(index) == [2,1]) % If index is N-by-2 or 2-by-N vector
                sidechains = self.sidechain_positions(index(1):index(2),:);
                output = sidechains - self.position_One;
            else
                fprintf("Error:Helix_Object:Get_Sidechain_Position_Helix_frame");
            end
            
        end
        function output = Get_Sidechain_Position(self,index)
            if all(size(index) == [1,1]) % If index is single value
                output = self.sidechain_positions(index,:);
            elseif all(size(index) == [1,2]) || all(size(index) == [2,1]) % If index is N-by-2 or 2-by-N vector
                output = self.sidechain_positions(index(1):index(2),:);
            else
                fprintf("Error:Helix_Object:Get_Sidechain_Position");
            end
        end
        function output = Get_Sidechain_Properties(self,property)
        end
        
        function output = Get_Electrostatic_Approximation(self)
        end
        
        function output = Get_Number_Sidechains(self)
            output = length(self.sidechain_properties(:,1));
        end
    end
    methods (Access = protected)
        function Update(self)
        end
    end
end

function output = Point_Line_Intersect(p1,p2,p0)
%Calculate position along line defined by points p1 and p2 which is closest
%to point p0
    t = -(dot(p1-p0,p2-p1))/(norm(p2-p1)^2);
    output = p1 + t*(p2-p1);
end

function PlotSphere(X,Y,Z,r,c)
    if (nargin < 5)
        c = 'b';
    end
    [x,y,z] = sphere();
    surf(X + (x*r),Y + (y*r),Z + (z*r), 'EdgeAlpha',0,'FaceColor',c);
end

