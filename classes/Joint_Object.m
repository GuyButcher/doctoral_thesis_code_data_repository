classdef Joint_Object < handle
    %JOINT_OBJECT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        parent_frame
        position_One            % [x,y,z]
        position_Two            % [x,y,z]
    end
    
    methods
        function self = Joint_Object(position)
            
            % Init Properties
            self.position_One = zeros(1,3);
            self.position_Two = zeros(1,3);
            
            % Arguement to Property assignment
            switch nargin
                case 0
                    % No input arguements, expecting position to be set
                    % later
                case 1
                    self.position_Two = Vector_Helper.Make_Row (position);
                otherwise
                    fprintf("ERROR: unexpected number of input arguements.\n");
            end
        end
        function Set_End_Position(self,position)
            self.position_Two = Vector_Helper.Make_Row(position);
            self.parent_frame.position_end = Vector_Helper.Make_Column(position);
        end
    end
end

