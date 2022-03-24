classdef Vector_Helper < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Static)
        function output = Make_Column(vector)
            if (class(vector) ~= 'double')
                error('input vector is not a number.')
            end
            dimensions = size(vector);
            if (~all(dimensions == [3 1]) && ~all(dimensions == [1 3]))
                error('vector is not column or row vector of lenght 3.')
            end
            output = [vector(1);vector(2);vector(3)];
        end
        function output = Make_Row(vector)
            if (class(vector) ~= 'double')
                error('input vector is not a number.')
            end
            dimensions = size(vector);
            if (~all(dimensions == [3 1]) && ~all(dimensions == [1 3]))
                error('vector is not column or row vector of lenght 3.')
            end
            output = [vector(1),vector(2),vector(3)];
        end
    end
end

