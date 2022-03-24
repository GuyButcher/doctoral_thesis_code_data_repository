function output = apply_rotation(position, angle)
% Remeber +ev is counter clockwise and -ev is clockwise.
    rotMat = [cosd(angle),-sind(angle);sind(angle),cosd(angle)];
    output = rotMat * position;
end