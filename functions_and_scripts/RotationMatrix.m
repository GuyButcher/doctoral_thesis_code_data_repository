function output = RotationMatrix(a_axis,a_angle)
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