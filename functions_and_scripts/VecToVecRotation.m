function output = VecToVecRotation(vecOne,VecTwo)
    %VECTOVECROTAION Generates a 3x3 rotation matrix which transforms and
    %axis from the direction of vecOne to that of vecTwo
    %The output rotation matrix functions when applied to column vectors.
    vec_one = vecOne;
    vec_one = vec_one / norm(vec_one);
    vec_two = VecTwo;
    vec_cross = cross(vec_one,vec_two);
    vec_dot = dot(vec_one,vec_two);
    output = Coordinate_Frame.rotationMatrix(vec_cross,acosd(vec_dot));
end