function output = Point_Line_Intersection(p1,p2,p0)
%Calculate position along line defined by points p1 and p2 which is closest
%to point p0
    t = -(dot(p1-p0,p2-p1))/(norm(p2-p1)^2);
    output = p1 + t*(p2-p1);
end

