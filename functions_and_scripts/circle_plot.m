function circle_plot(ax_handle, circle_centre, circle_radius, az, el, LineSpec)
%PLOTS A CIRCLE THAT'S SURFACE IS PERPENDICULAR TO FIGURES CAMERA AXIS.
    
    if(isempty(LineSpec))
        LineSpec = 'k-';
    end

    C = circle_centre;
    R = circle_radius;
    theta = linspace(0,2*pi);
    
    X = C(1) + R * cos(theta);
    Z = C(2) + R * sin(theta);
    Y = C(3) + 0 * theta;
    
    pnts = [X;Y;Z];
    
    rotaz = RotationMatrix([0 0 1], az);
    rotel = RotationMatrix([1 0 0], -el);
    
    pnts = rotel * (rotaz * pnts);
    
    % Makes sure the circle plot does not overwrite.
    hold(ax_handle,'on');
    plot3(ax_handle,pnts(1,:),pnts(2,:),pnts(3,:),LineSpec);
    hold(ax_handle,'off');
end

