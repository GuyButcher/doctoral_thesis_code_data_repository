function [forces, distance_output] = do_seperation_sim(h1Length, h2Length, distanceLimit)
    % Two Unfolding Simple Helix Models: rotation at helix base

    distance_output = 0:0.5:distanceLimit;

    % Helix Models
    force_constant = 10;
    distance_interval = 1;
    distance_limit = distanceLimit;
    
    h1length = h1Length;
    h1positions = zeros(2,h1length);
    h1positions(2,:) = 1:h1length;
    h1rot = 0;
    
    h1offset = [-2 2]';
    h1offsetRot = 0;
    
    h2length = h2Length;
    h2positions = zeros(2,h2length);
    h2positions(2,:) = 1:h2length;
    h2rot = 0;
    
    h2offset = [2 2]';
    h2offsetRot = 0;
    
    % Possible Sidechain Interactions
    if h1length == h2length
        interaction_indices = nchoosek(1:h1length,2);
        for i = 1:h1length
            interaction_indices = [interaction_indices; [i i]];
        end
    elseif h1length > h2length
        interaction_indices = nchoosek(1:h1length,2);
        for i = 1:h1length
            interaction_indices = [interaction_indices; [i i]];
        end
        
        % Remove indices that are too large for smaller helix
        % Smaller helix is vertical index 1
        index = interaction_indices(:,1) > h2length;
        interaction_indices(index,:) = [];
        % As h2 is smaller, put on right of list
        interaction_indices = fliplr(interaction_indices);
    else % h2length > h1length
        interaction_indices = nchoosek(1:h2length,2);
        for i = 1:h2length
            interaction_indices = [interaction_indices; [i i]];
        end
        
        % Remove indices that are too large for smaller helix
        % Smaller helix is vertical index 1
        index = interaction_indices(:,1) > h1length;
        interaction_indices(index,:) = [];
    end
    
    % Update Current Sidechain Positions
    h1current = h1offset + apply_rotation(h1positions,h1rot);
    h1current = apply_rotation(h1current,h1offsetRot);
    h2current = h2offset + apply_rotation(h2positions,h2rot);
    h2current = apply_rotation(h2current,h2offsetRot);
    
    % Do Force calculation
    forces = distance_output * 0;
    for k = 1:length(interaction_indices)
        pos1 = h1current(:,interaction_indices(k,1));
        pos2 = h2current(:,interaction_indices(k,2));
        distance = norm(pos2 - pos1);
        force = force_constant/(distance^2);
        forces(1) = forces(1) + force;
    end
    
    for i = 2:length(distance_output)
        % Update Sim
        displacement = distance_output(i)/2;
        h1offset(1) = -2 - displacement;
        h2offset(1) = 2 + displacement;
        % Update Current Sidechain Positions
        h1current = h1offset + apply_rotation(h1positions,h1rot);
        h1current = apply_rotation(h1current,h1offsetRot);
        h2current = h2offset + apply_rotation(h2positions,h2rot);
        h2current = apply_rotation(h2current,h2offsetRot);
        % Do Force calculation
        for k = 1:length(interaction_indices)
            pos1 = h1current(:,interaction_indices(k,1));
            pos2 = h2current(:,interaction_indices(k,2));
            distance = norm(pos2 - pos1);
            force = force_constant/(distance^2);
            forces(i) = forces(i) + force;
        end
    end
end