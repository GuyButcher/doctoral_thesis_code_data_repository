function data = generate_force_extension_plot_data()

    if(not(isfile('optim_collision_results_table.mat')))
        generate_results;
    end
    

    collision_results = load('results\force_sim_and_optim_collision_results.mat', 'table_data');
    collision_results = collision_results.table_data;

    
    load_digitised_experimental_data;

    unfolding_event_extension = 34.4;
    unfolding_event_force_diff = 0.0991;
    stretch_extension_per_pn = 7.6;

    % Original Simulation Force Plot

    original_force = collision_results.original_force;
    %              0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13;
    original_extension = [160,200,240,280,320,360,400,440,520,520,560,600,640,680];

    [sorted_original_forces,original_sortID] = sort(-original_force' * 10^(12));

    y_start = 160;
    x_start = -30;

    X_Original = x_start;
    Y_Original = y_start;

    for bundle_index = 1:13
        if(original_sortID(bundle_index) == 7)
            X_Original = [X_Original, sorted_original_forces(bundle_index), sorted_original_forces(bundle_index)];
            Y_Original = [Y_Original, Y_Original(end), Y_Original(end) + 80];
        elseif(original_sortID(bundle_index) == 8)
        else
            X_Original = [X_Original, sorted_original_forces(bundle_index), sorted_original_forces(bundle_index)];
            Y_Original = [Y_Original, Y_Original(end), Y_Original(end) + 40];
        end
    end

    X_Original = [X_Original, 50];
    Y_Original = [Y_Original, Y_Original(end)];

    data.original.x = X_Original;
    data.original.y = Y_Original;

    % Collision Simulation Plot

    collision_force = collision_results.optim_col_force;
    collision_extension = [160,200,240,280,320,360,400,440,520,520,560,600,640,680];

    [sorted_collision_forces,collision_sortID] = sort(-collision_force' * 10^(12));

    y_start = 160;
    x_start = 0;

    X_Collision = x_start;
    Y_Collision = y_start;

    for bundle_index = 1:13
        if(collision_sortID(bundle_index) == 7)
            collision_force_7 = sorted_collision_forces(bundle_index) + sorted_collision_forces(collision_sortID == 8);
            X_Collision = [X_Collision, collision_force_7, collision_force_7];
            Y_Collision = [Y_Collision, Y_Collision(end), Y_Collision(end) + 80];
        elseif(collision_sortID(bundle_index) == 8)
        else
            X_Collision = [X_Collision, sorted_collision_forces(bundle_index), sorted_collision_forces(bundle_index)];
            Y_Collision = [Y_Collision, Y_Collision(end), Y_Collision(end) + 40];
        end
    end

    X_Collision = [X_Collision,50];
    Y_Collision = [Y_Collision, Y_Collision(end)];

    data.collision.x = X_Collision;
    data.collision.y = Y_Collision;

    % Yoa Graph Plot with stretch params

    yao_force = experimental_results.force;
    yao_index = 1:13;
    yao_extension = experimental_results.extension;

    unfolding_event_extension = 34.4;
    unfolding_event_extension = 39;
    unfolding_event_force_diff = 0.0991 / 2;
    stretch_extension_per_pn = 7.6 / 2;

    y_start = 160;
    x_start = 0;

    X_Yao = x_start;
    Y_Yao = y_start;

    for bundle_index = 1:13
        if(yao_index(bundle_index) == 7)
            force_diff = yao_force(bundle_index) - X_Yao(end);
            X_Yao = [X_Yao, yao_force(bundle_index), (yao_force(bundle_index) + unfolding_event_force_diff)];
            y_end = Y_Yao(end) + (force_diff * stretch_extension_per_pn);
            Y_Yao = [Y_Yao, y_end, y_end + (2 * unfolding_event_extension)];
        elseif(yao_index(bundle_index) == 8)
        elseif(yao_index(bundle_index) == 1)
            force_diff = yao_force(bundle_index) - X_Yao(end);
            X_Yao = [X_Yao, yao_force(bundle_index), (yao_force(bundle_index) + unfolding_event_force_diff)];
            y_end = Y_Yao(end) + (force_diff * stretch_extension_per_pn);
            Y_Yao = [Y_Yao, y_end, y_end + (unfolding_event_extension / 2)];
        else
            force_diff = yao_force(bundle_index) - X_Yao(end);
            X_Yao = [X_Yao, yao_force(bundle_index), (yao_force(bundle_index) + unfolding_event_force_diff)];
            y_end = Y_Yao(end) + (force_diff * stretch_extension_per_pn);
            Y_Yao = [Y_Yao, y_end, y_end + unfolding_event_extension];
        end
    end

    force_diff = 30 - yao_force(end);
    X_Yao = [X_Yao, 30];
    Y_Yao = [Y_Yao, Y_Yao(end) + (force_diff * stretch_extension_per_pn)];

    data.yao_stretch.x = X_Yao;
    data.yao_stretch.y = Y_Yao;

end

