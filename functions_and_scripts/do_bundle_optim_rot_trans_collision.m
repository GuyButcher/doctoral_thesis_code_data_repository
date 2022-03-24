function [results, output_image_data] = do_bundle_optim_rot_trans_collision(domain_index, helix_indeces,keep_figures)
%DO_BUNDLE_OPTIMISATION_ROT Summary of this function goes here
%   Detailed explanation goes here

    switch nargin
        case 1
            error('Function needs either two or three input arguements.');
        case 2
            keep_figures = false;
        case 3            
    end
    
    % Generate Sim Structure for base 'grey' image.

    sim = Simulation.Load_Sim_Talin_Subdomain(domain_index);
    num_helices = length(sim.helices);
    attributes = sim.Generate_Attributes_Struct();
    attributes.mode = ["dots","primary_sidechain_emphesis","hide_normal_sidechains"];
    attributes.primary_sidechain_stick_colour = 'k';
    attributes.helix_colours = ones(5,3) * 0.6627;
    fig1 = sim.Make_Figure_With_Attributes(attributes);
    view(fig1.CurrentAxes,122,12);
    light(fig1.CurrentAxes);

    % Run Optimisation to Generate

    sim = Simulation();
    helices_to_include = 1:num_helices;
    helices_to_rotate = helices_to_include;
    inital_rotation_conditions_coefficient = 0;
    inital_translation_conditions_coefficient = 0;
    bounds_rotation_coefficient = 15;
    bounds_translation_coefficient = 10;
    sim.Initialise_Optim_RotTrans_Collision_Sim(domain_index,helix_indeces,helices_to_include,helices_to_rotate);
    results = sim.Run_Optim_RotTrans_Collision_Sim(inital_rotation_conditions_coefficient,inital_translation_conditions_coefficient,bounds_rotation_coefficient,bounds_translation_coefficient,100);
    attributes.primary_sidechain_stick_colour = '';
    attributes.helix_colours = ones(5,3) * -1;

    fig2 = sim.Make_Figure_With_Attributes(attributes);
    view(fig2.CurrentAxes,122,12);
    light(fig2.CurrentAxes);

    xl = xlim(fig2.CurrentAxes);
    yl = ylim(fig2.CurrentAxes);
    zl = zlim(fig2.CurrentAxes);

    xlim(fig1.CurrentAxes, xl);
    ylim(fig1.CurrentAxes, yl);
    zlim(fig1.CurrentAxes, zl);

    imgfig1 = print(fig1,'-RGBImage','-r300');
    imgfig2 = print(fig2,'-RGBImage','-r300');
    
    if (~keep_figures)
        close(fig1);
        close(fig2);
    end

    % image manipulation for transparency overlay

    
    imgfig3 = (imgfig1*0.5) + (imgfig2 * 0.6);
    output_image_data = imgfig3;    
%     fig3 = figure();
%     imshow(fig3);    
end

