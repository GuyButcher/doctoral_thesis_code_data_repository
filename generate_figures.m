%% Loading and Generating Data if not already done.

startup;

% Optimisation Force Data
if (not(isfile('results/force_sim_and_optim_collision_results.mat')) || ...
    not(isfile('results/hydrophobic_results.mat')) || ...
    not(isfile('results/backbone_validation_results.mat')))
    generate_results;
    clear;
    clc;
end

table_data = load('results\force_sim_and_optim_collision_results.mat');
table_data = table_data.table_data;

clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

%% Digitise Data from Yao et al.

data = generate_force_extension_plot_data();

fig = figure();
plot(data.yao_stretch.x,data.yao_stretch.y, 'Color','Blue');

xlabel("Force, pN");
ylabel("Extension, nm");
grid(fig.CurrentAxes, 'on');
grid(fig.CurrentAxes, 'minor');
legend(fig.CurrentAxes, 'Yao et al.', 'Location', 'southeast');

exportgraphics(fig.CurrentAxes, 'figures/yao_force_extension_digitised.pdf',"ContentType","vector");
close(fig);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

%% Hybrid Spring Unfolding Figures

force_threshold = 15; % pN
forces = 0:1:40;
f_length = length(forces);

extension = zeros(1,f_length);
e1 = zeros(1,f_length);
e2 = zeros(1,f_length);

for i = 1:f_length
    [extension(i), e1(i), e2(i)] = subDomainSpring(i,15);
end

fig1 = figure();
plot(forces,extension);
xlim(fig1.CurrentAxes, [-4 44]);
ylim(fig1.CurrentAxes, [-4.5 49.5]);

xlabel('Force, pN');
ylabel('Extension, nm');
exportgraphics(fig1.CurrentAxes,'figures/hybrid_spring_single_domain_extetnsion.pdf',"ContentType","vector");

fig2 = figure();
h1=subplot(2,1,1);
plot(forces,e1);
xlim(h1, [-4 44]);
ylim(h1, [0 3.3]);

h2=subplot(2,1,2);
plot(forces,e2);
xlim(h2, [-4 44]);
ylim(h2, [-4 44]);

xlabel('Force, pN');

p1 = get(h1,'position');
p2 = get(h2,'position');
height = p1(2)+p1(4)-p2(2);
h3=axes('position',[p2(1) p2(2) p2(3) height],'Visible','off');
h_label = ylabel('Extension, nm','Visible','on');
exportgraphics(fig2,'figures/hydrid_sping_dual_extension_subplot.pdf',"ContentType","vector");

testThresholds = [5, 10, 11, 12, 12.5, 14, 15, 15, 16, 17.5, 19, 21, 24];
forces = 0:0.1:40;
extension = zeros(1,size(forces,2));
extension = force_domain_unfolding(forces,testThresholds);
data = generate_force_extension_plot_data();
fig3 = figure();
axes(fig3);
plot(fig3.CurrentAxes,forces, extension + 160)

xlabel("Force, pN");
ylabel("Extension, nm");

hold on
plot(fig3.CurrentAxes,data.yao_stretch.x,data.yao_stretch.y);
legend({'Hybrid Spring Model','Yao Data'}, "Location",'nw')
hold off
xlim(fig3.CurrentAxes, [0 30]);
exportgraphics(fig3.CurrentAxes,'figures/hydrid_spring_rod_domain_model.pdf',"ContentType","vector");

close([fig1 fig2 fig3]);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

%% Simple Helix Model Seperating Figures

[forces_one, distance_output] = do_seperation_sim(25,25,25);
[forces_two, distance_output] = do_seperation_sim(25,30,25);
[forces_three, distance_output] = do_seperation_sim(25,20,25);
[forces_four, distance_output] = do_seperation_sim(30,30,25);
[forces_five, distance_output] = do_seperation_sim(20,20,25);
[forces_six, distance_output] = do_seperation_sim(30,20,25);

fig1 = figure('PaperPositionMode', 'auto');
axes(fig1);
hold(fig1.CurrentAxes,'on');
plot(fig1.CurrentAxes, distance_output, forces_one);
plot(fig1.CurrentAxes, distance_output, forces_two);
plot(fig1.CurrentAxes, distance_output, forces_three);
plot(fig1.CurrentAxes, distance_output, forces_four);
plot(fig1.CurrentAxes, distance_output, forces_five);
plot(fig1.CurrentAxes, distance_output, forces_six);
hold(fig1.CurrentAxes,'off');
legend(fig1.CurrentAxes, {'a','b','c','d','e','f'}, "Location","ne");
xlabel(fig1.CurrentAxes,'Distance Between Helix Models');
ylabel(fig1.CurrentAxes,'Interaction Force');
exportgraphics(fig1.CurrentAxes, 'figures/simple_model_seperating_overlay.pdf','ContentType',"vector");

close(fig1);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

%% Simple Helix Model Unfolding Figures

% Two Unfolding Simple Helix Models: rotation at helix base

% Helix Models
helix_lengths = 25;

h1length = helix_lengths;
h1positions = zeros(2,h1length);
h1positions(2,:) = 1:h1length;
h1rot = 0;

h1offset = [-2 2]';
h1offsetRot = 0;

h2length = helix_lengths;
h2positions = zeros(2,h2length);
h2positions(2,:) = 1:h2length;
h2rot = 0;

h2offset = [2 2]';
h2offsetRot = 0;

force_constant = 10;

% Possible Sidechain Interactions
interaction_indices = nchoosek(1:helix_lengths,2);
for i = 1:helix_lengths
    interaction_indices = [interaction_indices; [i i]];
end

% Update Current Sidechain Positions
h1current = h1offset + apply_rotation(h1positions,h1rot);
h1current = apply_rotation(h1current,h1offsetRot);
h2current = h2offset + apply_rotation(h2positions,h2rot);
h2current = apply_rotation(h2current,h2offsetRot);

% Do Force calculation
forces = zeros(90 + 1,1);
for k = 1:length(interaction_indices)
    pos1 = h1current(:,interaction_indices(k,1));
    pos2 = h2current(:,interaction_indices(k,2));
    distance = norm(pos2 - pos1);
    force = force_constant/(distance^2);
    forces(1) = forces(1) + force;
end

for i = 1:90
    % Update Sim
    h1rot = i;
    h2rot = -i;
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
        forces(1 + i) = forces(1 + i) + force;
    end
end

forces_one = forces;

% Two Unfolding Simple Helix Models: rotation at joint base

% Helix Models
helix_lengths = 25;

h1length = helix_lengths;
h1positions = zeros(2,h1length);
h1positions(2,:) = 1:h1length;
h1rot = 0;

h1offset = [-2 2]';
h1offsetRot = 0;

h2length = helix_lengths;
h2positions = zeros(2,h2length);
h2positions(2,:) = 1:h2length;
h2rot = 0;

h2offset = [2 2]';
h2offsetRot = 0;

force_constant = 10;

% Possible Sidechain Interactions
interaction_indices = nchoosek(1:helix_lengths,2);
for i = 1:helix_lengths
    interaction_indices = [interaction_indices; [i i]];
end

% Update Current Sidechain Positions
h1current = h1offset + apply_rotation(h1positions,h1rot);
h1current = apply_rotation(h1current,h1offsetRot);
h2current = h2offset + apply_rotation(h2positions,h2rot);
h2current = apply_rotation(h2current,h2offsetRot);

% Do Force calculation
forces = zeros(90 + 1,1);
for k = 1:length(interaction_indices)
    pos1 = h1current(:,interaction_indices(k,1));
    pos2 = h2current(:,interaction_indices(k,2));
    distance = norm(pos2 - pos1);
    force = force_constant/(distance^2);
    forces(1) = forces(1) + force;
end

for i = 1:90
    % Update Sim
    h1offsetRot = i;
    h2offsetRot = -i;
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
        forces(1 + i) = forces(1 + i) + force;
    end
end

forces_two = forces;

% Two Unfolding Simple Helix Models: rotation at helix and joint base

% Helix Models
helix_lengths = 25;

h1length = helix_lengths;
h1positions = zeros(2,h1length);
h1positions(2,:) = 1:h1length;
h1rot = 0;

h1offset = [-2 2]';
h1offsetRot = 0;

h2length = helix_lengths;
h2positions = zeros(2,h2length);
h2positions(2,:) = 1:h2length;
h2rot = 0;

h2offset = [2 2]';
h2offsetRot = 0;

force_constant = 10;

% Possible Sidechain Interactions
interaction_indices = nchoosek(1:helix_lengths,2);
for i = 1:helix_lengths
    interaction_indices = [interaction_indices; [i i]];
end

% Update Current Sidechain Positions
h1current = h1offset + apply_rotation(h1positions,h1rot);
h1current = apply_rotation(h1current,h1offsetRot);
h2current = h2offset + apply_rotation(h2positions,h2rot);
h2current = apply_rotation(h2current,h2offsetRot);

% Do Force calculation
forces = zeros(90 + 1,1);
for k = 1:length(interaction_indices)
    pos1 = h1current(:,interaction_indices(k,1));
    pos2 = h2current(:,interaction_indices(k,2));
    distance = norm(pos2 - pos1);
    force = force_constant/(distance^2);
    forces(1) = forces(1) + force;
end

for i = 1:90
    % Update Sim
    h1rot = i/2;
    h2rot = -i/2;
    h1offsetRot = i/2;
    h2offsetRot = -i/2;
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
        forces(1 + i) = forces(1 + i) + force;
    end
end

forces_three = forces;

% Figure

fig_c = figure();
axes(fig_c);
hold(fig_c.CurrentAxes,"on");
plot(fig_c.CurrentAxes,0:90, forces_two);
plot(fig_c.CurrentAxes,0:90, forces_one);
plot(fig_c.CurrentAxes,0:90, forces_three);
hold(fig_c.CurrentAxes,"off");
legend(fig_c.CurrentAxes, {'a','b','c'});
xlabel(fig_c.CurrentAxes,'Angle Between Helix Models, degrees');
ylabel(fig_c.CurrentAxes,'Interaction Force');
exportgraphics(fig_c.CurrentAxes,'figures/unfolding_overlayed.pdf',"ContentType","vector");

close(fig_c);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

%% Amino-Acid Structure Figure

sim = Simulation.Load_Sim_Talin_Subdomain(3);
sim.chain.selected_data = [];
sim.chain.Make_Figure();
fig1 = gca();
view(fig1,-160,-70);
axis off
exportgraphics(fig1,'figures/FullAminoAcidChainWithLabels.pdf',"ContentType","vector");

close(fig1.Parent);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

%% Example Simulation Set-up with R3 Helices 1 and 2

startup;

sim = Simulation();
sim.Load_Amino_Acid_PDB_File(pdb_filepaths(3));
sim.Load_Helix(3,1,helix_indeces);
sim.Load_Helix(3,2,helix_indeces);
sim.Load_PDB_Parallel_Chain();

force = sim.Run_Static_Force(1,2);

attrib = Simulation.Generate_Attributes_Struct();
attrib.mode = 'simple_primary_sidechain_emphesis';
fig2 = sim.Make_Figure_With_Attributes(attrib);

view(fig2.CurrentAxes, [78.1 20.0]);
print(fig2,'figures/two_ah_r3','-dpng', '-r700');

close(fig2);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

%% Backbone Approximation Validation Figures

% Load Validation Data

if(not(isfile('results\backbone_validation_results.mat')))
    generate_results;
end

gradients = load('results\backbone_validation_results.mat');
gradients = gradients.gradients;

indices = [
    1,1; 1,2; 1,3; 1,4; 1,5;
    2,1; 2,2; 2,3; 2,4;
    3,1; 3,2; 3,3; 3,4;
    4,1; 4,2; 4,3; 4,4;
    5,1; 5,2; 5,3; 5,4; 5,5;
    6,1; 6,2; 6,3; 6,4; 6,5;
    7,1; 7,2; 7,3; 7,4; 7,5;
    8,1; 8,2; 8,3; 8,4;
    9,1; 9,2; 9,3; 9,4; 9,5;
    10,1; 10,2; 10,3; 10,4; 10,5;
    11,1; 11,2; 11,3; 11,4; 11,5;
    12,1; 12,2; 12,3; 12,4; 12,5;
    13,1; 13,2; 13,3; 13,4; 13,5;];

% Meta Figures

fig1 = figure();
axes(fig1);
plot(fig1.CurrentAxes,1:length(gradients), sort(abs(gradients)), 'LineWidth',2);
mean_gradient = mean(abs(gradients));
hold(fig1.CurrentAxes,"on");
plot(fig1.CurrentAxes,1:length(gradients),ones(1,length(gradients))*mean_gradient, 'LineWidth',1, "LineStyle","--");
hold(fig1.CurrentAxes,"off");
xlabel(fig1.CurrentAxes,'Helix Index');
ylabel(fig1.CurrentAxes,'Gradient Magnitude, degrees');
axis(fig1.CurrentAxes,'tight');
legend(fig1.CurrentAxes,{'Gradients','Mean Line'},'Location','nw');
median_gradient = median(abs(gradients));
exportgraphics(fig1.CurrentAxes,'figures/SkewGradientDistributionPlot.pdf',"ContentType","vector");

fig2 = figure();
axes(fig2);
histfit(fig2.CurrentAxes,gradients,13);
xlabel(fig2.CurrentAxes,'Gradient, degrees');
ylabel(fig2.CurrentAxes,'Count');
axis(fig2.CurrentAxes,'tight');
legend(fig2.CurrentAxes,{'Gradient Count','Normal Fit'},'Location','nw')
exportgraphics(fig2.CurrentAxes,'figures/SkewGradientNormalPlot.pdf',"ContentType","vector");

fig3 = figure();
axes(fig3);
boxplot(fig3.CurrentAxes,abs(gradients));
ylabel(fig3.CurrentAxes,'Gradients, degrees');
fig3.CurrentAxes.XTickLabel{1} = '';
exportgraphics(fig3.CurrentAxes,'figures/SkewGradientBoxplot.pdf',"ContentType","vector");

% Figures

[~,gradients_sort_abs_min_to_max] = sort(abs(gradients));
grad_domain_helix = horzcat(gradients(gradients_sort_abs_min_to_max),indices(gradients_sort_abs_min_to_max,:));
table(grad_domain_helix(:,1), grad_domain_helix(:,2), grad_domain_helix(:,3), 'VariableNames',{'gradient','domain','helix'});

[fig71_1, fig71_2] = bp_atoms_dist(7,1,8);
title(fig71_1.CurrentAxes,'');
exportgraphics(fig71_1.CurrentAxes,'figures/bb_dist_scatter_g_7-1.pdf',"ContentType","vector");
title(fig71_2.CurrentAxes,'');
exportgraphics(fig71_2.CurrentAxes,'figures/bb_dist_normal_g_7-1.pdf',"ContentType","vector");
grad1 = bp_atoms_dist_grad(7,1);
[fig122_1, fig122_2] = bp_atoms_dist(12,2,7);
title(fig122_1.CurrentAxes,'');
exportgraphics(fig122_1.CurrentAxes,'figures/bb_dist_scatter_g_12-2.pdf',"ContentType","vector");
title(fig122_2.CurrentAxes,'');
exportgraphics(fig122_2.CurrentAxes,'figures/bb_dist_normal_g_12-2.pdf',"ContentType","vector");
grad2 = bp_atoms_dist_grad(12,2);
[fig92_1, fig92_2] = bp_atoms_dist(9,2);
title(fig92_1.CurrentAxes,'');
exportgraphics(fig92_1.CurrentAxes,'figures/bb_dist_scatter_b_9-2.pdf',"ContentType","vector");
title(fig92_2.CurrentAxes,'');
exportgraphics(fig92_2.CurrentAxes,'figures/bb_dist_normal_b_9-2.pdf',"ContentType","vector");
grad3 = bp_atoms_dist_grad(9,2);
[fig13_1, fig13_2] = bp_atoms_dist(1,3);
title(fig13_1.CurrentAxes,'');
exportgraphics(fig13_1.CurrentAxes,'figures/bb_dist_scatter_b_1-3.pdf',"ContentType","vector");
title(fig13_2.CurrentAxes,'');
exportgraphics(fig13_2.CurrentAxes,'figures/bb_dist_normal_b_1-3.pdf',"ContentType","vector");
grad4 = bp_atoms_dist_grad(1,3);

close([fig1 fig2 fig3 fig71_1 fig71_2 fig122_1 fig122_2 fig92_1 fig92_2 fig13_1 fig13_2]);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

%% Sidechain Approximation Validation Figures

% Generate Sidechain Validation Data

data_structs = sidechain_atom_extraction();
num_domains = length(data_structs);
input_data = process_sidechain_data(data_structs);

% Rotational Figures

amino_index = 10;
output_data = generate_distance_sidechain_data(data_structs, input_data, amino_index);

az = 0;
el = 0;

fig1_handle = generate_distance_sidechain_scatter_figure(input_data, output_data, amino_index, num_domains,az, el);

ax1 = fig1_handle.CurrentAxes;
set(ax1, 'DataAspectRatio', [1 1 1]);
view(ax1,az,el);
title(ax1,'');
zlabel(ax1,'Distance \r{A}');
ax1.XAxis.Visible = 'off';
ax1.YAxis.Visible = 'off';
exportgraphics(fig1_handle, 'figures/sc_rotation_1.pdf',"ContentType","vector");

az = 30;
el = 0;

fig2_handle = generate_distance_sidechain_scatter_figure(input_data, output_data, amino_index, num_domains,az, el);

ax2 = fig2_handle.CurrentAxes;
set(ax2, 'DataAspectRatio', [1 1 1]);
view(ax2,az,el);
% title(ax2,'30 Degrees Rotation');
title(ax2,'');
zlabel(ax2,'Distance \r{A}');
ax2.XAxis.Visible = 'off';
ax2.YAxis.Visible = 'off';
exportgraphics(fig2_handle, 'figures/sc_rotation_2.pdf',"ContentType","vector");

az = 60;
el = 0;

fig3_handle = generate_distance_sidechain_scatter_figure(input_data, output_data, amino_index, num_domains,az, el);

ax3 = fig3_handle.CurrentAxes;
set(ax3, 'DataAspectRatio', [1 1 1]);
view(ax3,az,el);
% title(ax3,'60 Degrees Rotation');
title(ax3,'');
zlabel(ax3,'Distance \r{A}');
ax3.XAxis.Visible = 'off';
ax3.YAxis.Visible = 'off';
exportgraphics(fig3_handle, 'figures/sc_rotation_3.pdf',"ContentType","vector");

az = 90;
el = 0;

fig4_handle = generate_distance_sidechain_scatter_figure(input_data, output_data, amino_index, num_domains,az, el);

ax4 = fig4_handle.CurrentAxes;
set(ax4, 'DataAspectRatio', [1 1 1]);
view(ax4,az,el);
% title(ax4,'90 Degrees Rotation');
title(ax4,'');
zlabel(ax4,'Distance \r{A}');
ax4.XAxis.Visible = 'off';
ax4.YAxis.Visible = 'off';
exportgraphics(fig4_handle, 'figures/sc_rotation_4.pdf',"ContentType","vector");

% Sidechain Figures

amino_index = 9;

output_data = generate_distance_sidechain_data(data_structs, input_data, amino_index);
fig4 = generate_distance_sidechain_scatter_figure(input_data, output_data, amino_index, num_domains, 0, el);
set(fig4.CurrentAxes, 'DataAspectRatio', [1 1 1]);
view(fig4.CurrentAxes,0,el);
title(fig4.CurrentAxes, '');
exportgraphics(fig4.CurrentAxes,'figures/sc_circle_scatter_1.pdf',"ContentType","vector");

amino_index = 10;

output_data = generate_distance_sidechain_data(data_structs, input_data, amino_index);
fig2 = generate_distance_sidechain_scatter_figure(input_data, output_data, amino_index, num_domains, 0, el);
set(fig2.CurrentAxes, 'DataAspectRatio', [1 1 1]);
view(fig2.CurrentAxes,0,el);
title(fig2.CurrentAxes, '');
exportgraphics(fig2.CurrentAxes,'figures/sc_circle_scatter_2.pdf',"ContentType","vector");

amino_index = 13;

output_data = generate_distance_sidechain_data(data_structs, input_data, amino_index);
fig3 = generate_distance_sidechain_scatter_figure(input_data, output_data, amino_index, num_domains, 0, el);
set(fig3.CurrentAxes, 'DataAspectRatio', [1 1 1]);
view(fig3.CurrentAxes,0,el);
title(fig3.CurrentAxes, '');
exportgraphics(fig2.CurrentAxes,'figures/sc_circle_scatter_3.pdf',"ContentType","vector");

amino_index = 18;

output_data = generate_distance_sidechain_data(data_structs, input_data, amino_index);
fig1 = generate_distance_sidechain_scatter_figure(input_data, output_data, amino_index, num_domains, 0, el);
ax1 = fig1.CurrentAxes;
set(fig1.CurrentAxes, 'DataAspectRatio', [1 1 1]);
view(fig1.CurrentAxes,0,el);
title(fig1.CurrentAxes, '');
exportgraphics(fig1.CurrentAxes,'figures/sc_circle_scatter_4.pdf',"ContentType","vector");

close([fig1 fig2 fig3 fig4 fig1_handle fig2_handle fig3_handle fig4_handle]);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

%% Static Force Sim Result Comparison With Digitised Data from Yao et al.

data = generate_force_extension_plot_data();
fig = figure();
plot(data.yao_stretch.x,data.yao_stretch.y, 'Color','Blue');

xlabel("Force, pN");
ylabel("Extension, nm");

hold on
plot(data.original.x,data.original.y, 'Color',[0.5 0.5 0.5]);
hold off

grid(fig.CurrentAxes, 'on');
grid(fig.CurrentAxes, 'minor');
legend(fig.CurrentAxes, 'Yao et al.', 'Sim Original', 'Location', 'southeast');
exportgraphics(fig.CurrentAxes, 'figures/yao_and_static_force_sim_results.pdf',"ContentType","vector");

close(fig);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

%% Helix to Helix Electrostatic Mapping Figures

% Figure One

sim = Simulation();
sim.Load_Amino_Acid_PDB_File(pdb_filepaths(6));
sim.Load_Helix(6,1,helix_indeces);
sim.Load_Helix(6,2,helix_indeces);
sim.Load_Helix(6,3,helix_indeces);
sim.Load_Helix(6,4,helix_indeces);
sim.Load_Helix(6,5,helix_indeces);
sim.Load_PDB_Parallel_Chain();
output = sim.Run_Interaction_Mapping(2,4);
sim.Run_Static_Force(2,4);

sim.Generate_Mapping_Figure(circshift(output,[length(output)/2,length(output)/2]));
xticklabels({'-180','-90','0','90','180'});
yticklabels({'-180','-90','0','90','180'});
view(-114,35);
xlabel('Helix Rotation, degrees');
ylabel('Helix Rotation, degrees');
zlabel('Interaction Force, N');
fig = gca();
exportgraphics(fig,'figures/mapping_figure.pdf',"ContentType","vector");

close(fig.Parent);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

% Figure Two

sim = Simulation();
sim.Load_Amino_Acid_PDB_File(pdb_filepaths(6));
sim.Load_Helix(6,1,helix_indeces);
sim.Load_Helix(6,2,helix_indeces);
sim.Load_Helix(6,3,helix_indeces);
sim.Load_Helix(6,4,helix_indeces);
sim.Load_Helix(6,5,helix_indeces);
sim.Load_PDB_Parallel_Chain();
output = sim.Run_Interaction_Mapping(2,4);
sim.Run_Static_Force(2,4);


sim.Generate_Mapping_Figure(circshift(output,[length(output)/2,length(output)/2]));
xticklabels({'-180','-90','0','90','180'});
yticklabels({'-180','-90','0','90','180'});
view(-49,35);
xlabel('Helix Rotation, degrees');
ylabel('Helix Rotation, degrees');
zlabel('Interaction Force, N');
fig2 = gca();
exportgraphics(fig2,'figures/mapping_figure2.pdf',"ContentType","vector");

close(fig2.Parent);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

%% Static Force Sim and Optimised Sim Results Comparison With Digitised Data from Yao et al.

data = generate_force_extension_plot_data();
fig = figure();

plot(data.yao_stretch.x,data.yao_stretch.y, 'Color','Blue');

xlabel("Force, pN");
ylabel("Extension, nm");

hold on
plot(data.collision.x,data.collision.y, 'Color', 'Red');
plot(data.original.x,data.original.y, 'Color',[0.5 0.5 0.5]);
hold off

grid(fig.CurrentAxes, 'on');
grid(fig.CurrentAxes, 'minor');
legend(fig.CurrentAxes, 'Yao et al.', 'Sim Optimised', 'Sim Original', 'Location', 'southeast');
exportgraphics(fig.CurrentAxes, 'figures/yao_and_optim_sim_results.pdf',"ContentType","vector");

close(fig);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

%% Hydrophobic Scatter Plots

% R4 Figure

domain_index = 4;
sim = Simulation.Load_Sim_Talin_Subdomain(domain_index);
sim.Calculate_Bundle_Axis();
data_raw = sim.Calculate_Bundle_Hydropathy_Data();

if length(data_raw) == 2
    fig_1 = sim.Make_Hydropathy_Sidechain_Scatter_Plot_3_Part();
    fig_1 = sim.Apply_Labels_Scatter_Plot(fig_1, data_raw(1));
    title(fig_1.CurrentAxes, sprintf("R%d, 3 Part - Hydropathy Data Scatter",domain_index));
    fig_2 = sim.Make_Hydropathy_Sidechain_Scatter_Plot_4_Part();
    fig_2 = sim.Apply_Labels_Scatter_Plot(fig_2, data_raw(2));
    title(fig_2.CurrentAxes, sprintf("R%d, 4 Part - Hydropathy Data Scatter",domain_index));
else
    fig_1 = sim.Make_Hydropathy_Sidechain_Scatter_Plot_4_Part();
    fig_1 = sim.Apply_Labels_Scatter_Plot(fig_1, data_raw(1));
    title(fig_1.CurrentAxes, sprintf("R%d, 4 Part - Hydropathy Data Scatter",domain_index));
end
title(fig_1.CurrentAxes,'');
exportgraphics(fig_1.CurrentAxes, 'figures/hydroscatter_1.pdf',"ContentType","vector");

close(fig_1);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

% R5 Figures

domain_index = 5;
sim = Simulation.Load_Sim_Talin_Subdomain(domain_index);
sim.Calculate_Bundle_Axis();
data_raw = sim.Calculate_Bundle_Hydropathy_Data();

if length(data_raw) == 2
    fig_1 = sim.Make_Hydropathy_Sidechain_Scatter_Plot_3_Part();
    fig_1 = sim.Apply_Labels_Scatter_Plot(fig_1, data_raw(1));
    title(fig_1.CurrentAxes, sprintf("R%d, 3 Part - Hydropathy Data Scatter",domain_index));
    fig_2 = sim.Make_Hydropathy_Sidechain_Scatter_Plot_4_Part();
    fig_2 = sim.Apply_Labels_Scatter_Plot(fig_2, data_raw(2));
    title(fig_2.CurrentAxes, sprintf("R%d, 4 Part - Hydropathy Data Scatter",domain_index));
else
    fig_1 = sim.Make_Hydropathy_Sidechain_Scatter_Plot_4_Part();
    fig_1 = sim.Apply_Labels_Scatter_Plot(fig_1, data_raw(1));
    title(fig_1.CurrentAxes, sprintf("R%d, 4 Part - Hydropathy Data Scatter",domain_index));
end
title(fig_1.CurrentAxes,'');
exportgraphics(fig_1.CurrentAxes, 'figures/hydroscatter_2a.pdf',"ContentType","vector");
title(fig_2.CurrentAxes,'');
exportgraphics(fig_2.CurrentAxes, 'figures/hydroscatter_2b.pdf',"ContentType","vector");

close([fig_1 fig_2]);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;

%%

hydrophobic_results = load('results\hydrophobic_results.mat');
hydrophobic_results = hydrophobic_results.hydrophobic_results;

bar_fig = figure();
bar_fig.CurrentAxes = axes();
bar(bar_fig.CurrentAxes,hydrophobic_results.bar_graph_n_close_hydrophobic_vs_hydrophilic_combined);
exportgraphics(bar_fig.CurrentAxes, 'figures/hydrophobic_meta_graph.pdf',"ContentType","vector");

close(bar_fig);
clearvars -except table_data helix_indeces pdb_filepaths residue_properties;