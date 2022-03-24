function [fig1,fig2] = bp_atoms_dist(pdb_file_index, helix_index, nbins, atom)
    silent_startup;
    
    if(nargin >= 3)
        atom = "CA";
    elseif (nargin < 3)
        nbins = 10;
        atom = "CA";
    end
    
    chain = Amino_Acid_Chain(pdb_filepaths(pdb_file_index));
    indeces = [helix_indeces(pdb_file_index,helix_index,1);helix_indeces(pdb_file_index,helix_index,2)];
    chain.Select_Data(indeces);
    bb_atoms = Amino_Acid_Chain.Get_Backbone_Atom_Positions(pdb_filepaths(pdb_file_index),indeces(1),indeces(2));
    backbone_positions = chain.backbone_positions(indeces(1):indeces(2),:);
    helix = chain.Make_Helix_Object_V2();

    bp_num = length(backbone_positions);
    distances = zeros(bp_num,1);

    
    atom = convertCharsToStrings(atom);
    if(atom == "CA")
         atom_positions = backbone_positions;
    elseif(atom == "C")
        atom_positions = bb_atoms.c_positions;
    elseif(atom == "O")
        atom_positions = bb_atoms.o_positions;
    elseif(atom == "N")
        atom_positions = bb_atoms.n_positions;
    end
    
    for i = 1:bp_num
        position = Point_Line_Intersection(helix.position_One,helix.position_Two,atom_positions(i,:));
        distances(i) = norm(backbone_positions(i,:) - position);
    end
    
    distance_count = length(distances);
    
    mean_distance = mean(distances);
    y_mean_line = ones(1,distance_count)*mean_distance;
    
    x = 1:distance_count;
    fit = polyfit(x,distances',1);
    fit_line = polyval(fit,x);
    
    fig1 = figure();
    axes(fig1);
    plot(fig1.CurrentAxes,x,distances,'--b','LineWidth',1.5);
    hold(fig1.CurrentAxes,'on');
    plot(fig1.CurrentAxes,x,fit_line,'r','LineWidth',2);
    plot(fig1.CurrentAxes,x,y_mean_line,'g','LineWidth',2);
    plot(fig1.CurrentAxes,x,distances,'xb','LineWidth',2,'MarkerSize', 12);
    legend(fig1.CurrentAxes,'Distance Spread','Regression Line','Mean Value')
    title(fig1.CurrentAxes,strcat("R",num2str(pdb_file_index),", Helix ", num2str(helix_index),", ",convertCharsToStrings(atom), "'s"));
    xlabel(fig1.CurrentAxes,"Backbone Index");
    ylabel(fig1.CurrentAxes,"Distance from Helix Axis");
    hold(fig1.CurrentAxes,'off');
%     disp(helix.radius);
    fig2 = figure();
    axes(fig2);
    histfit(fig2.CurrentAxes,distances,nbins);
    title(fig2.CurrentAxes,strcat("R",num2str(pdb_file_index),", Helix ", num2str(helix_index)));    
    xlabel(fig2.CurrentAxes,"Distance, \r{A}", 'Interpreter', 'latex');
    ylabel(fig2.CurrentAxes,"Instances per Bin");
end

