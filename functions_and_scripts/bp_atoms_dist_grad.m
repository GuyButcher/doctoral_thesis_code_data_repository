function gradient = bp_atoms_dist_grad(pdb_file_index, helix_index)
%BP_ATOMS_DIST_OUTPUT Summary of this function goes here
%   Detailed explanation goes here
    silent_startup;

    chain = Amino_Acid_Chain(pdb_filepaths(pdb_file_index));
    indeces = [helix_indeces(pdb_file_index,helix_index,1);helix_indeces(pdb_file_index,helix_index,2)];
    chain.Select_Data(indeces);
    bb_atoms = Amino_Acid_Chain.Get_Backbone_Atom_Positions(pdb_filepaths(pdb_file_index),indeces(1),indeces(2));
    backbone_positions = chain.backbone_positions(indeces(1):indeces(2),:);
    helix = chain.Make_Helix_Object_V2();

    bp_num = length(backbone_positions);
    distances = zeros(bp_num,1);
    
    for i = 1:bp_num
        position = Point_Line_Intersection(helix.position_One,helix.position_Two,backbone_positions(i,:));
        distances(i) = norm(backbone_positions(i,:) - position);
    end
    
    distance_count = length(distances);
    
    mean_distance = mean(distances);
    y_mean_line = ones(1,distance_count)*mean_distance;
    
    x = 1:distance_count;
    fit = polyfit(x,distances',1);
    fit_line = polyval(fit,x);
    
    gradient = atan2d(fit_line(end) - fit_line(1), distance_count);
end

