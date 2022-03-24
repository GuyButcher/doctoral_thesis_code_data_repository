experimental_results.force = [5,8,11,12,12.5,13,14.5,5,16,17,17,20,24];
experimental_results.x_end = [200,240,290,330,370,400,440,440,530,580,630,690,750];

x_start = 154;
x_finish = 1092;
x_diff = x_finish - x_start;

x_axis_pixel_coords = [309,489,504,517,548,591,615,0,637,686,720,809,902];
x_axis_force_values = round((((x_axis_pixel_coords - x_start) / x_diff) * 30),3,'significant');

y_start = 757;
y_finish = 84;
y_diff = y_start - y_finish;

y_axis_pixel_coords = [668,633,596,563,530,491,452,0,421,327,285,223,172];
y_axis_start_distance_values = round((((1 - ((y_axis_pixel_coords - y_finish)/y_diff)) * 700) + 100),3,'significant');

experimental_results.force = x_axis_force_values;
experimental_results.extension = y_axis_start_distance_values;