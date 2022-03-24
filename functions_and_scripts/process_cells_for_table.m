function output = process_cells_for_table(cell_array)
    num_cells = length(cell_array);
    string_array = strings(num_cells,1);
    for i = 1:num_cells
        temp_data = round(cell_array{i},3,'significant');
        temp_string = "{[}";
        for j = 1:length(temp_data)
            if (j == length(temp_data))
                temp_string = temp_string +  num2str(temp_data(j) + "{]}");
            else
                temp_string = temp_string +  num2str(temp_data(j) + ", ");
            end
        end
        string_array(i) = temp_string;
    end
    output = string_array;
end