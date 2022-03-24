function latex_string = tableLatex(input_table, with_table_environment, with_resize, table_name)

    export_as_file = false;
    with_resize = false;

    if (nargin < 2)
        with_table_environment = false;
        table_name = '';
    elseif (nargin < 3)
        table_name = '';
    elseif (nargin < 4)
        table_name = '';
    elseif (nargin == 4)
        export_as_file = true;
    end

    num_cols = width(input_table);
    num_rows = height(input_table);
    table_spec = "";
    for i = 1:num_cols
        table_spec = table_spec + "l";
    end
    table_spec = "@{}" + table_spec + "@{}";

    latex_string = "\begin{tabular}{" + table_spec + "}" + newline;
    latex_string = latex_string + "\toprule" + newline;
    col_names = "";
    for i = 1:num_cols
        if (i == num_cols)
            col_names = col_names + input_table.Properties.VariableNames{i} + " \\ \midrule";
        else
            col_names = col_names + input_table.Properties.VariableNames{i} + " & ";
        end
    end
    col_names = strrep(col_names, "_", "\_");
    latex_string = latex_string + col_names + newline;

    for i = 1:num_rows
        current_row = "";
        row_data = input_table(i,:).Variables;
        for j = 1:num_cols
            if(j == num_cols)
                current_row = current_row + row_data(j) + " \\" + newline;
            else
                current_row = current_row + row_data(j) + " & ";
            end
        end
        current_row = strrep(current_row, "_", "\_");
        latex_string = latex_string + current_row;
    end
    latex_string = latex_string + "\bottomrule" + newline + "\end{tabular}";
    
%     if (with_table_environment)
%         latex_string = "\begin{table}[h!]" + newline + "\centering" + newline + "\resizebox{\textwidth}{!}{%" + newline + ...
%             latex_string + "%" + newline + "}" + newline + "\caption{}" + newline + ...
%             "\label{tab:}" + newline + "\end{table}";
%     end
    
    if (with_table_environment)
        latex_string_temp = "\begin{table}[h!]" + newline + "\centering" + newline;
        if (with_resize)
           latex_string_temp = latex_string_temp + "\resizebox{\textwidth}{!}{%" + newline;
        end
        latex_string_temp = latex_string_temp + latex_string;
        if (with_resize)
            latex_string_temp = latex_string_temp + "%" + newline + "}";
        end
        latex_string_temp = latex_string_temp + newline + "\caption{}" + newline + ...
            "\label{tab:}" + newline + "\end{table}";
        latex_string = latex_string_temp;
    end

    clipboard('copy',latex_string);
        
    if export_as_file
        fileID = fopen(strcat(table_name,'.tex'),'w');
        fprintf(fileID,'%s', latex_string);
        fclose(fileID);  
    end
    
    fprintf('LaTex Table Copied to Clipboard.\n\n');
end

