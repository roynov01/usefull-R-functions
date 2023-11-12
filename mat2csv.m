%output_file_name = 'X:\Common\Roys plotter\Shani_liver_mouse_atlas_shalev_analysis\Smillie_colon.csv';
output_file_name = 'X:\roy\resources\pythonGUI\datasets\Smillie_colon.csv';

colnames = cell_types;
rownames = gene_name;
data = mean_mat_norm;

matrixTable = array2table(data, 'RowNames', rownames, 'VariableNames', colnames);
writetable(matrixTable, output_file_name, 'Delimiter', ',', 'WriteRowNames', true);

