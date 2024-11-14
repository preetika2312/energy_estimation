function plot_combined_psd800Spatch_Oct7(fileNames, labels, saveDir)
% This function reads surface data from .datx files, performs planarization,
% plots the combined Power Spectral Density (PSD) for different objectives,
% calculates Hrms, fits a curve to log-transformed PSD data using lsqcurvefit,
% and saves all results in a specified directory.

% Ensure that both fileNames and labels are provided and have the correct format
if nargin < 3
    error('Not enough input arguments. fileNames, labels, and saveDir are required.');
end

if ~iscell(fileNames) || ~iscell(labels)
    error('fileNames and labels should be cell arrays.');
end

if length(fileNames) ~= length(labels)
    error('fileNames and labels must have the same length.');
end

% Check and create the save directory if it does not exist
if ~isfolder(saveDir)
    mkdir(saveDir);
    disp(['Created directory: ' saveDir]);
end

% Define the resolutions for each objective in micrometers
resolutions = [0.0029632, 0.00081237, 0.00016296] * 1e3;

numFiles = length(fileNames);

% Preallocate for individual PSD plots and Hrms values
individual_k_bins = cell(numFiles, 1);
individual_PSD_radial_avg = cell(numFiles, 1);
hrms_values = zeros(numFiles, 1);

% Define the constants beta, sigma, and G
beta = 1;
sigma = 20/1000;  % N/m
%G = 560000;  % Pa
% G = 186667;  % Pa
% G = 81667;  % Pa
% G = 16667;  % Pa
% G = 333;  % Pa
G=[560000 186667 81667 16667 333  ];

% Initialize the integral result for Ws
Ws_total = 0;

% Loop through each file to process data and generate plots
for i = 1:numFiles
    % Read surface data from file
    data = ReadDatx(fileNames{i});

    % Extract and process the surface matrix
    surface_matrix = data.Surface.matrix;
    surface_matrix = fillmissing(surface_matrix, 'linear');

    % Save the original surface_matrix to a new sheet in the Excel file
    surface_matrix_filename = fullfile(saveDir, ['PSD_radial_avg_' labels{i} '.xlsx']);
    writematrix(surface_matrix, surface_matrix_filename, 'Sheet', 'Surface_Matrix');

    % Perform planarization to remove curvature
    [X, Y] = meshgrid(1:size(surface_matrix, 2), 1:size(surface_matrix, 1));
    A = [X(:), Y(:), ones(numel(surface_matrix), 1)];
    coeff = A \ surface_matrix(:);
    surface_planarized = reshape(A * coeff, size(surface_matrix));
    surface_flat = surface_matrix - surface_planarized;

    % Calculate Hrms (root mean square height)
    hrms = sqrt(mean(surface_flat(:).^2, 'omitnan'));
    hrms_values(i) = hrms;  % Store Hrms value for each file


    % Define patch size
    patch_size = 10;

    % Calculate number of patches in each dimension
    num_patches_x = floor(size(surface_matrix, 1) / patch_size);
    num_patches_y = floor(size(surface_matrix, 2) / patch_size);

    % Initialize matrix to store h_rms values of each patch
    h_rms_patches = zeros(num_patches_x, num_patches_y);

    % Loop over each patch
    for px = 1:num_patches_x
        for py = 1:num_patches_y
            % Extract patch
            patch = surface_matrix((px-1)*patch_size + 1:px*patch_size, (py-1)*patch_size + 1:py*patch_size);

            % Perform planarization to remove curvature from each patch
            [X_patch, Y_patch] = meshgrid(1:size(patch, 2), 1:size(patch, 1));
            A_patch = [X_patch(:), Y_patch(:), ones(numel(patch), 1)];
            coeff_patch = A_patch \ patch(:);
            patch_planarized = reshape(A_patch * coeff_patch, size(patch));
            patch_flat = patch - patch_planarized;

            % Calculate h_rms for the planarized patch
            h_rms_patches(px, py) = sqrt(mean(patch_flat(:).^2, 'omitnan'));

            % % Save patch_flat of each patch to a new sheet in the same Excel file
            % patch_flat_sheet_name = sprintf('Patch_Flat_%d_%d', px, py);
            % writematrix(patch_flat, fullfile(saveDir, ['Hrms_patches_' labels{i} '.xlsx']), 'Sheet', patch_flat_sheet_name);
        end
    end

    % Save the h_rms of each patch and the average h_rms to a new Excel sheet
    h_rms_patches_filename = fullfile(saveDir, ['Hrms_patches_' labels{i} '.xlsx']);
    writematrix(h_rms_patches, h_rms_patches_filename, 'Sheet', 'Hrms_Patches');

    % Calculate and save the average h_rms of all patches
    average_h_rms = mean(h_rms_patches(:), 'omitnan');
    writematrix(average_h_rms, h_rms_patches_filename, 'Sheet', 'Average_Hrms', 'Range', 'A1');


    % Get the resolution for the current objective
    dx = resolutions(i);

    % Compute 2D Fourier Transform
    [n, m] = size(surface_flat);
    F = fft2(surface_flat) * dx * dx;
    F_shifted = fftshift(F);

    % Compute Power Spectral Density (PSD)
    PSD = abs(F_shifted).^2 / (n * m) / (dx * dx);

    % Compute the spatial frequency axes
    kx = (-n/2:n/2-1) * pi * 2 / (n * dx);
    ky = (-m/2:m/2-1) * pi * 2 / (m * dx);
    [kx_grid, ky_grid] = meshgrid(kx, ky);
    k = sqrt(kx_grid.^2 + ky_grid.^2);

    % Radial average of the PSD
    k_bins = linspace(0, max(kx(:)), 100);
    PSD_radial_avg = zeros(size(k_bins));

    for j = 1:length(k_bins)-1
        mask = (k >= k_bins(j)) & (k < k_bins(j+1));
        PSD_radial_avg(j) = mean(PSD(mask), 'omitnan');
    end

    % Store individual PSD data
    individual_k_bins{i} = k_bins(1:end-1);
    individual_PSD_radial_avg{i} = PSD_radial_avg(1:end-1);

    % Save individual PSD data and Hrms to Excel in the first sheet
    PSD_radial_avg_with_hrms = [individual_k_bins{i}', individual_PSD_radial_avg{i}', hrms * ones(length(individual_k_bins{i}), 1)];
    PSD_radial_avg_filename = fullfile(saveDir, ['PSD_radial_avg_' labels{i} '.xlsx']);
    writematrix(PSD_radial_avg_with_hrms, PSD_radial_avg_filename, 'Sheet', 'Data');

    % Save surface_flat matrix to a new sheet
    writematrix(surface_flat, PSD_radial_avg_filename, 'Sheet', 'Surface_Flat');

    % Plot and save surface data
    fig1 = figure;
    surf(X, Y, surface_flat);
    title(['Surface with Slope/Curvature Removed for ' labels{i}]);
    shading flat;
    colorbar;
    filename1_png = fullfile(saveDir, ['Surface_' labels{i} '.png']);
    filename1_fig = fullfile(saveDir, ['Surface_' labels{i} '.fig']);
    saveas(fig1, filename1_png);
    saveas(fig1, filename1_fig);
    close(fig1);

    % Plot and save 2D PSD
    fig2 = figure;
    imagesc(kx, ky, log(1 + PSD));
    colormap('jet');
    colorbar;
    axis equal;
    title(['2D PSD for ' labels{i}]);
    xlabel('kx (1/\mum)');
    ylabel('ky (1/\mum)');
    filename2_png = fullfile(saveDir, ['2D_PSD_' labels{i} '.png']);
    filename2_fig = fullfile(saveDir, ['2D_PSD_' labels{i} '.fig']);
    saveas(fig2, filename2_png);
    saveas(fig2, filename2_fig);
    close(fig2);

    % Plot and save radial PSD
    fig3 = figure;
    loglog(individual_k_bins{i}, individual_PSD_radial_avg{i}, 'LineWidth', 2, 'DisplayName', labels{i});
    title(['Radial Power Spectral Density for ' labels{i}]);
    xlabel('Spatial Frequency (\mum^{-1})');
    ylabel('Power');
    legend show;
    grid on;
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    filename3_png = fullfile(saveDir, ['Radial_PSD_' labels{i} '.png']);
    filename3_fig = fullfile(saveDir, ['Radial_PSD_' labels{i} '.fig']);
    saveas(fig3, filename3_png);
    saveas(fig3, filename3_fig);
    close(fig3);
end

% Combine PSD data for all files
x_together = [];
y_together = [];
for i = 1:numFiles
    x_together = [x_together; individual_k_bins{i}];
    y_together = [y_together; individual_PSD_radial_avg{i}];
end

% Remove non-positive values from the combined data
valid_indices = (x_together > 0) & (y_together > 0);
x_together = x_together(valid_indices)*(10^6);
y_together = y_together(valid_indices)*(10^-24);

% Transform data to log-log scale
log_x_together = log(x_together);
log_y_together = log(y_together);
log_z_together = log_y_together./log_x_together;

% Define the linear model for log-log data
linearModel = @(params, x) (params(1) + params(2) * x)/x;  % Model: log(y)/log(x) = (log(a) + b * log(x))/log(x)

% Set initial guesses for the parameters [log(a), b]
initialGuess = [0, -1]; % Adjust these values as needed

% Perform the fitting using lsqcurvefit
options = optimset('Display', 'off');  % Suppress output display
[fittedParams, resnorm] = lsqcurvefit(linearModel, initialGuess, log_x_together(8,:), log_z_together(8,:), [], [], options);

% Extract fitted parameters
log_a_fitted = fittedParams(1);
b_fitted = fittedParams(2);
a_fitted = exp(log_a_fitted);

% Generate fitted values in the original scale
y_fit = (a_fitted * x_together.^b_fitted);

% Prepare data to be saved in the new sheet
fit_data = table(x_together, y_together, y_fit, 'VariableNames', {'SpatialFrequency', 'OriginalPSD', 'FittedPSD'});
fit_params = table(log_a_fitted, a_fitted, b_fitted, 'VariableNames', {'log_a_fitted','a_fitted', 'b_fitted'});

% Write the fitted data and parameters to a new sheet in the Excel file
writetable(fit_data, PSD_radial_avg_filename, 'Sheet', 'Fit Results');
writetable(fit_params, PSD_radial_avg_filename, 'Sheet', 'Fit Results', 'Range', 'E1');

% Plot combined data and fitted curve
fig5 = figure;
loglog(x_together, y_together, 'o', 'DisplayName', 'Combined PSD Data');
hold on;
box on;
loglog(x_together, y_fit, '-', 'LineWidth', 2, 'DisplayName', sprintf('Fitted Curve: y = %.2f * x^{%.2f}', a_fitted, b_fitted));
title('Combined Radial PSD Data with Fitted Curve');
%xlabel('Spatial Frequency (\mum^{-1})');
%ylabel('Power');
box on;
xlabel('q (m^{-1})', 'FontSize',14)
ylabel('C_i_s_o (m^4)', 'FontSize',14)
legend show;
grid on;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
hold off;
combinedFilename_fit_png = fullfile(saveDir, 'Combined_Radial_fit_PSD.png');
combinedFilename_fit_fig = fullfile(saveDir, 'Combined_Radial_fit_PSD.fig');
saveas(fig5, combinedFilename_fit_png);
saveas(fig5, combinedFilename_fit_fig);
close(fig5);

% Plot and save combined radial PSDs for all files
fig4 = figure;
hold on;
colors = {'r', 'g', 'b'};
for i = 1:numFiles
    loglog(individual_k_bins{i}, individual_PSD_radial_avg{i}, 'LineWidth', 2, 'Color', colors{i}, 'DisplayName', labels{i});
end
xlabel('q (\mum^{-1})')
ylabel('C_i_s_o (\mum^4)')
legend show;
grid on;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
hold off;
combinedFilename_png = fullfile(saveDir, 'Combined_Radial_PSD.png');
combinedFilename_fig = fullfile(saveDir, 'Combined_Radial_PSD.fig');
saveas(fig4, combinedFilename_png);
saveas(fig4, combinedFilename_fig);
close(fig4);

% Plot and save combined radial PSDs for all files in SI units
fig6 = figure;
hold on;
for i = 1:numFiles
    loglog(individual_k_bins{i}*(10^6), individual_PSD_radial_avg{i}*(10^-24), 'LineWidth', 2, 'Color', colors{i}, 'DisplayName', labels{i});
end
xlabel('q (m^{-1})', 'FontSize',14)
ylabel('C_i_s_o (m^4)', 'FontSize',14)
legend show;
box on;
grid on;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
hold off;
combinedFilename_SI_png = fullfile(saveDir, 'Combined_Radial_PSD_SI.png');
combinedFilename_SI_fig = fullfile(saveDir, 'Combined_Radial_PSD_SI.fig');
saveas(fig6, combinedFilename_SI_png);
saveas(fig6, combinedFilename_SI_fig);
close(fig6);

% Confirmation message
disp(['All figures, PSD data, and Hrms values have been saved in the directory: ' saveDir]);
% Initialize the integral result
integral_result = 0;

% Loop through each file's data
for i = 1:numFiles
    % Get the values of q (spatial frequency) and C(q) (PSD) for this file
    q_values = individual_k_bins{i}*(10^6);  % q is the spatial frequency
    C_values = individual_PSD_radial_avg{i}*(10^-24);  % C(q) is the PSD radial average
end

% Define the integrand q^2 * C(q)
integrand = @(q, Cq) q.^2 .* Cq/4/pi/pi;

% Numerically integrate q^2 * C(q) over the range of q using the trapezoidal method
UEL_integrate = trapz(q_values, integrand(q_values, C_values));

% Add the contribution from this file to the total integral result
integral_result = integral_result + UEL_integrate;

% Define the integrand for Ws
integrand_Ws = @(q, Cq, G) (q.^4 .* Cq/4/pi/pi) ./ (1 + (beta * sigma / (2 * G)) * q).^2;
for j=1:5

    % Numerically integrate the integrand over the range of q using the trapezoidal method
    Ws_integrate(j) = trapz(q_values, integrand_Ws(q_values, C_values, G(j)));

    % Add the contribution from this file to the total Ws result
    %Ws_total = Ws_total + Ws_integrate;


    % Display the final result of the integral
    disp(['Total UEL_integrate = ', num2str(integral_result)]);

    % Save the total UEL value in an Excel file
    output_filename = fullfile(saveDir, 'Estimations.xlsx');  % Define the output file path
    writematrix(integral_result, output_filename, 'Sheet', 'UEL_integrate', 'Range', 'A1');

    % Compute the final value of Ws
    Ws(j) = sigma * pi * (beta * sigma / (2 * G(j))) * Ws_integrate(j);
    %
    % Display the final result of Ws
    disp(['Total Ws = ', num2str(Ws)]);


    % Save the total Ws value in an Excel file
    output_filename = fullfile(saveDir, 'Estimations.xlsx');  % Define the output file path
    writematrix(Ws(j), output_filename, 'Sheet', 'Ws_Result', 'Range', 'A1');
end


end


% fileNames = {'/Volumes/Extreme SSD/Project1/Zegage/2021_2023/October2023/Prep2_FrostdGlass_Oct3/2p75x.datx', '/Volumes/Extreme SSD/Project1/Zegage/2021_2023/October2023/Prep2_FrostdGlass_Oct3/10x.datx', '/Volumes/Extreme SSD/Project1/Zegage/2021_2023/October2023/Prep2_FrostdGlass_Oct3/50x.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Volumes/Extreme SSD/Project1/Zegage/2021_2023/October2023/Prep2_FrostdGlass_Oct3/plotspsd800Spatch_integrand_Oct7_update_oct18_try2'
% plot_combined_psd800Spatch_refit_Oct18(fileNames, labels, saveDir);

% fileNames = {'/Volumes/Extreme SSD/Project1/Zegage/2024/FrostedGlassSlide_Jul_2024/2p75x_10scans.datx', '/Volumes/Extreme SSD/Project1/Zegage/2024/FrostedGlassSlide_Jul_2024/10x_10scans.datx', '/Volumes/Extreme SSD/Project1/Zegage/2024/FrostedGlassSlide_Jul_2024/50x_10scans.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Volumes/Extreme SSD/Project1/Zegage/2024/FrostedGlassSlide_Jul_2024/plotspsd800Spatch'
% plot_combined_psd800Spatch(fileNames, labels, saveDir);

% fileNames = {'/Volumes/Extreme SSD/Project1/Zegage/2024/FrostedGlassSlide_May30_2023/2p75x.datx', '/Volumes/Extreme SSD/Project1/Zegage/2024/FrostedGlassSlide_May30_2023/10x.datx', '/Volumes/Extreme SSD/Project1/Zegage/2024/FrostedGlassSlide_May30_2023/50x.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Volumes/Extreme SSD/Project1/Zegage/2024/FrostedGlassSlide_May30_2023/plotspsd800Spatch'
% plot_combined_psd800Spatch(fileNames, labels, saveDir);

% fileNames = {'/Volumes/Extreme SSD/Project1/Zegage/2024/FrostedGlassSlide_Sep13/2p75x.datx', '/Volumes/Extreme SSD/Project1/Zegage/2024/FrostedGlassSlide_Sep13/10x.datx', '/Volumes/Extreme SSD/Project1/Zegage/2024/FrostedGlassSlide_Sep13/50x.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Volumes/Extreme SSD/Project1/Zegage/2024/FrostedGlassSlide_Sep13/plotspsd800Spatch'
% plot_combined_psd800Spatch_Oct4_2024_integrand(fileNames, labels, saveDir);


% fileNames = {'/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/10_1_March2024_rough/Rough_2p5x.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/10_1_March2024_rough/Rough_10x.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/10_1_March2024_rough/Rough_50x.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/10_1_March2024_rough/plotspsd800Spatch_integrand_Oct7_update_oct14_3'
% plot_combined_psd800Spatch_Oct7(fileNames, labels, saveDir);


% fileNames = {'/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/10_1_March2024_rough/Smooth_Aug26/2p75x_spot2.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/10_1_March2024_rough/Smooth_Aug26/10x_spot2.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/10_1_March2024_rough/Smooth_Aug26/50x_spot2.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/10_1_March2024_rough/Smooth_Aug26/plotspsd800S_smooth_patch'
% plot_combined_psd800Spatch(fileNames, labels, saveDir);

% fileNames = {'/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/10_1_March2024_rough/Smooth_2p5x.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/10_1_March2024_rough/Smooth_10x.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/10_1_March2024_rough/50x.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/10_1_March2024_rough/plotspsd800S_smooth_patch'
% plot_combined_psd800Spatch(fileNames, labels, saveDir);

% fileNames = {'/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/20_1_PDMS_March2024_rough/2p75x_10scans.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/20_1_PDMS_March2024_rough/10x_spot2.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/20_1_PDMS_March2024_rough/50x.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/20_1_PDMS_March2024_rough/plotspsd800Spatch_integrand_Oct7_update_oct14_2'
% plot_combined_psd800Spatch_Oct7(fileNames, labels, saveDir);

% fileNames = {'/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/20_1_PDMS_March2024_rough/Smooth_Aug28/2p75x.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/20_1_PDMS_March2024_rough/Smooth_Aug28/10x.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/20_1_PDMS_March2024_rough/Smooth_Aug28/50x.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/20_1_PDMS_March2024_rough/Smooth_Aug28/plotspsd800Spatch'
% plot_combined_psd800Spatch(fileNames, labels, saveDir);

%
% fileNames = {'/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/30_1_PDMS_rough/2p75x_spot2_jul30.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/30_1_PDMS_rough/10x_spot3.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/30_1_PDMS_rough/50x_jul30_10scans.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/30_1_PDMS_rough/plotspsd800Spatch_integrand_Oct7_update_oct14_2'
% plot_combined_psd800Spatch_Oct7(fileNames, labels, saveDir);

% fileNames = {'/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/30_1_PDMS_rough/smooth/2p75x.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/30_1_PDMS_rough/smooth/10x.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/30_1_PDMS_rough/smooth/50x.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/30_1_PDMS_rough/smooth/plotspsd800Spatch'
% plot_combined_psd800Spatch(fileNames, labels, saveDir);

% fileNames = {'/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/40_1_PDMS_March2024/rough_2p75x_spot2.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/40_1_PDMS_March2024/rough_10x_spot2.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/40_1_PDMS_March2024/rough_50x_spot2.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/40_1_PDMS_March2024/plotspsd800Spatch_integrand_Oct7_updatedoct14_2'
% plot_combined_psd800Spatch_Oct7(fileNames, labels, saveDir);

% fileNames = {'/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/40_1_PDMS_March2024/smooth/2p75x_spot2.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/40_1_PDMS_March2024/smooth/10x_spot2.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/40_1_PDMS_March2024/smooth/50x_spot2.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/40_1_PDMS_March2024/smooth/plotspsd800Spatch'
% plot_combined_psd800Spatch(fileNames, labels, saveDir);

% fileNames = {'/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/GelestPDMS_9_1/Jul26sample_rough/2p75x_20Aug24.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/GelestPDMS_9_1/Jul26sample_rough/10x_spot1.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/GelestPDMS_9_1/Jul26sample_rough/50x_spot4.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/GelestPDMS_9_1/Jul26sample_rough/plotspsd800Spatch_integrand_Oct7_updatedoct14_2'
% plot_combined_psd800Spatch_Oct7(fileNames, labels, saveDir);

% fileNames = {'/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/GelestPDMS_9_1/August25_smooth/2p75x_10scans_spot2.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/GelestPDMS_9_1/August25_smooth/10x_10scans_spot2.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/GelestPDMS_9_1/August25_smooth/50x_10scans_spot2.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/GelestPDMS_9_1/August25_smooth/plotspsd800S_smooth_patch'
% plot_combined_psd800Spatch(fileNames, labels, saveDir);

% fileNames = {'/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/GlassSlide/Smooth/2p75x.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/GlassSlide/Smooth/10x.datx', '/Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/GlassSlide/Smooth/50x.datx'};
% labels = {'2.75x', '10x', '50x'};
% saveDir = '//Users/preetikakarnal/Documents/Work/Projects/Effect_of_surface_stress/Data/Zegage/PreetikaZegage/GlassSlide/Smooth/plotspsd800S_smooth_patch'
% plot_combined_psd800Spatch(fileNames, labels, saveDir);
