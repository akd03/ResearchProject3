%% PLOT RADIAL BASIS FUNCTIONS
clear; clc; close all;

% --- 1. DEFINE DOMAIN ---
% Normalized radius from 0 to 1
r = linspace(0, 1, 500);

% --- 2. CALCULATE BASIS FUNCTIONS ---
% Norm (Linear / Distance)
phi_norm = r;

% Wendland C0
% Mathematically: (1-r)^2 for r <= 1
phi_WC0 = (1 - r).^2;

% Wendland C2
% Mathematically: (1-r)^4 * (4r + 1) for r <= 1
phi_WC2 = ((1 - r).^4) .* (4 * r + 1);

% Inverse Multiquadric (IMQ)
% Mathematically: 1 / sqrt(1 + (e*r)^2). Global support, decaying.
epsilon_imq = 3;
phi_IMQ = 1 ./ sqrt(1 + (epsilon_imq * r).^2);

% Gaussian
% Mathematically: exp(-(e*r)^2). Global support, decaying.
epsilon_gauss = 3; 
phi_gauss = exp(-(epsilon_gauss * r).^2);

% --- 3. PLOT FUNCTIONS ---
fig_rbf = figure('Name', 'Radial Basis Functions', 'Position', [100, 100, 800, 600]);
hold on; grid on;

%plot(r, phi_norm, '-', 'LineWidth', 2, 'DisplayName', 'Norm (r)');
plot(r, phi_WC0, '-', 'LineWidth', 1.5, 'DisplayName', 'Wendland C^0');
plot(r, phi_WC2, '-', 'LineWidth', 1.5, 'DisplayName', 'Wendland C^2');
plot(r, phi_IMQ, '-', 'LineWidth', 1.5, 'DisplayName', sprintf('Inverse Multiquadric (\\epsilon=%d)', epsilon_imq));
plot(r, phi_gauss, '-', 'LineWidth', 1.5, 'DisplayName', sprintf('Gaussian (\\epsilon=%d)', epsilon_gauss));

% Formatting
xlabel('Normalized Radius (r)');
ylabel('\phi(r)');
title('Comparison of Radial Basis Functions');
legend('Location', 'best', 'FontSize', 11);
xlim([0, 1]);
hold off;

% --- 4. SAVE FIGURE ---
% Ensure the target directory exists
output_dir = '04_figures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Generate dynamic filename starting with F1_
date_str = datestr(now, 'yyyymmdd');
save_filename = sprintf('F1_BasisFunctions_%s', date_str);
save_filepath = fullfile(output_dir, [save_filename, '.fig']);

% Save as .fig for future editing
savefig(fig_rbf, save_filepath);

% Optional: Also save as a PDF for easy insertion into reports
%exportgraphics(fig_rbf, fullfile(output_dir, [save_filename, '.pdf']), 'ContentType', 'vector');

fprintf('Figure successfully saved to: %s\n', save_filepath);