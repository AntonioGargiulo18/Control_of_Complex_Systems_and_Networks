clc; clear all; close all; 
currentDir = pwd;
savePath = fullfile(currentDir, 'figures');


italy = shaperead('gadm41_ITA_0.shp','UseGeoCoords',true);
italy_reg =  shaperead('gadm41_ITA_1.shp','UseGeoCoords',true);

%% 1. Preliminary data

region_number = 20; 
switch region_number 
    case 3
        % Existing code for 3 regions (optional adjustments similar to 20)
    case 20
        % All Italian regions
        regions = { ...
            'Sicily', 'Calabria', 'Sardegna', 'Basilicata', 'Campania', ...
            'Apulia', 'Molise', 'Lazio', 'Abruzzo', 'Umbria', 'Marche', ...
            'Toscana', 'Liguria', 'Emilia-Romagna', 'Piemonte', 'Veneto', ...
            'Lombardia', 'Friuli-Venezia Giulia', 'Valle d Aosta', 'Trentino-Alto Adige'};

        % Extract region centroids from shapefile
        latitudes = zeros(1, 20);
        longitudes = zeros(1, 20);
        shapeRegions = {italy_reg.NAME_1};

        for i = 1:length(regions)
            idx = find(strcmpi(shapeRegions, regions{i}));
            if isempty(idx)
                error('Region %s not found in shapefile.', regions{i});
            end
            % Get polygon coordinates (handle multipart regions)
            lat = italy_reg(idx).Lat;
            lon = italy_reg(idx).Lon;
            if iscell(lat)
                lat = vertcat(lat{:});
                lon = vertcat(lon{:});
            end
            % Approximate centroid by mean
            latitudes(i) = mean(lat);
            longitudes(i) = mean(lon);
        end

        % Compute distance matrix
        D = get_distance_matrix(regions, latitudes, longitudes);
end


% Node color gradient from South to North
colorSouth = [1, 0, 0];  % Red
colorNorth = [0, 0, 1];  % Blue
colors = interp1([min(latitudes), max(latitudes)], [colorSouth; colorNorth], latitudes);

% Create graph
G = digraph(D, regions);

% Plot settings
figure('Position', [100 100 1200 800]);
latlim = [35 48]; 
lonlim = [6 19];
axesm('mercator', 'MapLatLimit', latlim, 'MapLonLimit', lonlim);
geoshow(italy, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'k');
hold on;

% Project node coordinates to map projection
[x, y] = mfwdtran(latitudes, longitudes);

% Plot graph
h = plot(G, 'XData', x, 'YData', y, ...
    'NodeLabel', regions, ...
    'EdgeLabel', [], ...
    'MarkerSize', 8, ...
    'NodeColor', colors, ...
    'LineWidth', 0.5, ...
    'EdgeColor', [0.4 0.4 0.4]);

% Enhance plot
title('Graph of Italian Regions');
h.EdgeAlpha = 0.3;
axis equal off;
grid off;
colormap(jet);
colorbar('southoutside', 'Ticks', linspace(min(latitudes), max(latitudes), 5), ...
         'TickLabels', arrayfun(@num2str, linspace(min(latitudes), max(latitudes), 5), 'UniformOutput', false));
%% 2. Non-linear model and constitutive equations
% Initial conditions for the number of population
n_0 = [4779371; 5575025; 10035481];  %year 2024 from ISTAT website

alpha_initial = [0.5; 0.5; 0.5];

T = 0:20; 
observed_rates  = [-2.8; -3.3; 1.3];  %year 2024 from ISTAT website  
reference_rates = [-2.8; -1.5; 1.3];  %reference alpha used for the first control attempt with three regions

lb = zeros(size(alpha_initial)) + eps; % Avoid alpha = 0
ub = ones(size(alpha_initial));
options = optimoptions('lsqnonlin', 'Display', 'iter', 'FunctionTolerance', 1e-12, 'StepTolerance',0.001);
alpha_opt = lsqnonlin(@(alpha_initial) objective_function(alpha_initial, n_0, D, T, observed_rates), ...
                      alpha_initial, lb, ub, options);

alpha_ref = lsqnonlin(@(alpha_initial) objective_function(alpha_initial, n_0, D, T, reference_rates), ...
                      alpha_initial, lb, ub, options);
%%
%alpha_opt = [0.3732; 0.3065; 0.6669];

%alpha_opt = [0.3734; 0.3059; 0.6674]; %Open - loop for observed rates

[n, n_m, r, J, outflux, influx, total_flux, n_net] = solve_continuity_equation (n_0, alpha_opt, D, T);
%[n_pin, n_m_pin, r_pin, J_pin, outflux_pin, influx_pin, total_flux_pin, n_net_pin] = solve_continuity_equation (n_0, alpha_ref, D, T);

[n_pin, n_m_pin, r_pin, J_pin, outflux_pin, influx_pin, total_flux_pin, n_net_pin, pinning_input] = solve_continuity_equation_pinning (n_0, alpha_opt, D, T);



%% PLOTS
save_flag =1; 
pin_flag = 1; 

main_name = '3_regions_pinning_camp_'; 
figure;
name_fig = strcat(main_name,'Population_Change_Over_Time_Relative_to_Initial_Value'); 
fileName = fullfile(savePath, name_fig);
grid on; 
hold on;

legend_entries = gobjects(length(regions),1); % Preallocate for legend handles

for i = 1:length(regions)
    % Plot main (solid) line and save handle for legend
    h_main = plot(T, ((n(i,:) - n_0(i)) ./ n_0(i)) * 100, ...
                  'LineWidth', 2, 'Color', colors(i,:));
    legend_entries(i) = h_main; % Save handle only for solid lines

    % Optional: plot pin (dashed) line without legend
    if pin_flag
        plot(T, ((n_pin(i,:) - n_0(i)) ./ n_0(i)) * 100, ...
             'LineStyle', '--', 'LineWidth', 2, 'Color', colors(i,:), ...
             'HandleVisibility', 'off');
    end
end

legend(legend_entries, regions, 'Location', 'best', 'FontSize', 10);
title('Population Change Over Time Relative to Initial Value [%]'); 
xlabel('k', 'FontSize', 10, 'FontWeight', 'bold');
set(gca, 'FontSize', 10); 
set(gcf, 'PaperPositionMode', 'auto'); 

if save_flag
    print(gcf, fileName, '-dsvg', '-r300');
    print(gcf, fileName, '-depsc', '-r0');
end

% NET MIGRATION FLUX 
figure; 
name_fig = strcat(main_name,'net_migration_flux_over_time'); 
fileName = fullfile(savePath, name_fig);
grid on; 
hold on;

legend_entries = gobjects(length(regions), 1);  % Preallocate legend handles

for i = 1:length(regions)
    % Plot solid line and store handle for legend
    h_main = plot(T, total_flux(i,:), 'LineWidth', 2, 'Color', colors(i,:));
    legend_entries(i) = h_main;

    % Plot dashed line without adding to legend
    if pin_flag
        plot(T, total_flux_pin(i,:), 'LineStyle', '--', 'LineWidth', 2, ...
             'Color', colors(i,:), 'HandleVisibility', 'off');
    end
end

legend(legend_entries, regions, 'Location', 'best', 'FontSize', 10);
title('Net migration flux over time'); 
xlabel('k', 'FontSize', 10, 'FontWeight', 'bold');
set(gca, 'FontSize', 10); 
set(gcf, 'PaperPositionMode', 'auto'); 

if save_flag
    print(gcf, fileName, '-dsvg', '-r300');
    print(gcf, fileName, '-depsc', '-r0');
end


% Net population

figure;
name_fig = strcat(main_name,'net_population_bar'); 
fileName = fullfile(savePath, name_fig);
bar(n_net,  'LineWidth', 2); 
title('Net Population Variation'); 
xticks(1:length(regions)); 
xticklabels(regions); 
xtickangle(45); 
hold on;
positive_diff = n_net > 0;
negative_diff = n_net < 0;
bar(find(positive_diff), n_net(positive_diff), 'FaceColor', 'g'); % Green for growth
bar(find(negative_diff), n_net(negative_diff), 'FaceColor', 'r'); % Red for decline
hold off;
if pin_flag
hold on;
    positive_diff = n_net_pin > 0;
    negative_diff = n_net_pin < 0;
    b1=bar(find(positive_diff), n_net_pin(positive_diff), 'FaceColor', 'g'); % Green for growth
    b2=bar(find(negative_diff), n_net_pin(negative_diff), 'FaceColor', 'r'); % Red for decline
    b1.FaceColor = 'none';       % No fill color (transparent)
    b1.EdgeColor = 'k';          % Edge color, e.g., black ('k')
    b1.LineWidth = 2;          % Optional: make the edges thicker
    b2.FaceColor = 'none';       % No fill color (transparent)
    b2.EdgeColor = 'k';          % Edge color, e.g., black ('k')
    b2.LineWidth = 2;          % Optional: make the edges thicker
    hold off;
end
set(gca, 'FontSize', 10); 
set(gcf, 'PaperPositionMode', 'auto'); 
if save_flag
    print(gcf, fileName, '-dsvg', '-r300');
    print(gcf, fileName, '-depsc', '-r0');
end

% pinning input
if pin_flag
    figure; 
name_fig = strcat(main_name,'pinning_input_'); 
fileName = fullfile(savePath, name_fig);
grid on; 
hold on;

legend_entries = gobjects(length(regions), 1);  % Preallocate legend handles

for i = 1:length(regions)
    % Plot solid line and store handle for legend
    h_main = plot(T, pinning_input(i,:), 'LineWidth', 2, 'Color', colors(i,:));
    legend_entries(i) = h_main;
end

legend(legend_entries, regions, 'Location', 'best', 'FontSize', 10);
title('Evolution over time of the Pinning Input '); 
xlabel('k', 'FontSize', 10, 'FontWeight', 'bold');
set(gca, 'FontSize', 10); 
set(gcf, 'PaperPositionMode', 'auto'); 
end
if save_flag
    print(gcf, fileName, '-dsvg', '-r300');
    print(gcf, fileName, '-depsc', '-r0');
end

% mobile populations
name_fig = strcat(main_name,'mobile_populations_'); 
fileName = fullfile(savePath, name_fig);
figure;
grid on; 
hold on;

legend_entries = gobjects(length(regions), 1);  % Preallocate legend handles

for i = 1:length(regions)
    % Plot solid line and store handle for legend
    h_main = plot(T, 100.*( n_m(i,:) - n_m(i,1) )./ (n_m(i,1)), 'LineWidth', 2, 'Color', colors(i,:));
    legend_entries(i) = h_main;

    % Plot dashed line without adding to legend
    if pin_flag
        plot(T, 100.*( n_m_pin(i,:) - n_m_pin(i,1) )./ (n_m_pin(i,1)), 'LineStyle', '--', 'LineWidth', 2, ...
             'Color', colors(i,:), 'HandleVisibility', 'off');
    end
end

legend(legend_entries, regions, 'Location', 'best', 'FontSize', 10);
title('Mobile Population Change Over Time Relative to Initial Value [%]'); 
xlabel('k', 'FontSize', 10, 'FontWeight', 'bold');
set(gca, 'FontSize', 10); 
set(gcf, 'PaperPositionMode', 'auto'); 

if save_flag
    print(gcf, fileName, '-dsvg', '-r300');
    print(gcf, fileName, '-depsc', '-r0');
end
%% 3.Calculations

[net_internal_migration_x1000_without_pin, net_balance_without_pin, n_average ]= get_net_internal_migration_x1000 (total_flux, n); 


[net_internal_migration_x1000_pin, net_balance_pin, n_average_pin ]= get_net_internal_migration_x1000 (total_flux_pin, n_pin); 

% net_internal_migration_x1000_without_pin/1000 .* population_average_ss
% net_internal_migration_x1000_pin/1000 .* population_average_ss

% observed_rates/1000 .* population_average_ss

