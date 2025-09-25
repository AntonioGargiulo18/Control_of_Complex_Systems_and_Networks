clc; clear all; close all; 
currentDir = pwd;
savePath = fullfile(currentDir, 'figures');


%italy = shaperead('gadm41_ITA_0.shp','UseGeoCoords',true);
%italy_reg =  shaperead('gadm41_ITA_1.shp','UseGeoCoords',true);

%% 1. Preliminary data

region_number = 3; 
switch region_number 
    case 3
        % Only 3 Italian regions
        regions = {'Sicily', 'Campania', 'Lombardy'};
        latitudes  = [38.11, 40.85,  45.46];  % Palermo, Napoli,  Milano
        longitudes = [13.36, 14.27,  9.19]; 

        D = get_distance_matrix (regions, latitudes, longitudes); 
    case 20
        % All italian regions
       regions = { ...
    'Sicily', ...                % Palermo
    'Calabria', ...              % Catanzaro
    'Sardinia', ...              % Cagliari
    'Basilicata', ...            % Potenza
    'Campania', ...              % Naples
    'Apulia', ...                % Bari
    'Molise', ...                % Campobasso
    'Lazio', ...                 % Rome
    'Abruzzo', ...               % L''Aquila
    'Umbria', ...                % Perugia
    'Marche', ...                % Ancona
    'Tuscany', ...               % Florence
    'Liguria', ...               % Genoa
    'Emilia-Romagna', ...        % Bologna
    'Piedmont', ...              % Turin
    'Veneto', ...                % Venice
    'Lombardy', ...              % Milan
    'Friuli-Venezia Giulia', ...  % Trieste
    'Aosta Valley', ...          % Aosta
    'Trentino-Alto Adige'};       % Trento

latitudes = [38.11, ...   % Sicily (Palermo)
             38.91, ...   % Calabria (Catanzaro)
             39.22, ...   % Sardinia (Cagliari)
             40.64, ...   % Basilicata (Potenza)
             40.85, ...   % Campania (Naples)
             41.12, ...   % Apulia (Bari)
             41.56, ...   % Molise (Campobasso)
             41.90, ...   % Lazio (Rome)
             42.35, ...   % Abruzzo (L'Aquila)
             43.11, ...   % Umbria (Perugia)
             43.62, ...   % Marche (Ancona)
             43.77, ...   % Tuscany (Florence)
             44.41, ...   % Liguria (Genoa)
             44.50, ...   % Emilia-Romagna (Bologna)
             45.07, ...   % Piedmont (Turin)
             45.44, ...   % Veneto (Venice)
             45.46, ...   % Lombardy (Milan)
             45.65, ...   % Friuli-Venezia Giulia (Trieste)
             45.74, ...   % Aosta Valley (Aosta)
             46.07];      % Trentino-Alto Adige (Trento)

longitudes = [13.36, ...   % Sicily (Palermo)
              16.59, ...   % Calabria (Catanzaro)
              9.11,  ...   % Sardinia (Cagliari)
              15.80, ...   % Basilicata (Potenza)
              14.27, ...   % Campania (Naples)
              16.87, ...   % Apulia (Bari)
              14.66, ...   % Molise (Campobasso)
              12.48, ...   % Lazio (Rome)
              13.40, ...   % Abruzzo (L'Aquila)
              12.39, ...   % Umbria (Perugia)
              13.52, ...   % Marche (Ancona)
              11.25, ...   % Tuscany (Florence)
              8.93,  ...   % Liguria (Genoa)
              11.34, ...   % Emilia-Romagna (Bologna)
              7.69,  ...   % Piedmont (Turin)
              11.88, ...   % Veneto (Venice)
              9.19,  ...   % Lombardy (Milan)
              13.77, ...   % Friuli-Venezia Giulia (Trieste)
              7.32,  ...   % Aosta Valley (Aosta)
              11.12];      % Trentino-Alto Adige (Trento)

        D = get_distance_matrix (regions, latitudes, longitudes); 
end

colorSouth = [1, 0, 0];  % Red (for Sicily)
colorNorth = [0, 0, 1];  % Blue (for Lombardy)
colors = zeros(length(latitudes), 3);
minLat = min(latitudes);
maxLat = max(latitudes);
for i = 1:length(latitudes)
    t = (latitudes(i) - minLat) / (maxLat - minLat);  % Normalize to [0, 1]
    colors(i, :) = (1 - t) * colorSouth + t * colorNorth;
end



% Distances in km between any capital of the italian regions
% D = [ 0   312  886;  % Sicilia -> Campania, Lombardia
%      312    0  658;  % Campania -> Sicilia, Lombardia
%      886  658    0]; % Lombardia -> Sicilia, Campania
% Create and visualize the graph
G = digraph(D, regions); %Note that without the "upper" option, there will be considered two edges for each nodes in both directions. 


name_fig = strcat(num2str(region_number),'_regions_graph'); 
fileName = fullfile(savePath, name_fig);

switch region_number
    case 3
        figure;
        h = plot(G, 'XData', longitudes, 'YData', latitudes, 'EdgeLabel', G.Edges.Weight, ...
            'NodeLabel', regions, 'MarkerSize', 10, 'NodeColor', 'c', 'LineWidth', 1.5);
    case 20
        figure('Position',[100 100 1200 800]); % bigger figure
        latlim = [35 48]; % approximate lat limits for Italy
        lonlim = [6 19];  % approximate lon limits for Italy
        
        % Set up Mercator map projection
        %axesm('mercator','MapLatLimit',latlim,'MapLonLimit',lonlim);
        
        % Plot the country boundary
        %geoshow(italy, 'FaceColor',[0.9 0.9 0.9], 'EdgeColor','k');
        %hold on;
        
        % Plot the graph on top of the map
        h = plot(G, 'XData', longitudes, 'YData', latitudes, ...
            'NodeLabel', regions, ...
            'EdgeLabel', [], ...  % Hide edge labels for clarity
            'MarkerSize', 8, ...
            'NodeColor', 'c', ...
            'LineWidth', 0.5, ...
            'EdgeColor', [0.4 0.4 0.4]);
        
        % Title and formatting
        title('Graph of Italian Regions (Lat-Long Positions)');
        axis equal;
        axis off;
        grid on;
        box on;
        
        % Make edges partially transparent
        h.EdgeAlpha = 0.3;
        
        % Color nodes by latitude
        h.NodeCData = latitudes;
        colormap jet;
        colorbar; % To see the scale of latitudes
end

title('Graph of all Italian regions'); 
set(gca, 'FontSize', 10); % Set axis font size for better visibility
set(gcf, 'PaperPositionMode', 'auto'); 
print(gcf, fileName, '-dsvg', '-r300');
print(gcf, fileName, '-depsc', '-r0');
%% 2. Non-linear model and constitutive equations
% Initial conditions for the number of population

switch region_number
    case 3
        n_0 = [4779371; 5575025; 10035481];  %year 2024 from ISTAT website
        observed_rates  = [-2.8; -3.3; 1.3];  %year 2024 from ISTAT website  
        reference_rates = [-2.8; -1.5; 1.3];  %reference alpha used for the first control attempt with three regions
    case 20
        n_0 = [
                    4779371;   % Sicily
                    1832147;   % Calabria
                    1561339;   % Sardinia
                    529897;    % Basilicata
                    5575025;   % Campania
                    3874166;   % Apulia
                    287966;    % Molise
                    5710272;   % Lazio
                    1268430;   % Abruzzo
                    851954;    % Umbria
                    1481252;   % Marche
                    3660834;   % Tuscany
                    1509908;   % Liguria
                    4465678;   % Emilia-Romagna
                    4255702;   % Piedmont
                    4851851;   % Veneto
                    10035481;  % Lombardy
                    1194095;   % Friuli-Venezia Giulia
                    122714;    % Aosta Valley
                    1086095    % Trentino-Alto Adige
              ];
        observed_rates = [
                    -2.3;   % Sicily
                    -4.3;   % Calabria
                    -0.1;   % Sardinia
                    -4.8;    % Basilicata
                    -3.3;   % Campania
                    -2;   % Apulia
                    -3.8;    % Molise
                    -0.1;   % Lazio
                    +0.3;   % Abruzzo
                    +0.4;    % Umbria
                    +0.7;   % Marche
                    +1.0;   % Tuscany
                    +1.8;   % Liguria
                    +2.4;   % Emilia-Romagna
                    +2.2;   % Piedmont
                    +1.2;   % Veneto
                    +1.7;  % Lombardy
                    +1.9;   % Friuli-Venezia Giulia
                    +1.9;    % Aosta Valley
                    +1.8    % Trentino-Alto Adige
              ];
            
        reference_rates = [
                    -2.3;   % Sicily
                    -4.3;   % Calabria
                    -0.1;   % Sardinia
                    -4.8;    % Basilicata
                    -3.3+1.7;   % Campania
                    -2;   % Apulia
                    -3.8;    % Molise
                    -0.1;   % Lazio
                    +0.3;   % Abruzzo
                    +0.4;    % Umbria
                    +0.7;   % Marche
                    +1.0;   % Tuscany
                    +1.8;   % Liguria
                    +2.4;   % Emilia-Romagna
                    +2.2;   % Piedmont
                    +1.2;   % Veneto
                    +1.7;  % Lombardy
                    +1.9;   % Friuli-Venezia Giulia
                    +1.9;    % Aosta Valley
                    +1.8    % Trentino-Alto Adige
              ];

end


N = size(n_0,1); 
alpha_initial = 0.5.*ones(N,1);
T = 0:20; 

lb = zeros(size(alpha_initial)) + eps; % Avoid alpha = 0
ub = ones(size(alpha_initial));
options = optimoptions('lsqnonlin', 'Display', 'iter', 'FunctionTolerance', 1e-12, 'StepTolerance',0.001);
alpha_obs = lsqnonlin(@(alpha_initial) objective_function(alpha_initial, n_0, D, T, observed_rates), ...
                      alpha_initial, lb, ub, options);

% Let's control the migration tendency by calculating optimal alpha values
alpha_opt = lsqnonlin(@(alpha_initial) objective_function(alpha_initial, n_0, D, T, reference_rates), ...
                      alpha_initial, lb, ub, options);

% single tuning
% Index of the component to optimize
i_free     = 2;             
alpha_full = alpha_obs;  % store the fixed "base" vector

% Build a 1-D objective that only takes the varying scalar and reinjects it:
fun_wrap = @(a_scalar) ...
    objective_function( ...
        [ alpha_full(1:i_free-1);
          a_scalar;
          alpha_full(i_free+1:end) ], ...
        n_0, D, T, reference_rates );

% Initial guess + bounds for that one free parameter
x0   = alpha_full(i_free);
lb_i = lb;
ub_i = ub;

% Run the scalar optimization
a_opt = lsqnonlin(fun_wrap, x0, lb_i, ub_i, options);

% Reconstruct the full optimized vector
alpha_opt        = alpha_full;
alpha_opt(i_free)= a_opt;


%% RUN MODEL

alpha_obs=alpha_opt;
[n, n_m, r, J, outflux, influx, total_flux, n_net] = solve_continuity_equation (n_0, alpha_obs, D, T);

%[n_pin, n_m_pin, r_pin, J_pin, outflux_pin, influx_pin, total_flux_pin, n_net_pin] = solve_continuity_equation (n_0, alpha_opt, D, T);

pinned_nodes = 1:5; 
n_ref = zeros(N,numel(T)); 
for ii=1:numel(T)
    n_ref(pinned_nodes,ii) = n_0(pinned_nodes); 
    
end
[n_pin, n_m_pin, r_pin, J_pin, outflux_pin, influx_pin, total_flux_pin, n_net_pin, pinning_input] = solve_continuity_equation_pinning (n_0, alpha_obs, D, T, n_ref, pinned_nodes);



%% PLOTS
save_flag =0; 
pin_flag = 0; 
lw = 2; 
lw = 1.5;
legend_flag = 0;
main_name = '20_regions_pinning_control_campania_k_0_5'; 
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
             'LineStyle', '--', 'LineWidth', lw, 'Color', colors(i,:), ...
             'HandleVisibility', 'off');
    end
end

if legend_flag == 1
    legend(legend_entries, regions, 'Location', 'best', 'FontSize', 10);
end
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
    h_main = plot(T, total_flux(i,:), 'LineWidth', lw, 'Color', colors(i,:));
    legend_entries(i) = h_main;

    % Plot dashed line without adding to legend
    if pin_flag
        plot(T, total_flux_pin(i,:), 'LineStyle', '--', 'LineWidth', lw, ...
             'Color', colors(i,:), 'HandleVisibility', 'off');
    end
end
if legend_flag == 1
    legend(legend_entries, regions, 'Location', 'best', 'FontSize', 10);
end
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
ylim([-2*1e4,1.5*1e4]);
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
    h_main = plot(T, pinning_input(i,:), 'LineWidth', lw, 'Color', colors(i,:));
    legend_entries(i) = h_main;
end
if legend_flag == 1
    legend(legend_entries, regions, 'Location', 'best', 'FontSize', 10);
end
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
    h_main = plot(T, 100.*( n_m(i,:) - n_m(i,1) )./ (n_m(i,1)), 'LineWidth', lw, 'Color', colors(i,:));
    legend_entries(i) = h_main;

    % Plot dashed line without adding to legend
    if pin_flag
        plot(T, 100.*( n_m_pin(i,:) - n_m_pin(i,1) )./ (n_m_pin(i,1)), 'LineStyle', '--', 'LineWidth', 2, ...
             'Color', colors(i,:), 'HandleVisibility', 'off');
    end
end
if legend_flag == 1
    legend(legend_entries, regions, 'Location', 'best', 'FontSize', 10);
end
title('Mobile Population Change Over Time Relative to Initial Value [%]'); 
xlabel('k', 'FontSize', 10, 'FontWeight', 'bold');
set(gca, 'FontSize', 10); 
set(gcf, 'PaperPositionMode', 'auto'); 

if save_flag
    print(gcf, fileName, '-dsvg', '-r300');
    print(gcf, fileName, '-depsc', '-r0');
end
%% 3.Rates comparis
% n0=n_0; 
% n0_norm = (n0 - min(n0)) ./ (max(n0) - min(n0));
% figure;
% plot(alpha_obs, 'LineWidth', 2, 'DisplayName','Alpha observed');
% hold on; 
% plot(n0_norm,'DisplayName','Normalized Initial Population');
% legend show;
% %xlabel('FontSize', 10, 'FontWeight', 'bold');
% title('Alpha coefficients and normalized initial population for all regions'); 
% set(gca, 'FontSize', 10); 
% set(gcf, 'PaperPositionMode', 'auto'); 
% xticks(1:length(regions)); 
% xticklabels(regions); 
% xtickangle(45); 
% print(gcf, fileName, '-dsvg', '-r300');
% print(gcf, fileName, '-depsc', '-r0');

[predicted_rates_without_pin, net_balance_without_pin, n_average ]= get_net_internal_migration_x1000 (total_flux, n); 
figure;
name_fig = strcat(main_name,'net_migration_balance_bar'); 
fileName = fullfile(savePath, name_fig);
bar(predicted_rates_without_pin,  'LineWidth', 2); 
title('Net Internal Migration Balance'); 
xticks(1:length(regions)); 
xticklabels(regions); 
xtickangle(45); 
hold on;
positive_diff = predicted_rates_without_pin > 0;
negative_diff = predicted_rates_without_pin < 0;
bar(find(positive_diff), predicted_rates_without_pin(positive_diff), 'FaceColor', 'g'); % Green for growth
bar(find(negative_diff), predicted_rates_without_pin(negative_diff), 'FaceColor', 'r'); % Red for decline
hold off;
if pin_flag
hold on;
    positive_diff = observed_rates > 0;
    negative_diff = observed_rates < 0;
    b1=bar(find(positive_diff), observed_rates(positive_diff), 'FaceColor', 'g'); % Green for growth
    b2=bar(find(negative_diff), observed_rates(negative_diff), 'FaceColor', 'r'); % Red for decline
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
    display('totor')
end

%[net_internal_migration_x1000_pin, net_balance_pin, n_average_pin ]= get_net_internal_migration_x1000 (total_flux_pin, n_pin); 

% net_internal_migration_x1000_without_pin/1000 .* population_average_ss
% net_internal_migration_x1000_pin/1000 .* population_average_ss

% observed_rates/1000 .* population_average_ss

