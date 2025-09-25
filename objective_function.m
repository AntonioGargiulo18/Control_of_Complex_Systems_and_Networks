function error = objective_function(alpha, n_0, D, T, internal_migration_reference)

% Run the model with current alpha values
[n, ~, ~, ~, outflux, influx] = solve_continuity_equation(n_0, alpha, D, T);

% Compute net migration rates
total_flux_over_time = sum(influx - outflux, 2);
population_average_ss = (n(:, end) + n_0) / 2;
net_migration_x1000 = (total_flux_over_time ./ population_average_ss) * 1000;
predicted_rates = net_migration_x1000; 

% Compute error
error = predicted_rates - internal_migration_reference;

end