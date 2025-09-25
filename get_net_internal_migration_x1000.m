function [a,b,c] = get_net_internal_migration_x1000(total_flux, n)

total_flux_over_time = sum (total_flux,2);  
population_average_ss   = (n(:,end)+n(:,1) )./2; 
a = (total_flux_over_time./population_average_ss).*1000;

b= (a/1000) .* population_average_ss; 

c = population_average_ss; 
end