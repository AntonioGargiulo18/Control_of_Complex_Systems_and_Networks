function [n, n_m, r, J, outflux, influx, total_flux, n_net, pin_input] = solve_continuity_equation_pinning (n_0, alpha, D, T, n_ref, pinned_nodes)
% Inputs:
%   n_0        : Initial population vector 
%   alpha      : Coefficients affecting mobile population
%   D          : Distance matrix between each region
%
% Outputs:
%   n_steady   : Steady-state population (Nx1)
%   n_history  : Population history (NxNt)

Nt = size(T,2); 
N  = size(D,2); 

r = zeros(N,N,Nt);
J = zeros(N,N,Nt);

n_m   = zeros(N,Nt);
n     = zeros(N,Nt);
n (:,1) = n_0;   

outflux = zeros(N,Nt);
influx  = zeros(N,Nt); 
delta = zeros(N,Nt); 
delta(pinned_nodes,:)=1;

for k=1:Nt
    % Mobile population 
    n_m(:,k) = n(:,k) - alpha.*n_0 ; 
    for i=1:N  
        % Migration rate and flux 
        for j=1:N
            % Migration rate/flux i -> j
            r (i,j,k) = get_migration_rate_ij_k(n(:,k), D, i, j);
            J (i,j,k) = n_m(i,k) * r(i,j,k); 
        end
    end
    outflux (:,k) = sum(J(:,:,k),2);
    influx  (:,k) = sum(J(:,:,k),1);
    
    pin_input(:,k) = delta(:,k).* (5e-1).* (n(:,k)-n_ref(:,k)); %(n(5,k)-5568025);  %K(0)=1e-1, loss of 7k people, not predicted 14k
    if k<Nt
        n(:,k+1) = n(:,k) + influx(:,k) - outflux(:,k) - pin_input(:,k); 
    end
end

total_flux = influx(:,:)-outflux(:,:); 
n_net = n(:,end) - n_0; 

end




