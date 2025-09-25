function [J_ij_k] = get_migration_flux_ij_k (n_im_k, r_ij_k)
%get_migration_flux_ij_k computes the migration flux from origin i to
%destination j in the time interval [k,k+1).
%
% Input:
%   n_im_k       - mobile population of the i-th region at time k 
%   r_ij_k       - migration rate of the regions (i-th origin, j-th destination) 
%                  at time k. 
% Output
%   J_ij_k       - migration flux from origin i to destination j in the
%                  interval [k,k+1). 
%

J_ij_k = n_im_k * r_ij_k; 

end