function [n_im_k] = get_mobile_population_i_k (n_i_k, n_i_0, alpha_i)
%get_mobile_population_i_k computes the mobile population for the i-th
%region at time k, i.e., n_im_k
%
% Input:
%   n       - population of the i-th region at time k 
%   alpha_i - 
%
% Output
%   n_im_k - Mobile population for the i-th node at time k
%

n_im_k = n_i_k - alpha_i * n_i_0;  

end