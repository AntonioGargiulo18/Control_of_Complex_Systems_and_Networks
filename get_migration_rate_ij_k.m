function [r_ij_k] = get_migration_rate_ij_k (a, D, i, j)
%get_migration_rate_ij_k computes the migration rate r_ij_k (a) for a given
%input vector a, e.g., n_k.
%
% Input:
%   a   - Vector on which apply the calculation at time k.
%   D   - Matrix containing the distance in km between any nodes, e.g.,
%         D(1,:) contains the distance between the 1st node and all the
%         other.
%   i   - origin node.
%   j   - destination node.
%
% Output
%   r_ij_k - Migration Rate for a given couple of origin and destination.
%


% Build the Set N_ij
d_ij   = D(i,j);
n_l    = D(i,:)<d_ij;
n_l(i) = 0;             % To avoid considering l=i

% Calculate r_ij
if i ~= j
    num = a(i) * a(j);
    den = (a(i) + a(j) + sum(a(n_l))) * (a(i) + sum(a(n_l)));
    r_ij_k = num / den;
else
    r_ij_k = a(i) / (sum(a));
end

end
