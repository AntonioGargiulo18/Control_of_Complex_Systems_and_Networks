function D =  get_distance_matrix (regions, latitudes, longitudes)

% Function to compute great-circle distance using haversine formula
radius_earth = 6371; % Earth radius in km
N = length(regions);
D = zeros(N, N);

for i = 1:N
    for j = 1:N
        if i ~= j
            lat1 = deg2rad(latitudes(i));
            lon1 = deg2rad(longitudes(i));
            lat2 = deg2rad(latitudes(j));
            lon2 = deg2rad(longitudes(j));
            
            dlat = lat2 - lat1;
            dlon = lon2 - lon1;
            
            a = sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlon/2)^2;
            c = 2 * atan2(sqrt(a), sqrt(1 - a));
            D(i,j) = round(radius_earth * c, 2);
        end
    end
end



end