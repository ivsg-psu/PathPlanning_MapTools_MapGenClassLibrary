clear polytopes
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];
polytopes = fillPolytopeFieldsFromVerticies(polytopes)

function filled_polytopes = fillPolytopeFieldsFromVerticies(polytopes)

%% main code §


% Apply the stretch
num_poly = length(polytopes);
for ith_poly = 1:num_poly % pull each polytope
    
    % adjust polytopes
    filled_polytopes(ith_poly).xv        = (polytopes(ith_poly).vertices(1:end-1,1)');
    filled_polytopes(ith_poly).yv        = (polytopes(ith_poly).vertices(1:end-1,2)');
    filled_polytopes(ith_poly).distances = ...
        sum((polytopes(ith_poly).vertices(1:end-1,:) - ...
        polytopes(ith_poly).vertices(2:end,:)).^2,2).^0.5;
    
    % Calculate the mean and area
    [filled_polytopes(ith_poly).mean,filled_polytopes(ith_poly).area] = ...
        fcn_MapGen_polytopeCentroidAndArea(polytopes(ith_poly).vertices);
    
    % Find max radius
    radii = sum(...
        (polytopes(ith_poly).vertices(1:end-1,:) - ...
        ones(length(polytopes(ith_poly).xv),1)*polytopes(ith_poly).mean).^2,2).^0.5;
    filled_polytopes(ith_poly).max_radius = ...
        max(radii);
end

% §
% Debug
%
% Functions §
end
