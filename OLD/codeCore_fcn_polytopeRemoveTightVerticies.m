


function cleaned_polytope = codeCore_fcn_polytopeRemoveTightVerticies(polytope,tolerance)
%% main code §

% pull values
vertices = polytope.vertices;
centroid = polytope.mean;


% Work with tolerance squared, since it avoids square-root calculations
tol_squared = tolerance.^2;

% find all vertices that are larger away from next one, relative to
% tolerance
good_ind = ...
    sum((vertices(1:end-1,:)-vertices(2:end,:)).^2,2)>(tol_squared);

% Check how many are good. If only 2 points are left, then the
% result is just a line. If 1 or zero, then result is just a point.
if sum(good_ind)>2 % sufficient good points to make a shape
    new_vert = vertices(good_ind,:);
elseif sum(good_ind)==2 % line shape
    new_vert = vertices(good_ind,:);
    
    % The line may not go through the centroid, which is odd. We force
    % this by removing the point closest to the centroid
    distances_to_centroid = sum((new_vert-centroid).^2,2).^0.5;
    if distances_to_centroid(1,1)>distances_to_centroid(2,1)
        new_vert(2,:)=centroid;
    else
        new_vert(1,:)=centroid;
    end
    
    new_vert = [new_vert; flipud(new_vert)];
else % singular shape (i.e. point) or no shape
    new_vert = [centroid; centroid; centroid];
end

% adjust polytopes
cleaned_polytope.vertices = [new_vert; new_vert(1,:)];
cleaned_polytope = fcn_MapGen_fillPolytopeFieldsFromVerticies(cleaned_polytope);


% §
% Debug
%
% Functions §