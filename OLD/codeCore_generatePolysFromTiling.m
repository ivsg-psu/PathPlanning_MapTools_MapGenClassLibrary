%% create tiling
function polytopes = INTERNAL_fcn_MapGen_generatePolysFromTiling(seed_points,V,C,stretch)

%% main code §
num_poly = size(seed_points,1);
polytopes(num_poly) = ...
    struct(...
    'vertices',[],...
    'xv',[],...
    'yv',[],...
    'distances',[],...
    'mean',[],...
    'area',[],...
    'max_radius',[]);


remove = 0; % keep track of how many cells to be removed
for poly = 1:num_poly % pull each cell from the voronoi diagram
    % x and y values from this cell
    xv = V(C{poly},1)'; 
    yv = V(C{poly},2)';
    
    verticies = V(C{poly},:);
    interior_point = seed_points(poly,:);
    
    % Are any vertices outside the [0,1] range?
    if any(xv>1) || any(yv>1) || any(xv<0) || any(yv<0)

        % Crop vertices to allowable range
        [cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies, interior_point);
        xv = cropped_vertices(1:end-1,1)';
        yv = cropped_vertices(1:end-1,2)'; 
    else
        
    end 

    % Are polytopes not trivial in length? (This may not be needed)
    if length(xv)>2                
    
        % make sure cw
        vec1 = [xv(2)-xv(1),yv(2)-yv(1),0]; % vector leading into point
        vec2 = [xv(3)-xv(2),yv(3)-yv(2),0]; % vector leading out of point
        xing = cross(vec1,vec2); % cross product of two vectors
        if sign(xing(3)) == -1 % points ordered in wrong direction
            xv = fliplr(xv);
            yv = fliplr(yv);
        end
        
        % enter info into polytope structure
        polytopes(poly-remove).vertices = [[xv xv(1)]' [yv yv(1)]']; % repeat first vertice for easy plotting
        
        polytopes(poly-remove) = fcn_MapGen_fillPolytopeFieldsFromVerticies(polytopes(poly-remove));

    else % if 2 or less points in cell 
        remove = remove+1; % skip cell and remove later
    end
end

% remove extra empty polytopes
polytopes = polytopes(1:(num_poly-remove));

% Apply the stretch
num_poly = length(polytopes);
for poly = 1:num_poly % pull each cell from the voronoi diagram
    
    polytopes(poly).vertices  = polytopes(poly).vertices.*stretch;
    polytopes(poly) = fcn_MapGen_fillPolytopeFieldsFromVerticies(polytopes(poly));
    
end % Ends for loop for stretch

% §
% Debug
%
% Functions §
end


