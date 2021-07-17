function [ ...
polytopes ...
] = ...
fcn_MapGen_generatePolysFromTiling( ...
seed_points, ...
V, ...
C, ...
AABB, ...
stretch, ...
varargin...
)
% fcn_MapGen_generatePolysFromTiling
% creates polytopes given seed points, V and C matrices from Voronoi 
% tiling, and stretch matrix
% 
% 
% 
% FORMAT:
% 
%    [ ...
%    polytopes ...
%    ] = ...
%    fcn_MapGen_generatePolysFromTiling( ...
%    seed_points, ...
%    V, ...
%    C, ...
%    AABB, ...
%    stretch, ...
%    (fig_num) ...
%    )
% 
% INPUTS:
% 
%     seed_points: the list of seed points in [x y] format, where x and y 
%     are columns
% 
%     V: the V matrix resulting from Voronoi calculations
% 
%     C: the C matrix resulting from Voronoi calculations
%
%     AABB: the axis-aligned bounding box, in format of 
%     [xmin ymin xmax ymax], wherein the resulting polytopes must be
%     bounded.
% 
%     stretch: the factor to stretch the polytopes after they are
%     calculated, in format of [xmagnification ymagnification]
% 
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a 
%     figure to be drawn showing results.
% 
% 
% OUTPUTS:
% 
%     polytopes: the resulting polytopes after converting to polytope form.
% 
% 
% DEPENDENCIES:
% 
%     fcn_MapGen_checkInputsToFunctions
%     fcn_MapGen_cropPolytopeToRange
%     fcn_MapGen_fillPolytopeFieldsFromVertices
% 
% 
% EXAMPLES:
% 
% See the script: script_test_fcn_MapGen_generatePolysFromTiling
% for a full test suite.
% 
% This function was written on 2021_07_02 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

% 
% REVISION HISTORY:
% 
% 2021_07_02 by Sean Brennan
% -- first write of function

% 
% TO DO:
% 
% -- fill in to-do items here.

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments 
flag_do_plot = 0;      % Set equal to 1 for plotting 
flag_do_debug = 0;     % Set equal to 1 for debugging 

if flag_do_debug
    fig_for_debug = 846;
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end 

%% check input arguments?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if 1 == flag_check_inputs

    % Are there the right number of inputs?
    if nargin < 5 || nargin > 6
        error('Incorrect number of input arguments')
    end

    % Check the seed_points input, make sure it is '2column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        seed_points, '2column_of_numbers');
 
    % Check the stretch input, make sure it is '2column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        stretch, '2column_of_numbers',1);

end

% Does user want to show the plots?
if  6== nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
        flag_do_plot = 1;
    end
end

%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%See: http://patorjk.com/software/taag/#p=display&f=Big&t=Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§




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

Npolys = length(polytopes);
Nvertices_per_poly = 20; % Maximum estimate
Nvertices_per_map = Npolys*Nvertices_per_poly;
all_vertices = nan(Nvertices_per_map,3);
all_neighbors = nan(Nvertices_per_map,1);
flag_contains_corner = zeros(Npolys,1);

% Loop through the polytopes, filling verticies matrix
for ith_poly = 1:Npolys
    vertices_open = V(C{ith_poly},:); 
    vertices = [vertices_open; vertices_open(1,:)]; % Close off the vertices
    Nvertices = length(vertices(:,1));
    if Nvertices>Nvertices_per_poly
        error('Need to resize the number of allowable vertices');
    else
        row_offset = (ith_poly-1)*Nvertices_per_poly;
        all_vertices(row_offset+1:row_offset+Nvertices,1) = ith_poly;
        all_vertices(row_offset+1:row_offset+Nvertices,2:3) = vertices;
    end
    
    for jth_neighbor = 1:Npolys
        if jth_neighbor~=ith_poly % Make sure not checking self
            vertices_neighbors = V(C{jth_neighbor},:); % Grab verticies
            vertices_neighbors = vertices_neighbors(~isinf(vertices_neighbors(:,1)),:);
            if any(ismember(vertices_open,vertices_neighbors)) % Any shared?
                neighbors_vector = all_neighbors(row_offset+1:row_offset+Nvertices_per_poly,:);
                if ~any(ismember(neighbors_vector,jth_neighbor))
                    next_index = find(isnan(neighbors_vector),1,'first');
                    if ~isempty(next_index)
                        neighbors_vector(next_index) = jth_neighbor;
                        all_neighbors(row_offset+1:row_offset+Nvertices_per_poly,:) = ...
                            neighbors_vector;
                    else
                        error('Need to make neighbors bigger');
                    end
                end
            end
        end
    end % Ends loop through neighbors
    
end

% % Convert cartesian to polar
% [theta_values,~] = cart2pol(all_vertices(:,2),all_vertices(:,3));
% 
% % Reshape by polytope
% theta_matrix = reshape(theta_values,Nvertices_per_poly,Npolys);
% 
% 
% all_vertices_column_no_nan = all_vertices_column(~isnan(all_vertices_column));
% 
% exterior_polytopes = INTERNAL_fcn_MapGen_findBoundaryPolytopes(all_vertices);

remove = 0; % keep track of how many cells to be removed
for poly = 1:num_poly % pull each cell from the voronoi diagram
    % Grab verticies    
    vertices_open = V(C{poly},:); 
    vertices = [vertices_open; vertices_open(1,:)]; % Close off the vertices
    interior_point = seed_points(poly,:);

    %     tolerance = 0.001;
    %     location = [0.2102 0.0506];
    %     if (...
    %             (interior_point(1,1)<location(1)+tolerance) && ...
    %             (interior_point(1,1)>location(1)-tolerance) && ...
    %             (interior_point(1,2)<location(2)+tolerance) && ...
    %             (interior_point(1,2)>location(2)-tolerance))
    %         disp('stop here');
    %     end
    
   
    % Are any vertices outside the AABB?
    if ~all(fcn_MapGen_isWithinABBB(AABB,vertices)==1)
        
        % Crop vertices to allowable range        
        [cropped_vertices] = ...
            fcn_MapGen_cropPolytopeToRange(vertices, interior_point, AABB);
        %         xv = cropped_vertices(1:end-1,1)';
        %         yv = cropped_vertices(1:end-1,2)';
    else
        cropped_vertices = [vertices; vertices(1,:)]; % Close off the vertices
        %         % x and y values from this cell
        %         xv = verticies(:,1)';
        %         yv = verticies(:,2)';
    end 

    % Are polytopes not trivial in length? (This may not be needed)
    if length(cropped_vertices(:,1))>2                
    
        % make sure cw
        vec1 = [cropped_vertices(2,:)-cropped_vertices(1,:),0]; % vector leading into point
        vec2 = [cropped_vertices(3,:)-cropped_vertices(2,:),0]; % vector leading out of point
        xing = cross(vec1,vec2); % cross product of two vectors
        if sign(xing(3)) == -1 % points ordered in wrong direction
            cropped_vertices = flipud(cropped_vertices);
        end
        
        % enter info into polytope structure
        polytopes(poly-remove).vertices = cropped_vertices; 

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
end % Ends for loop for stretch

% Fill in all the other fields
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes);
    
%§
%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if flag_do_plot
    
    figure(fig_num);
    clf;
    hold on
    scale = max(AABB,[],'all') - min(AABB,[],'all');
    new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];
    axis(new_axis);
    
    % plot the polytopes
    fcn_MapGen_plotPolytopes(polytopes,fig_num,'b',2);
    
    % plot all vertices
    plot(all_vertices(:,2),all_vertices(:,3),'c','Linewidth',1);
    
    
    % plot the seed points in red
    plot(seed_points(:,1),seed_points(:,2),'r.','Markersize',10);
    
    % plot the means in black
    temp = zeros(length(polytopes),2);
    for ith_poly = 1:length(polytopes)
        temp(ith_poly,:) = polytopes(ith_poly).mean;
    end
    plot(temp(:,1),temp(:,2),'ko','Markersize',3);
    
    % plot the connections between the polytopes
    data = [];
    for ith_poly = 1:Npolys
        row_offset = (ith_poly-1)*Nvertices_per_poly;
        neighbors_vector = all_neighbors(row_offset+1:row_offset+Nvertices_per_poly,:);
        neighbors = neighbors_vector(~isnan(neighbors_vector));
        for jth_neighbor = 1:length(neighbors)
            data = [data; seed_points(ith_poly,:); seed_points(neighbors(jth_neighbor),:); nan nan]; %#ok<AGROW>
        end
    end
    plot(data(:,1),data(:,2),'-','Linewidth',0.5);
    
end % Ends the flag_do_plot if statement    

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends the function

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _                 
%  |  ____|              | | (_)                
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___ 
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%                                               
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§




