function [ ...
polytopes, ...
all_vertices...
] = ...
fcn_MapGen_generatePolysFromVoronoiAABB( ...
seed_points, ...
V, ...
C, ...
AABB, ...
stretch, ...
varargin...
)
% fcn_MapGen_generatePolysFromVoronoiAABB
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
%    fcn_MapGen_generatePolysFromVoronoiAABB( ...
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
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_cropPolytopeToRange
%     fcn_MapGen_fillPolytopeFieldsFromVertices
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_generatePolysFromVoronoiAABB
% for a full test suite.
%
% This function was written on 2021_07_02 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

%
% REVISION HISTORY:
%
% 2021_07_02 by Sean Brennan
% -- first write of function
% 2021_07_30 by Sean Brennan
% -- fixed errors due to corners being omitted
% 2023_03_13
% -- converted over to narginchk
% -- converted over to Debug input checking (via fcn_DebugTools_checkInputsToFunctions)
% -- fixed error with empty figure number inputs
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions


% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS");
    MATLABFLAG_MAPGEN_FLAG_DO_DEBUG = getenv("MATLABFLAG_MAPGEN_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_MAPGEN_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_MAPGEN_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
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

if (0==flag_max_speed)
    if 1 == flag_check_inputs

        % Are there the right number of inputs?
        narginchk(5,6);

        % Check the seed_points input, make sure it is '2column_of_numbers' type
        fcn_DebugTools_checkInputsToFunctions(...
            seed_points, '2column_of_numbers');

        % Check the stretch input, make sure it is '2column_of_numbers' type
        fcn_DebugTools_checkInputsToFunctions(...
            stretch, '2column_of_numbers',1);

    end
end


% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  6 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
    temp = varargin{end}; % Last argument is always figure number
    if ~isempty(temp) % Make sure the user is not giving empty input
        fig_num = temp;
        flag_do_plot = 1; % Set flag to do plotting
    end
else
    if flag_do_debug % If in debug mode, do plotting but to an arbitrary figure number
        fig = figure;
        fig_for_debug = fig.Number; %#ok<NASGU>
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



%% Initiate data structures
num_poly = size(seed_points,1);
polytopes(num_poly) = ...
    struct(...
    'vertices',[],...
    'xv',[],...
    'yv',[],...
    'distances',[],...
    'mean',[],...
    'area',[],...
    'max_radius',[],...
    'min_radius',[],...
    'mean_radius',[],...
    'radii',[],...
    'cost',[]);

Npolys = length(polytopes);
Nvertices_per_poly = 20; % Maximum estimate
Nvertices_per_map = Npolys*Nvertices_per_poly;
all_vertices = nan(Nvertices_per_map,3);
% all_neighbors = nan(Nvertices_per_map,1);

%% Loop through the polytopes, filling verticies and neighbors matrix
for ith_poly = 1:Npolys

    %     if ith_poly==50
    %         disp('Stop here');
    %     end

    vertices_open = V(C{ith_poly},:);

    % Remove ill-conditioned points by setting the ones really, really far
    % away to infinity
    scale = max(AABB,[],'all') - min(AABB,[],'all');
    center = [AABB(1)+AABB(3), AABB(2)+AABB(4)]/2;
    distances_from_center = sum((vertices_open-center).^2,2).^0.5;
    near_infinite = (distances_from_center/scale)>1E5;
    if any(near_infinite)
        vertices_open(near_infinite,:) = inf;
        % Remove repeated infinities
        vertices_open = unique(vertices_open,'rows','stable');
    end

    % Append results to close off the vector loop
    vertices = [vertices_open; vertices_open(1,:)]; % Close off the vertices

    Nvertices = length(vertices(:,1));
    if Nvertices>Nvertices_per_poly
        error('Need to resize the number of allowable vertices');
    else
        row_offset = (ith_poly-1)*Nvertices_per_poly;
        all_vertices(row_offset+1:row_offset+Nvertices,1) = ith_poly;
        all_vertices(row_offset+1:row_offset+Nvertices,2:3) = vertices;
    end



end


%% Remove infinite vertices
[bounded_vertices] = ...
    fcn_MapGen_removeInfiniteVertices(...
    all_vertices,seed_points,AABB,Nvertices_per_poly);

%% Crop vertices
remove = 0; % keep track of how many cells to be removed
for poly = 1:num_poly % pull each cell from the voronoi diagram
    % Grab verticies
    %     row_offset = (poly-1)*Nvertices_per_poly;
    %     vertices = bounded_vertices(row_offset+1:row_offset+Nvertices,2:3);
    %     vertices = vertices(~isnan(vertices(:,1)));
    indices = bounded_vertices(:,1)==poly;
    vertices = bounded_vertices(indices,2:3);
    interior_point = seed_points(poly,:);

    %     % For debugging
    %     tolerance = 0.001;
    %     location = [0.02970 0.12467];
    %     if (...
    %             (interior_point(1,1)<location(1)+tolerance) && ...
    %             (interior_point(1,1)>location(1)-tolerance) && ...
    %             (interior_point(1,2)<location(2)+tolerance) && ...
    %             (interior_point(1,2)>location(2)-tolerance))
    %         disp('stop here');
    %     end
    %
    %     if poly==50
    %         disp('stop here');
    %     end


    % Are any vertices outside the AABB? If so, must crop them
    if ~all(fcn_MapGen_isWithinABBB(AABB,vertices)==1)

      % Crop vertices to allowable range
        [cropped_vertices] = ...
            fcn_MapGen_cropPolytopeToRange(vertices, interior_point, AABB);
    else
        cropped_vertices = vertices;
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
        % polytopes(poly-remove).seed_point = interior_point;

    else % if 2 or less points in cell
        remove = remove+1; %#ok<NASGU> % skip cell and remove later
        error('2 points or less found in a polytope, must exit');
    end
end



% remove extra empty polytopes
polytopes = polytopes(1:(num_poly-remove));

% Check that all the wall corners are inside polytopes
polytopes = INTERNAL_fcn_addCorners(polytopes,seed_points,AABB);

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

    % number the polytopes at seed points
    for ith_poly = 1:length(polytopes)
        text_location = seed_points(ith_poly,:);
        text(text_location(1,1),text_location(1,2),sprintf('%.0d',ith_poly));
    end

    %     % number the polytopes at means
    %     for ith_poly = 1:length(polytopes)
    %         text_location = polytopes(ith_poly).mean;
    %         text(text_location(1,1),text_location(1,2),sprintf('%.0d',ith_poly));
    %     end

    % plot the connections between the polytope neighbors
    if 1==0
        % Clean up and sort the vertices so that we can associate neighbors
        all_vertices_no_nan = all_vertices(~isnan(all_vertices(:,1)),:);
        sorted_all_vertices = sortrows(all_vertices_no_nan,[2 3]);

        % Remove repeats
        sorted_all_vertices = unique(sorted_all_vertices,'rows','stable');

        % Remove infinities
        sorted_all_vertices = sorted_all_vertices(~isinf(sorted_all_vertices(:,2)));

        Nrealvertices = floor(length(sorted_all_vertices(:,1))/3);
        data = zeros(Nrealvertices*6,2);
        for ith_poly = 1:Nrealvertices
            row_offset = (ith_poly-1)*3;
            neighbors = sorted_all_vertices(row_offset+1:row_offset+3,1);

            for jth_neighbor = 2:length(neighbors)
                neigh_offset = (ith_poly-1)*6 + ((jth_neighbor-2)*3);
                data(neigh_offset+1:neigh_offset+3,:) = [seed_points(neighbors(1),:); seed_points(neighbors(jth_neighbor),:); nan nan];
            end
        end
        plot(data(:,1),data(:,2),'-','Linewidth',0.5);
    end

end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends the function




% Fill in neighbors - Note: this is SLOW. A much better way to do this
% would be to sort the all-verticies by the point values, then cluster
% the polytope indices together.

%         if jth_neighbor~=ith_poly % Make sure not checking self
%             vertices_neighbors = V(C{jth_neighbor},:); % Grab verticies
%             vertices_neighbors = vertices_neighbors(~isinf(vertices_neighbors(:,1)),:);
%             if any(ismember(vertices_open,vertices_neighbors)) % Any shared?
%                 neighbors_vector = all_neighbors(row_offset+1:row_offset+Nvertices_per_poly,:);
%                 if ~any(ismember(neighbors_vector,jth_neighbor))
%                     next_index = find(isnan(neighbors_vector),1,'first');
%                     if ~isempty(next_index)
%                         neighbors_vector(next_index) = jth_neighbor;
%                         all_neighbors(row_offset+1:row_offset+Nvertices_per_poly,:) = ...
%                             neighbors_vector;
%                     else
%                         error('Need to make neighbors bigger');
%                     end
%                 end
%             end
%         end
% end % Ends loop through neighbors

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

function polytopes_with_corners = INTERNAL_fcn_addCorners(polytopes,seed_points,AABB)
% This function loops through the corners of the AABB, and checks to see
% that all are within polytopes. If they are not, it finds the closest
% polytope to each missing corner (based on seed point location), finds the
% wall of the polytope that needs to be added, adds the corner in, then
% trims that polytope down to appropriate size with the corner included.

% Check that all the wall corners are inside polytopes
walls = fcn_MapGen_convertAABBtoWalls(AABB);
test_points = walls(1:4,:);
all_found = zeros(length(test_points(:,1)),1); % keep track of which vertices are hit
for poly = 1:length(polytopes)
    vertices = polytopes(poly).vertices;
    in_polytope = inpolygon(test_points(:,1),test_points(:,2),vertices(:,1),vertices(:,2));
    all_found = all_found + in_polytope;
end
missing_vertices = test_points(all_found==0,:);


% Loop through vertices that are NOT in polytopes, putting them in closest
% one
polytopes_with_corners = polytopes; % Initialize output
for ith_missing = 1:length(missing_vertices(:,1))
    missing_point = missing_vertices(ith_missing,:);

    % Find closest seed point
    distances_squared = sum((missing_point - seed_points).^2,2);
    [~,closest_poly] = min(distances_squared);
    vertices = polytopes(closest_poly).vertices;
    interior_point = seed_points(closest_poly,:);

    % Find the polytope wall that is closest to the missing point
    [~,~,wall_that_was_hit] = ...
        fcn_MapGen_findIntersectionOfSegments(...
        vertices(1:end-1,:),...  % wall start
        vertices(2:end,:),...    % wall end
        missing_point,...        % sensor_vector_start
        interior_point,...       % sensor_vector_end
        0);

    % Put the point into the vertices
    shoved_vertices = [...
        vertices(1:wall_that_was_hit,:); ...
        missing_point; ...
        vertices(wall_that_was_hit+1:end,:)];

    % Crop vertices to allowable range
    [cropped_vertices] = ...
        fcn_MapGen_cropPolytopeToRange(shoved_vertices, interior_point, AABB);

    % make sure cw
    vec1 = [cropped_vertices(2,:)-cropped_vertices(1,:),0]; % vector leading into point
    vec2 = [cropped_vertices(3,:)-cropped_vertices(2,:),0]; % vector leading out of point
    xing = cross(vec1,vec2); % cross product of two vectors
    if sign(xing(3)) == -1 % points ordered in wrong direction
        cropped_vertices = flipud(cropped_vertices);
    end

    % re-enter info into polytope structure
    polytopes_with_corners(closest_poly).vertices = cropped_vertices;
end
end
