function [bounded_vertices] = ...
    fcn_MapGen_removeInfiniteVertices(...
    all_vertices,seed_points,AABB,Nvertices_per_poly, varargin)

% fcn_MapGen_removeInfiniteVertices
% removes infinite vertices created by Voronoi tiling steps
%
% FORMAT:
%
%    [bounded_vertices] = ...
%    fcn_MapGen_removeInfiniteVertices(...
%    all_vertices,seed_points,AABB,Nvertices_per_poly, (fig_num))
%
% INPUTS:
%
%     all_vertices: the list of vertex points defining the polytope map, in
%     [ID x y] format, where ID, x, and y are columns. Each polytope is
%     organized into clusters of rows, with each cluster containing
%     Nvertices_per_poly rows. Unused rows are filled with NaN.
%
%     seed_points: all seed points inside the polytopes, in [x y]
%     format, where x and y are scalars, sorted by ID
%
%     AABB: the axis-aligned bounding box, in format of
%     [xmin ymin xmax ymax], wherein the resulting polytopes must be
%     bounded.
%
%     Nvertices_per_poly: a constant specifying how many polytope
%     verticies, maximum, are in each polytope.
%
%     (optional inputs)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%     bounded_vertices: the resulting vertices of the polytope, forced to
%     fit within the AABB
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_convertAABBtoWalls
%     fcn_MapGen_isWithinABBB
%     fcn_MapGen_snapToAABB
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_removeInfiniteVertices
% for a full test suite.
%
% This function was written on 2021_07_17 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

%
% REVISION HISTORY:
%
% 2021_07_17 by Sean Brennan
% -- first write of function
% 2021_07_30 by Sean Brennan
% -- fixed errors due to corners being omitted
% 2023_02_22 by Sean Brennan
% -- switched over to DebugTools
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_07_13 by Sean Brennan
% -- added error catching for single vertex triangles

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
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
        narginchk(4,5);

        % Check the all_vertices input, make sure it has 3 columns, and can be
        % mixed NaN with numeric values
        fcn_DebugTools_checkInputsToFunctions(...
            all_vertices, '3column_of_mixed');

        % Check the Nvertices_per_poly input, make sure it has 1 column, 1 row
        fcn_DebugTools_checkInputsToFunctions(...
            Nvertices_per_poly, '1column_of_numbers',1);

        % Check the seed_points input, make sure it is '2column_of_numbers'
        % type, with correct number of rows
        fcn_DebugTools_checkInputsToFunctions(...
            seed_points, '2column_of_numbers',...
            round(length(all_vertices(:,1))/Nvertices_per_poly));

        % Check the AABB input, make sure it is '4column_of_numbers' type, with
        % 1 row
        fcn_DebugTools_checkInputsToFunctions(...
            AABB, '4column_of_numbers',1);

    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  (5 == nargin) && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
    temp = varargin{end}; % Last argument is always figure number
    if ~isempty(temp) % Make sure the user is not giving empty input
        fig_num = temp;
        flag_do_plot = 1; % Set flag to do plotting
    end
else
    if flag_do_debug % If in debug mode, do plotting but to an arbitrary figure number
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
% Goes through the vertices and removes infinite values by inserting
% points prior, and after the infinite one that "close" the polytope.

% Convert axis-aligned bounding box to wall format
walls = fcn_MapGen_convertAABBtoWalls(AABB, -1);


if flag_do_debug
    figure(fig_for_debug);
    clf;
    hold on;
    
    scale = max(AABB,[],'all') - min(AABB,[],'all');
    new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];
    axis(new_axis);
    
    
    % Plot the vertices
    plot(all_vertices(:,2),all_vertices(:,3),'r.-');
    
    % Plot the walls
    plot(walls(:,1),walls(:,2),'k-');
end



% Clean up and sort the vertices so that we can associate neighbors
all_vertices_no_nan = all_vertices(~isnan(all_vertices(:,1)),:);
sorted_all_vertices = sortrows(all_vertices_no_nan,[2 3]);

% Remove repeats
sorted_all_vertices = unique(sorted_all_vertices,'rows','stable');



bounded_vertices = all_vertices; % Initialize the array
if any(isinf(all_vertices),'all') % Are there any infinite vertices
    
    % Find the bad indices
    bad_indices = isinf(all_vertices(:,2));
    
    % Grab the polytope IDs
    bad_polytopes = all_vertices(bad_indices,1);
        
    % Loop through bad indices
    for ith_poly = 1:length(bad_polytopes)
        bad_poly = bad_polytopes(ith_poly);
        
        %         % FOR DEBUGGING:
        %         interior_point = seed_points(bad_poly,:);
        %         tolerance = 0.001;
        %         location = [0.6128 0.9867];
        %         if (...
        %                 (interior_point(1,1)<location(1)+tolerance) && ...
        %                 (interior_point(1,1)>location(1)-tolerance) && ...
        %                 (interior_point(1,2)<location(2)+tolerance) && ...
        %                 (interior_point(1,2)>location(2)-tolerance))
        %             disp('stop here');
        %         end

        
        % Grab this bad polytope
        row_offset = (bad_poly-1)*Nvertices_per_poly;
        this_poly = all_vertices(row_offset+1:row_offset+Nvertices_per_poly,:);

        % Plot the bad polytope?
        if flag_do_debug           
            % Plot the vertices
            plot(this_poly(:,2),this_poly(:,3),'b.-');
        end
        
        % Grab vertices out of this polytope, removing nan values
        vertices = this_poly(:,2:3);
        vertices = vertices(~isnan(vertices(:,1)),:);
        
        % Check how many infinities are in this polytope
        bad_indices = find(isinf(vertices(:,1)));

        
        % Check to see if this value has already been fixed
        if ~any(isinf(this_poly),'all')
            break;
        end

        % Warn user if 2 infinite values found. This may cause the code to fail
        % because it searches for the points before and after infinity,
        % assuming these points are NOT infinite.
        if length(bad_indices)>1
            % Check to see if the infinity is at start and end, artificially
            % repeated
            if isequal(bad_indices,[1; length(vertices(:,1))])
                if length(bad_indices)==3
                    disp('stop here');
                else
                    vertices_no_repeats = vertices(1:end-1,:);
                end
            else
                warning('on','backtrace');
                warning('More than 2 infinities found in one vector. Code may break');
            end
        else
            vertices_no_repeats = vertices(1:end-1,:);
        end
        
        bad_index = bad_indices(1);
        % Rearrange the points so that the infinite index is the first one.
        % Makes things easier in later steps since we don't have to carry
        % around two snips of data, just one
        vertex_string = ...
            [vertices_no_repeats(bad_index+1:end,:); ...
            vertices_no_repeats(1:bad_index-1,:)];

        if length(vertex_string(:,1))==1
            warning('on','backtrace');
            warning('A single vertex polytope was encountered. Unable to continue');
            error('single vertex "triangle" encountered');
        end
        
        if flag_do_debug
            % Plot the vertex_string
            plot(vertex_string(:,1),vertex_string(:,2),'b.-');
        end
        
        %         [cropped_vertices,~] = ...
        %             fcn_MapGen_cropVerticesByWallIntersections(vertex_string,walls);
        %
        %         if flag_do_debug
        %             % Plot the cropped_vertices
        %             plot(cropped_vertices(:,1),cropped_vertices(:,2),'g.-');
        %         end
        %
        %         prior_point = cropped_vertices(end,:);
        %         prior_point_lead_in = cropped_vertices(end-1,:);
        %         next_point = cropped_vertices(1,:);
        %         start_data = [];
        %         end_data = cropped_vertices(1:end,:);
        
        % Find the prior and next points relative to the bad index point
        % The prior_point and next_point are not used as vertices themselves,
        % but rather to create new vertices by snapping to the bounding polygon
        % or ABB. Thus, the data before and after these points, including these
        % points, must be kept.
        prior_point = vertex_string(end,:);
        prior_point_lead_in = vertex_string(end-1,:);
        next_point = vertex_string(1,:);
        next_point_lead_in = vertex_string(2,:);
        start_data = [];
        end_data = vertex_string(1:end,:);
        
        % If prior or next points are inside the AABB, need to push them
        % out to the boundaries. This is a bit challenging since these
        % represent cases where the Vornoi boundary is missing, so
        % basically it means we have to calculate the boundary ourselves.
        isInside = fcn_MapGen_isWithinABBB(AABB, [prior_point; next_point], -1);        
        if isInside(1)
            new_prior = INTERNAL_fcn_MapGen_findMissingBoundaryPoint(...
                prior_point,prior_point_lead_in, bad_poly,AABB,sorted_all_vertices,seed_points);
        else
            new_prior = prior_point;
        end
        
        if isInside(2)
            new_next = INTERNAL_fcn_MapGen_findMissingBoundaryPoint(...
                next_point,next_point_lead_in, bad_poly,AABB,sorted_all_vertices,seed_points);
        else
            new_next = next_point;
        end
        
        
        %         % Substitute data in, removing the infinite value
        %         resulting_vertices = [start_data; new_prior; new_next; end_data];
        %
        %         if flag_do_debug
        %             % Plot the resulting_vertices
        %             plot(resulting_vertices(:,1),resulting_vertices(:,2),'g.-');
        %         end
        
        if flag_do_debug
            % Plot the resulting_vertices
            plot(new_prior(:,1),new_prior(:,2),'g.');
            plot(new_next(1,1),new_next(1,2),'go');
        end
        
        % Is the point is in a corner? If so, add an extra point for the
        % corner
        extra_point = [];
        % COMMENTED OUT BECAUSE CORNER CAPTURE IS NOW IN LATER SECTIONS
        %         if all(isinf(vertices(bad_index,:))) && any(isInside)
        %             % Calculate the angles covered by the vertices
        %             end_points_offset = [new_prior; new_next] - mean([AABB(1:2);AABB(3:4)],1);
        %             [angles, ~] = cart2pol(end_points_offset(:,1),end_points_offset(:,2));
        %             min_angle = min(angles);
        %             max_angle = max(angles);
        %
        %             % Check to see if angle crosses over -180 degrees
        %             if (max_angle-min_angle)>pi
        %                 angles = mod(angles,2*pi);
        %                 min_angle = min(angles);
        %                 max_angle = max(angles);
        %             end
        %
        %             box_angle_right = abs(atan2(AABB(4)-AABB(2),AABB(3)-AABB(1)));
        %             box_angle_left  = pi-box_angle_right;
        %
        %             % Find location of the corners, and add them
        %             if min_angle<-box_angle_left && max_angle>= -box_angle_left
        %                 % Bottom left corner
        %                 extra_point = [AABB(1) AABB(2)];
        %             elseif min_angle<-box_angle_right && max_angle>= -box_angle_right
        %                 % Bottom right corner
        %                 extra_point = [AABB(3) AABB(2)];
        %             elseif min_angle<box_angle_left && max_angle>= box_angle_left
        %                 % Top left corner
        %                 extra_point = [AABB(1) AABB(4)];
        %             elseif min_angle<box_angle_right && max_angle>= box_angle_right
        %                 % Top right corner
        %                 extra_point = [AABB(3) AABB(4)];
        %             end
        %         end
        
        % Shuffle data into result
        % Make a fresh nan_vector
        non_repeating_resulting_vertices = nan(Nvertices_per_poly,3);
        
        % Create data
        fill_vertices = [start_data; new_prior; extra_point; new_next; end_data];
        
        % Close off vertices
        closed_fill_vertices = ...
            [fill_vertices; fill_vertices(1,:)];
        
        % Push into format for all_vertices
        non_repeating_resulting_vertices(1:length(closed_fill_vertices(:,1)),1) = ...
            bad_poly;
        non_repeating_resulting_vertices(1:length(closed_fill_vertices(:,1)),2:3) = ...
            closed_fill_vertices;


        bounded_vertices(row_offset+1:row_offset+Nvertices_per_poly,:) = ...
            non_repeating_resulting_vertices;        
    end % Ends loop through bad indices
     
end

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
    hold on;
    grid on;
    grid minor;
    
    scale = max(AABB,[],'all') - min(AABB,[],'all');
    new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];
    axis(new_axis);
    
    
    % Plot the vertices
    plot(all_vertices(:,2),all_vertices(:,3),'r.-','Linewidth',3);
    
    % Plot the walls
    plot(walls(:,1),walls(:,2),'k-');
    
    % Plot the seed_points
    plot(seed_points(:,1),seed_points(:,2),'r.','Markersize',10);
    
    % Plot the cropped_vertices locations
    plot(bounded_vertices(:,2),bounded_vertices(:,3),'g-','Linewidth',2);
    
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


function new_point = INTERNAL_fcn_MapGen_findMissingBoundaryPoint(...
    test_point,incoming_point, bad_poly,AABB,sorted_all_vertices,seed_points)
% The goal of this funtion is to replace a missing bounary point on a
% polytope. This usually occurs because a polytope is specified as having a
% point at infinity.

flag_do_debug = 0;
fig_for_debug = 846; %#ok<NASGU>

% Find the current and neighbor seed points by finding which ones
% have the same x values
current_seed = seed_points(bad_poly,:);

% For debugging:
tolerance = 0.001;
location = [0.9935 0.1354];
if (...
        (current_seed(1,1)<location(1)+tolerance) && ...
        (current_seed(1,1)>location(1)-tolerance) && ...
        (current_seed(1,2)<location(2)+tolerance) && ...
        (current_seed(1,2)>location(2)-tolerance))
    disp('stop here');
end


if flag_do_debug
    figure(fig_for_debug);
    clf;
    hold on;
    
    % Plot the current_seed
    plot(current_seed(:,1),current_seed(:,2),'g.','Markersize',20);
    plot(current_seed(:,1),current_seed(:,2),'b.','Markersize',10);
    plot(test_point(:,1),test_point(:,2),'bo','Markersize',10);
    plot(incoming_point(:,1),incoming_point(:,2),'bo','Markersize',5);

end 

% Each point can be connected to only three other polytopes. If there's an
% infinity there, then one of the connections is missing. To find the
% missing connections, we trace out the connections.

% Grab the neighbors by first finding all the vertices that match the test
% point
index_neighbors = ismember(sorted_all_vertices(:,2),test_point(1,1));

% Grab the ID numbers of all neighbors
neighbors = sorted_all_vertices(index_neighbors,1);

% Remove the neighbor that matches the current polytope
valid_neighbors = neighbors(neighbors~=bad_poly);


% Grab the incoming by first finding all the vertices that match the
% incoming point
index_incoming = ismember(sorted_all_vertices(:,2),incoming_point(1,1));

% Grab the ID numbers of all incoming
incoming = sorted_all_vertices(index_incoming,1);

% Remove the incoming that matches the current polytope
valid_incoming = incoming(incoming~=bad_poly);

% We can identify which vertex to create by finding which neighbor is NOT
% incoming
missing_neighbor = valid_neighbors(~ismember(valid_neighbors,valid_incoming));

% Grab the seed points of the missing_neighbor
missing_neighbors_seed = seed_points(missing_neighbor,:);

% Take midpoint between seeds. This will be used to create a vector that
% creates the missing wall between the current polytope and the other one
% that contains the missing wall.
midpoint = (current_seed+missing_neighbors_seed)/2;

if flag_do_debug
    figure(fig_for_debug);
       
    % Plot the neighbors_prior_seeds
    plot(missing_neighbors_seed(:,1),missing_neighbors_seed(:,2),'g.','Markersize',20);
    
    % Plot the midpoint
    plot(midpoint(:,1),midpoint(:,2),'k.','Markersize',10);

end

% Check which is closer - midpoint or the current point? - to where the
% snap point would have been:
points_for_vector = [midpoint;test_point];
nominal_new_boundary_point = fcn_MapGen_snapToAABB(AABB,test_point,[],-1);
distances_squared = ...
    sum((points_for_vector-nominal_new_boundary_point).^2,2);
[~,sort_index] = sort(distances_squared,'descend');
sorted_points_for_vector = points_for_vector(sort_index,:);

if length(sorted_points_for_vector(:,1))<2
    error('Insufficient points for a vector!');
end
if flag_do_debug
    figure(fig_for_debug);
    
    % Plot the vote
    plot(nominal_new_boundary_point(:,1),nominal_new_boundary_point(:,2),'mo','Markersize',10);

    % Plot the tail of the vector small
    plot(sorted_points_for_vector(1,1),sorted_points_for_vector(1,2),'k.','Markersize',20);
    
    % Plot the head of the vector big
    plot(sorted_points_for_vector(2,1),sorted_points_for_vector(2,2),'k.','Markersize',30);

end


snap_type = 2;
new_point = fcn_MapGen_snapToAABB(AABB, sorted_points_for_vector, snap_type,-1);


if flag_do_debug
    
    % Plot the new_point
    plot(new_point(:,1),new_point(:,2),'k.','Markersize',40);

end


end % Ends INTERNAL_fcn_MapGen_findMissingBoundaryPoint