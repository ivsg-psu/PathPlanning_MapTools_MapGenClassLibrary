function [new_vertices, new_projection_vectors, cut_distance] = ...
    fcn_MapGen_polytopeFindVertexSkeleton(vertices, varargin)
% Calculates the VertexSkeleton for a polytope, i.e. where the vertices
% would land if the polytope were shrunk.
%
% FORMAT:
%
% [new_vertices, new_projection_vectors, cut_distance] = ...
% fcn_MapGen_polytopeFindVertexSkeleton(vertices, varargin)
%
% INPUTS:
%
%     vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%         the number of the individual polytope vertices
%
%    (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
% OUTPUTS:
%
%     new_vertices: a cell array of dimension M, where each index 1:M
%     stores a N x 2 array of the coordinates of the nested shape inside,
%     with M = 1 being the starting shape (and with dimension K+1, where K
%     is number of vertices given), and N being smaller and smaller for
%     each M value.
%
%     new_projection_vectors: a cell array of M, where each index 1:M
%     stores a N x 2 array of the unit vectors that point
%     away from the vertices of the nested shape inside, with M = 1 being
%     the starting unit vectors and N being smaller and smaller for each M value.
%
%     cut_distance: an array of 1 x M, starting from 0 for M(1) to the
%     maximum cut distance that can be used, at M(end)
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     INTERNAL_fcn_findUnitDirectionVectors
%     fcn_MapGen_polytopeFindVertexAngles
%
% % EXAMPLES:
% For additional examples, see: script_test_fcn_MapGen_polytopeFindVertexSkeleton
%
% This function was written on 2022_02_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%

% Revision History:
% 2022_02_13 - S. Brennan
% -- first write of code
% -- pulled the function out of edge shrinking code
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
if (nargin==2 && isequal(varargin{end},-1))
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
    if flag_check_inputs
        % Are there the right number of inputs?
        if nargin < 1 || nargin > 2
            error('Incorrect number of input arguments')
        end

        % Check the vertices input
        fcn_DebugTools_checkInputsToFunctions(...
            vertices, '2column_of_numbers');

    end
end


% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  2 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize results
iteration = 1; % Initialize iterations to 1
flag_stop_loop = 0;  % Set stop flag to 0 (e.g. "keep going!")
total_cut = 0;  % The total cut distance is set to zero
working_vertices = vertices;  % The vertices we are now using are the starting vertices

% Loop until the stop flag is set, e.g. loop "inward" to calculate the
% skeleton
while 0 == flag_stop_loop
    % Save the current iteration's results
    new_vertices{iteration} = working_vertices; %#ok<AGROW>
    cut_distance(iteration,1) = total_cut; %#ok<AGROW>

    % Check how many vertices we have
    Nvertices = length(working_vertices(:,1));
    if 2==Nvertices  % If 2, then no way to project, since this is the same point - so this needs to be last iteration
        new_projection_vectors{iteration} = [0 0; 0 0]; %#ok<AGROW>
        flag_stop_loop = 1;
    else % Need to calculate inward vector projection
        % Find the unit vectors that point inward from each vertex point.
        % While at it, also pull out the half-angles (half of each vertex's
        % internal angle), the edge distance from vertex to vertex, and the
        % unit vectors from vertex to vertex. The last 2 outputs are used
        % to calculate the projection location inward, after we find
        % intersection distances, to find the actual intersection points.
        if flag_do_debug  % Show what is happening
            [vector_direction_of_unit_cut, ...
                half_angles,...
                distances_vertex_to_vertex,...
                unit_vectors_vertex_to_vertex] = ...
                INTERNAL_fcn_findUnitDirectionVectors(working_vertices,fig_for_debug);
        else % Dont show what is happening
            [vector_direction_of_unit_cut, ...
                half_angles,...
                distances_vertex_to_vertex,...
                unit_vectors_vertex_to_vertex] = ...
                INTERNAL_fcn_findUnitDirectionVectors(working_vertices);
        end

        % Save the results into the projection vectors, being sure to
        % repeat the first row to last to match the points.
        new_projection_vectors{iteration} = [vector_direction_of_unit_cut; vector_direction_of_unit_cut(1,:)] ; %#ok<AGROW>

        % Next, do we need to find where the points move?
        if 3==Nvertices % If 3, then this is a line segment with only 2 points and only need to average points.


            % Fill in values for last iteration, which will just be the average
            % of the last points remaining (the first two rows of points). Note
            % that, since this is a single point, need to repeat it twice since
            % first point must match the last point.
            iteration = iteration + 1;
            new_vertices{iteration} = [mean(working_vertices(1:2,:),1); mean(working_vertices(1:2,:),1)]; %#ok<AGROW>

            % Calculate the final cut distance using the first row to represent
            % the last point.
            last_cut_dist = sum((new_vertices{iteration}(1,:) - working_vertices(1,:)).^2,2).^0.5; % Calculate this last cut distance
            cut_distance(iteration) = total_cut+last_cut_dist; %#ok<AGROW>

            % Check that final cut distance is as expected, within tolerance
            if cut_distance(iteration)>(max_cut_depth + eps*1000) || cut_distance(iteration)<(max_cut_depth - eps*1000)
                st = dbstack;
                warning('on','backtrace');
                warning('Within function: %s, in file: %s',st(1).name,st(1).file);
                warning('The predicted cut distance of %f does not match actual cut distance, %f. The results may be in error.\n',cut_distance(iteration),max_cut_depth);
            end
            new_projection_vectors{iteration} = [0 0; 0 0]; %#ok<AGROW>

            flag_stop_loop = 1;
        else  % More than 2 points, so need to shrink inward

            % Find the projection points. To do this, we use the half-angle
            % results, do some trig (see documentation), then figure out where
            % on each edge of the polytope the projection of the sides inward
            % would "hit" each other at the same location as they would hit the
            % vertex, when the edges meet. We call this the "projection point"
            % since this is where we want to project the edges inward to
            % calculate the intersections of edges with each other and the
            % vertex simultaneously.
            looped_half_angles = [half_angles; half_angles(1)]; %These are all the half-angles
            theta2s = looped_half_angles(2:end); % These are the half angles from point N+1
            theta1s = looped_half_angles(1:end-1); % These are the half angles from point N

            % Calculate A and L1 values (see documentation)
            A = tan(theta2s)./tan(theta1s);
            L1 = distances_vertex_to_vertex.*(A./(A+1));
            projection_points = working_vertices(1:end-1,:) + L1.*unit_vectors_vertex_to_vertex;

            if flag_do_debug
                figure(fig_for_debug);
                grid on
                grid minor
                hold on
                axis equal

                plot(projection_points(:,1),projection_points(:,2),'b.','Markersize',20);
            end

            % We next rotate the unit vectors from one corner to another by 90
            % degrees, and starting from the projection point, move inward by
            % Lcut. This predicts the intersection points of adjacent vertices.
            % Thus, these points predect the shape of the polytope once the
            % vertices first merge.
            Lcuts = distances_vertex_to_vertex.*tan(theta1s).*tan(theta2s)./(tan(theta1s)+tan(theta2s));
            projection_directions = unit_vectors_vertex_to_vertex*[0 1; -1 0]; % Rotate by 90 degrees so that the direction is orthogonal
            intersection_points = projection_points + Lcuts.*projection_directions;

            if flag_do_debug
                plot(intersection_points(:,1),intersection_points(:,2),'g.','Markersize',20);
                % this could be used to find interior vertex normal vectors
            end

            % Find the tightest cut that is possible, and use this to update
            % the total cut. Only the smallest shrink distance should be used
            % since this is the first internal "new" polytope created once the
            % previous one is shrunk.
            min_cut = min(Lcuts);
            total_cut = total_cut + min_cut;

            % If on the first iteration, check the maximum cut depth
            if 1 == iteration
                max_cut_depth = max(Lcuts);
            end

            %% Merge points together, if they are close to each other
            % The minimum cut distance defines how far we can cut inward before
            % vertices merge. The indices that merge have to have the same cut
            % distance and they have follow one after the other. So we search
            % for vertices that have both properties - these will have to be
            % merged (which causes a loss of vertices in the shrunk polytope).

            % Sort by the cut size, finding the cuts that are nearly exactly
            % the same cut length
            indices_repeated = find(Lcuts<(min_cut+1E5*eps));

            % Tag the vertices that are merged
            all_indices = (1:Nvertices)';
            indices_following = indices_repeated+1;
            % indices_following = mod(indices_repeated,Nvertices-1)+1;

            vertices_merged = union(indices_repeated,indices_following);
            vertices_not_merged = ~ismember(all_indices(1:end-1),vertices_merged);

            % Associate the indices with the repeats. This is done by saying
            % that any that is a "following" index should have the same value
            % as the one prior to it.
            for ith_repeat = 1:length(indices_following)
                current_index = indices_following(ith_repeat);
                all_indices(current_index) = all_indices(current_index-1);
            end

            % Check for the rollover condition, since the last point is first
            % point. In other words, the for loop above won't catch the
            % situation where the last value is the "previous" value to the
            % first value
            if all_indices(end)~=Nvertices
                first_indices = all_indices==1;
                all_indices(first_indices) = all_indices(end);
            end

            % Crop back all indices to avoid rollover end point. At this point,
            % all the vertices merged are indexed only from 1 to N-1 (in other
            % words, the polytope must shrink by at least one vertex).
            all_indices = all_indices(1:(Nvertices-1),:);
            vertices_merged = vertices_merged(vertices_merged<Nvertices);

            % Calculate the movements of the vertices, e.g. map old vertices to
            % where the new ones are located
            moved_vertex_locations = intersection_points;
            moved_vertex_locations(vertices_merged,:) = intersection_points(all_indices(vertices_merged),:);
            moved_vertex_locations(vertices_not_merged,:) = ...
                working_vertices(vertices_not_merged,:)+...
                vector_direction_of_unit_cut(vertices_not_merged,:)*min_cut;

            if flag_do_debug
                plot(moved_vertex_locations(:,1),moved_vertex_locations(:,2),'r.','Markersize',20);
            end

            % Group unique unit_direction_of_cut vectors
            unique_new_points = unique(indices_repeated);
            for ith_unique = 1:length(unique_new_points)
                current_point = unique_new_points(ith_unique);
                vertices_same = all_indices==current_point;
                vector_direction_of_unit_cut(current_point,:) = ...
                    mean(vector_direction_of_unit_cut(vertices_same,:),1);
            end

            % Move the vertices
            moved_vertices = [];
            moved_unit_vectors = [];
            for ith_vertex = 1:(Nvertices-1)
                if all_indices(ith_vertex)==ith_vertex
                    moved_vertices = [moved_vertices; moved_vertex_locations(ith_vertex,:)]; %#ok<AGROW>
                    moved_unit_vectors = [moved_unit_vectors; vector_direction_of_unit_cut(ith_vertex,:)];                 %#ok<AGROW>
                end
            end

            % Assign new vertices, remembering to wrap around
            working_vertices = [moved_vertices; moved_vertices(1,:)];

        end % Ends if-statement check if only 3 vertices, or more
    end % Ends Nvertices if statement
    iteration = iteration+1;

end % Ends while loop

%% Plot results?
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
    fcn_MapGen_polytopePlotSkeleton(new_vertices, new_projection_vectors, cut_distance,fig_num);
end % Ends flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % ends fucntion fcn_MapGen_polytopeFindVertexSkeleton


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


function [...
    vector_direction_of_unit_cut, ...
    half_angles,...
    distances_vertex_to_vertex,...
    unit_vectors_vertex_to_vertex] = ...
    INTERNAL_fcn_findUnitDirectionVectors(vertices,varargin)
% find the vector_direction_of_unit_cut to use out of each vertex point,
% e.g. the direction and distance needed to move each point, given a
% unit edge cut

% Revision History:
% 2021_08_06 - S. Brennan
% -- first write of code

% TO DO
% -- none

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    fig_for_debug = 5168;
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

if flag_check_inputs
    % Are there the right number of inputs?
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end

    % Check the vertices input
    fcn_DebugTools_checkInputsToFunctions(...
        vertices, '2column_of_numbers');

end


% Does user want to show the plots?
if  2 == nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig_for_debug = 1595;
        flag_do_plot = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate the angles. There will be 1 less angle than vertices, since
% last vertex is repeated
[angles, in_vectors, out_vectors] = fcn_MapGen_polytopeFindVertexAngles(...
    vertices);

% Confirm that all angles are positive
if ~all(angles>=0)
    error('All vertices must be organized counter-clockwise, e.g. with positive cross-products');
end

% Do preliminary calculations to determine distances and internal angles
internal_angles = pi - angles;
half_angles = internal_angles/2; % used to calculate distance to move the point
scale_factor = 1./sin(half_angles); % distance needed to move each point, given a unit edge cut

% find projection vectors pointing toward each new vertex
mean_vectors = (out_vectors-in_vectors)/2;
length_mean_vectors = sum(mean_vectors.^2,2).^0.5;
unit_direction_of_cut = mean_vectors./length_mean_vectors;

% find distances and unit vectors from vertex to vertex
difference_vertex_to_vertex = vertices(2:end,:)-vertices(1:end-1,:);
distances_vertex_to_vertex = sum(difference_vertex_to_vertex.^2,2).^0.5;
unit_vectors_vertex_to_vertex = difference_vertex_to_vertex./distances_vertex_to_vertex;

vector_direction_of_unit_cut = unit_direction_of_cut.*scale_factor;

%% Plot results?
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
    grid on
    grid minor
    hold on
    axis equal

    % Plot the polytope in red
    plot(vertices(:,1),vertices(:,2),'r-','Linewidth',2);

    % Find size of vertices
    size = max(max(vertices)) - min(min(vertices));
    nudge = size*0.003;

    % Number the midpoints of vertices with distances
    midpoints = (vertices(2:end,:)+vertices(1:end-1,:))/2;
    for ith_midpoint = 1:length(midpoints(:,1))
        text(midpoints(ith_midpoint,1)+nudge,midpoints(ith_midpoint,2),...
            sprintf('%.2f',distances_vertex_to_vertex(ith_midpoint)));
    end

    % Label the vertices with their angles
    for ith_angle = 1:length(angles(:,1))
        ith_vertex = ith_angle;
        text(vertices(ith_vertex,1)+nudge,vertices(ith_vertex,2),...
            sprintf('%.0d = %.0f deg',ith_vertex, angles(ith_angle,1)*180/pi));
    end

    % Draw the unit vectors in the cut direction
    quiver(vertices(1:end-1,1),vertices(1:end-1,2),unit_direction_of_cut(:,1),unit_direction_of_cut(:,2),0);
    % this could give interior vertex normal vectors
    % Draw the unit vectors in the vertex direction
    quiver(vertices(1:end-1,1),vertices(1:end-1,2),unit_vectors_vertex_to_vertex(:,1),unit_vectors_vertex_to_vertex(:,2),0);
    % this could give interior vertex normal vectors
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends INTERNAL_fcn_findUnitDirectionVectors



function [h_fig] = ...
    fcn_MapGen_polytopePlotSkeleton(vertices, projection_vectors, cut_distance,varargin)
% plots the skelton of a polytope, e.g. where the polytope will shrink if
% the edges are all brought in at exactly the same rate.
%
% FORMAT:
%
% fcn_MapGen_polytopePlotSkeleton(vertices, projection_vectors, cut_distance, (fig_num))
%
% INPUTS:
%
%     vertices: a cell array of dimension M, where each index 1:M
%     stores a N x 2 array of the coordinates of the nested shape inside,
%     with M = 1 being the starting shape (and with dimension K+1, where K
%     is number of vertices given), and N being smaller and smaller for
%     each M value.
%
%     projection_vectors: a cell array of M, where each index 1:M
%     stores a N x 2 array of the unit vectors that point
%     away from the vertices of the nested shape inside, with M = 1 being
%     the starting unit vectors and N being smaller and smaller for each M value.
%
%     cut_distance: an array of 1 x M, starting from 0 for M(1) to the
%     maximum cut distance that can be used, at M(end)
%
%    (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%     h_fig: a handle to the resulting figure
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% % EXAMPLES:
% For additional examples, see: script_test_fcn_MapGen_polytopePlotSkeleton
%
% This function was written on 2022_02_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%

% Revision History:
% 2022_02_16 - S. Brennan
% -- first write of code
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
if (nargin==4 && isequal(varargin{end},-1))
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
    if flag_check_inputs
        % Are there the right number of inputs?
        if nargin < 3 || nargin > 4
            error('Incorrect number of input arguments')
        end

        % Check the cut_distance input
        fcn_DebugTools_checkInputsToFunctions(...
            cut_distance, '1column_of_numbers');

    end
end


% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  (4 == nargin) && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_fig = figure(fig_num);
grid on
grid minor
hold on
axis equal

% Grab the original vertices
original_vertices = vertices{1};

% Plot the original polytope in red using the vertices
h_plot = plot(original_vertices(:,1),original_vertices(:,2),'r-','Linewidth',2);
last_color = get(h_plot,'Color');

% Find size of vertices so we can figure out nudging
size = max(max(original_vertices)) - min(min(original_vertices));
nudge = size*0.003;

% Number the original vertices with labels, nudging a litte so they don't
% land right on top of the points
for ith_vertex = 1:(length(original_vertices(:,1))-1)
    text(original_vertices(ith_vertex,1)+nudge,original_vertices(ith_vertex,2),...
        sprintf('%.0d',ith_vertex));
end

% Plot each contraction
Ncontractions = length(cut_distance);
for ith_contraction = 2:Ncontractions

    % Do calculations to determine vector start points, length of vectors,
    % and vectors themselves
    cut_length = cut_distance(ith_contraction)-cut_distance(ith_contraction-1);
    starting_points = vertices{ith_contraction-1}(1:end-1,:);
    vectors_from_starting_points = ...
        projection_vectors{ith_contraction-1}(1:end-1,:)*cut_length;

    % Plot the vectors as arrows going out from the starting point to
    % ending point. Use the color from the previous plot of the points, so
    % the vectors are same color as the points they originate from.
    quiver(starting_points(:,1),starting_points(:,2),vectors_from_starting_points(:,1),vectors_from_starting_points(:,2),0,'Color',last_color);

    % Plot the ending points, and save the color since they will be the new
    % start points on the next iteration.
    ending_points = vertices{ith_contraction}(1:end-1,:);
    h_plot = plot(ending_points(:,1),ending_points(:,2),'.','Markersize',20);
    last_color = get(h_plot,'Color');

end

%% Plot results?
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
    % Nothing to do here - it's a plotting function
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends fcn_MapGen_polytopePlotSkeleton
