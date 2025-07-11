function [croppedVertices] = ...
    fcn_MapGen_cropPolytopeToRange(...
    vertices, ...
    interiorPoint,...
    AABB,...
    varargin...
    )
% fcn_MapGen_cropPolytopeToRange
% crops the given vertices of a polytope to an axis-aligned bounding box
%
% FORMAT:
%
%    [croppedVertices] = ...
%     fcn_MapGen_cropPolytopeToRange(...
%     vertices, ...
%     interiorPoint,...
%     AABB,...
%    (fig_num) ...
%    )
%
% INPUTS:
%
%     vertices: the list of vertex points defining the polytope, in [x y]
%     format, where x and y are columns
%
%     interiorPoint: a point inside the polytope, in [x y]
%     format, where x and y are scalars
%
%     AABB: the axis-aligned bounding box, in format of
%     [xmin ymin xmax ymax], wherein the resulting polytopes must be
%     bounded.
%
%     (optional inputs)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
%
% OUTPUTS:
%
%     croppedVertices: the resulting vertices of the polytope, forced to
%     fit within the AABB
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_convertAABBtoWalls
%     fcn_MapGen_polytopeProjectVerticesOntoWalls
%     fcn_MapGen_polytopeRemoveColinearVertices
%     fcn_Path_findSensorHitOnWall
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_cropPolytopeToRange
% for a full test suite.
%
% This function was written on 2021_07_11 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

%
% REVISION HISTORY:
%
% 2021_07_11 by Sean Brennan
% -- first write of function
% 2021_07_30 by Sean Brennan
% -- fixed errors due to corners being omitted
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_07_10 by Sean Brennan
% -- changed fcn_MapGen_findIntersectionOfSegments to use
% fcn_Path_findSensorHitOnWall instead, as the Path function is much more
% tested/debugged and regularly updated
% -- renamed variables for clarity
% -- improved plotting

% TO DO:
%
% -- allow user to enter the allowable range (hard-coded now to 0 to 1)
% -- check that interior point is actually inside vertices!

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
    flag_do_debug = 0; %     % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; %     % Flag to plot the results for debugging
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

if 0==flag_max_speed
    if 1 == flag_check_inputs

        % Are there the right number of inputs?
        if nargin < 3 || nargin > 4
            error('Incorrect number of input arguments')
        end

        % Check the vertices input, make sure it is '2column_of_numbers' type
        fcn_DebugTools_checkInputsToFunctions(...
            vertices, '2column_of_numbers');

        % Check the interiorPoint input, make sure it is '2column_of_numbers'
        % type, with 1 row
        fcn_DebugTools_checkInputsToFunctions(...
            interiorPoint, '2column_of_numbers',1);

        % Check the AABB input, make sure it is '4column_of_numbers' type, with
        % 1 row
        fcn_DebugTools_checkInputsToFunctions(...
            AABB, '4column_of_numbers',1);

    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  4 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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

% % For debugging
% tolerance = 0.001;
% location = [0.6128 0.9867];
% if (...
%         (interiorPoint(1,1)<location(1)+tolerance) && ...
%         (interiorPoint(1,1)>location(1)-tolerance) && ...
%         (interiorPoint(1,2)<location(2)+tolerance) && ...
%         (interiorPoint(1,2)>location(2)-tolerance))
%     disp('stop here');
% end

% Convert axis-aligned bounding box to wall format
walls = fcn_MapGen_convertAABBtoWalls(AABB, -1);

% Open the figure if doing debugging
if flag_do_debug
    figure(1);
    clf;
    hold on;

    % Set the axis
    scale = max(AABB,[],'all') - min(AABB,[],'all');
    new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];
    axis(new_axis);

    % Plot the original vertices
    plot(...
        [vertices(:,1); vertices(1,1)],...
        [vertices(:,2); vertices(1,2)],...
        '.-','Linewidth',3);

    % Plot the walls
    plot(walls(:,1),walls(:,2),'k-');

    % Plot the interior point
    plot(interiorPoint(:,1),interiorPoint(:,2),'ro');

end

% Nudge the interior point inward, if it is on a border
interiorPoint = ...
    fcn_INTERNAL_nudgeInteriorPointInward(interiorPoint, AABB);

if flag_do_debug
    % Plot the new interior point
    figure(1);
    plot(interiorPoint(:,1),interiorPoint(:,2),'ro');

end

vertices_no_infinite = vertices;

% Sometimes the polytopes intersect the box boundaries. We can artificially
% add these border crossings as extra points so that we can project the
% polytope correctly back onto walls (in a later step).
[all_points, flag_was_intersection] = ...
    fcn_INTERNAL_findAllPoints(vertices_no_infinite, walls);

if flag_do_debug
    figure(1);
    % Plot the all_points locations
    plot(all_points(:,1),all_points(:,2),'kx');
end

% Check for the enclosing case where the polytope goes completely around
% the bounding box (e.g. bounding box is INSIDE the polytope!?!). In this
% case, there will be no projection, and so we should just exit.
flag_vertices_outside = ((vertices_no_infinite(:,1)>=AABB(1,3)) + ...
    (vertices_no_infinite(:,1)<=AABB(1,1))).*((vertices_no_infinite(:,2)>=AABB(1,4)) + ...
    (vertices_no_infinite(:,2)<=AABB(1,2)));
if all(flag_vertices_outside) && (flag_was_intersection==0)
    croppedVertices = walls;
else

    % From the interior point, project all_points back onto the wall to create
    % a polytope limited by the bounding box.    
    [projected_points] = ...
        fcn_MapGen_polytopeProjectVerticesOntoWalls(...,
        interiorPoint,...
        all_points,...
        walls(1:end-1,:),...
        walls(2:end,:), -1);

    if flag_do_debug
        figure(1);
        % Plot the projected_points locations
        plot(projected_points(:,1),projected_points(:,2),'go-');
    end

    % Use the cross-product to eliminate co-linear points, as sometimes the
    % above process generates multiple points in a line, which is technically
    % not a polytope.
    [croppedVertices] = ...
        fcn_MapGen_polytopeRemoveColinearVertices(...,
        projected_points, -1);

    % Sometimes the cross-product step above removes the repeated last vertex.
    % So we may have to fix this
    if ~isempty(croppedVertices)
        if ~isequal(croppedVertices(1,:),croppedVertices(end,:))
            croppedVertices = [croppedVertices; croppedVertices(1,:)];
        end
    else
        error('Error in cropPolytopeToRange');
    end



end % Ends if statement to check if all points are enclosed


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
    hold on;

    % Plot the input interiorPoint
    plot(interiorPoint(:,1),interiorPoint(:,2),'g.','MarkerSize', 30, 'DisplayName','Input: interiorPoint');

    % Plot the original vertices
    plot(...
        [vertices(:,1); vertices(1,1)],...
        [vertices(:,2); vertices(1,2)],...
        '.-', 'DisplayName','Input: vertices');

    % Plot the AABB
    plot(walls(:,1),walls(:,2),'k.-','DisplayName','Input: AABB');

    % Plot the croppedVertices locations
    plot(croppedVertices(:,1),croppedVertices(:,2),'mo-','DisplayName','Output: croppedVertices');

    legend('Interpreter','none');

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

%% fcn_INTERNAL_findAllPoints
function [all_points, flag_was_intersection] = fcn_INTERNAL_findAllPoints(vertices,walls)
% Finds all vertices where the polytope intersects the walls, as well as
% any wall start locations that are within the polytope

% Fill in blank all_points as starter
all_points = [];

% Check if any of the corners, e.g. where the walls start, are inside the
% polytope. If so, need to add these.
test_points = walls(1:4,:);
in_polytope = inpolygon(test_points(:,1),test_points(:,2),vertices(:,1),vertices(:,2));
if any(in_polytope)
    all_points = [all_points; test_points(in_polytope,:)];
end

% Pad the vertices to wrap around, so we don't miss the last wall
vertices = [vertices; vertices(1,:)];
start_vertices = vertices(1:end-1,:);
end_vertices = vertices(2:end,:);



flag_was_intersection = 0;
% Find all intersection points
for ith_point = 1:length(start_vertices(:,1))
    all_points = [all_points; start_vertices(ith_point,:)]; %#ok<AGROW>

    % Define start and end of the sensor
    sensor_vector_start = start_vertices(ith_point,:);
    sensor_vector_end = end_vertices(ith_point,:);

    % Define the start and end of the walls
    wall_start = walls(1:end-1,:);
    wall_end = walls(2:end,:);

    % Call a function that determines where and if the sensor crosses the
    % walls
    [distance, location] = ...
        fcn_Path_findSensorHitOnWall(...
        wall_start,...           % wall start
        wall_end,...             % wall end
        sensor_vector_start,...  % sensor_vector_start
        sensor_vector_end,...    % sensor_vector_end
        (1), ...                 % (flag_search_return_type) -- 1 means ALL hits of any results,
        (0), ...                 % (flag_search_range_type)  -- 0 means only if overlapping wall/sensor, ...
        ([]),...                 % (tolerance) -- default is eps * 1000,
        (-1));                   % (fig_num) -- -1 means to use "fast mode")

    if ~isnan(distance)
        all_points = [all_points; location]; %#ok<AGROW>
        flag_was_intersection = 1;
    end

end

% Get rid of duplicates (occurs when two points are both on edges)
indices_not_repeated = [~all(abs(diff(all_points))<eps*10,2); 1];
all_points = all_points(indices_not_repeated>0,:);

end % Ends fcn_INTERNAL_findAllPoints

%% fcn_INTERNAL_nudgeInteriorPointInward
function interiorPoint = fcn_INTERNAL_nudgeInteriorPointInward(interiorPoint,box)
% If the interior point is on the edge of the box, or even outside, this
% nudges the interior point to the true interior.

nudge = 1e-8;

if interiorPoint(1,1)<=box(1,1)
    interiorPoint(1,1) = nudge;
elseif interiorPoint(1,1)>=box(1,3)
    interiorPoint(1,1) = 1 - nudge;
end
if interiorPoint(1,2)<=box(1,2)
    interiorPoint(1,2) = nudge;
elseif interiorPoint(1,2)>=box(1,4)
    interiorPoint(1,2) = 1 - nudge;
end
end % Ends fcn_INTERNAL_nudgeInteriorPointInward



