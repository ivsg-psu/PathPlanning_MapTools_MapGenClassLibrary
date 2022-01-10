function [cropped_vertices] = ...
    fcn_MapGen_cropPolytopeToRange(...
    vertices, ...
    interior_point,...
    AABB,...
    varargin...
    )
% fcn_MapGen_cropPolytopeToRange
% crops the given vertices of a polytope to an axis-aligned bounding box
%
% FORMAT:
%
%    [cropped_vertices] = ...
%     fcn_MapGen_cropPolytopeToRange(...
%     vertices, ...
%     interior_point,...
%     AABB,...
%    (fig_num) ...
%    )
%
% INPUTS:
%
%     vertices: the list of vertex points defining the polytope, in [x y]
%     format, where x and y are columns
%
%     interior_point: a point inside the polytope, in [x y]
%     format, where x and y are scalars
%
%     AABB: the axis-aligned bounding box, in format of
%     [xmin ymin xmax ymax], wherein the resulting polytopes must be
%     bounded.
%
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a
%     figure to be drawn showing results.
%
%
% OUTPUTS:
%
%     cropped_vertices: the resulting vertices of the polytope, forced to
%     fit within the AABB
%
%
% DEPENDENCIES:
%
%     fcn_MapGen_checkInputsToFunctions
%     fcn_MapGen_convertAABBtoWalls
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

%
% TO DO:
%
% -- allow user to enter the allowable range (hard-coded now to 0 to 1)
% -- check that inrerior point is inside vertices

%% Debugging and Input checks
% set an environment variable on your machine with the getenv function...
% in the Matlab console.  Char array of '1' will be true and '0' will be false.
flag_check_inputs = getenv('ENV_FLAG_CHECK_INPUTS');  % '1' will check input args
flag_do_plot = getenv('ENV_FLAG_DO_PLOT'); % '1' will make plots
flag_do_debug = getenv('ENV_FLAG_DO_DEBUG'); % '1' will enable debugging

% if the char array has length 0, assume the env var isn't set and default to...
% dipslaying more information rather than potentially hiding an issue
if length(flag_check_inputs) == 0
    flag_check_inputs = '1';
end
if length(flag_do_plot) == 0
    flag_do_plot = '1';
end
if length(flag_do_debug) == 0
    flag_do_debug = '1';
end

% convert flag from char string to logical
flag_check_inputs = flag_check_inputs == '1';
flag_do_plot = flag_do_plot == '1';
flag_do_debug = flag_do_debug == '1';

if flag_do_debug
    fig_for_debug = 27564;
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
    if nargin < 3 || nargin > 4
        error('Incorrect number of input arguments')
    end

    % Check the vertices input, make sure it is '2column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        vertices, '2column_of_numbers');

    % Check the interior_point input, make sure it is '2column_of_numbers'
    % type, with 1 row
    fcn_MapGen_checkInputsToFunctions(...
        interior_point, '2column_of_numbers',1);

    % Check the AABB input, make sure it is '4column_of_numbers' type, with
    % 1 row
    fcn_MapGen_checkInputsToFunctions(...
        AABB, '4column_of_numbers',1);

end

% Does user want to show the plots?
if  4== nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig_num = fig_for_debug;
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
%         (interior_point(1,1)<location(1)+tolerance) && ...
%         (interior_point(1,1)>location(1)-tolerance) && ...
%         (interior_point(1,2)<location(2)+tolerance) && ...
%         (interior_point(1,2)>location(2)-tolerance))
%     disp('stop here');
% end

% Convert axis-aligned bounding box to wall format
walls = fcn_MapGen_convertAABBtoWalls(AABB);

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
    plot(interior_point(:,1),interior_point(:,2),'ro');

end

% Nudge the interior point inward, if it is on a border
interior_point = ...
    INTERNAL_fcn_nudgeInteriorPointInward(interior_point,AABB);

if flag_do_debug
    % Plot the new interior point
    figure(1);
    plot(interior_point(:,1),interior_point(:,2),'ro');

end

vertices_no_infinite = vertices;

% Sometimes the polytopes intersect the box boundaries. We can artificially
% add these border crossings as extra points so that we can project the
% polytope correctly back onto walls (in a later step).
[all_points, flag_was_intersection] = ...
    INTERNAL_fcn_findAllPoints(vertices_no_infinite,walls);

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
    cropped_vertices = walls;
else

    % From the interior point, project all_points back onto the wall to create
    % a polytope limited by the bounding box.
    [projected_points] = ...
    fcn_MapGen_polytopeProjectVerticesOntoWalls(...,
    interior_point,...
    all_points,...
    walls(1:end-1,:),...
    walls(2:end,:));

    if flag_do_debug
        figure(1);
        % Plot the projected_points locations
        plot(projected_points(:,1),projected_points(:,2),'go-');
    end

    % Use the cross-product to eliminate co-linear points, as sometimes the
    % above process generates multiple points in a line, which is technically
    % not a polytope.
    [cropped_vertices] = ...
    fcn_MapGen_polytopeRemoveColinearVertices(...,
    projected_points);

    % Sometimes the cross-product step above removes the repeated last vertex.
    % So we may have to fix this
    if ~isempty(cropped_vertices)
        if ~isequal(cropped_vertices(1,:),cropped_vertices(end,:))
            cropped_vertices = [cropped_vertices; cropped_vertices(1,:)];
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
    clf;
    hold on;
    grid on;
    grid minor;

    % Plot the original vertices
    plot(...
        [vertices(:,1); vertices(1,1)],...
        [vertices(:,2); vertices(1,2)],...
        '.-');

    % Plot the walls
    plot(walls(:,1),walls(:,2),'k-');

    % Plot the interior point
    plot(interior_point(:,1),interior_point(:,2),'ro');

    % Plot the cropped_vertices locations
    plot(cropped_vertices(:,1),cropped_vertices(:,2),'mo-');

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


function [all_points, flag_was_intersection] = INTERNAL_fcn_findAllPoints(vertices,walls)
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
    [distance,location,~] = ...
        fcn_MapGen_findIntersectionOfSegments(...
        wall_start,...
        wall_end,...
        sensor_vector_start,...
        sensor_vector_end,2);

    if ~isnan(distance)
        all_points = [all_points; location]; %#ok<AGROW>
        flag_was_intersection = 1;
    end

end

% Get rid of duplicates (occurs when two points are both on edges)
indices_not_repeated = [~all(abs(diff(all_points))<eps*10,2); 1];
all_points = all_points(indices_not_repeated>0,:);

end % Ends INTERNAL_fcn_findAllPoints


function interior_point = INTERNAL_fcn_nudgeInteriorPointInward(interior_point,box)
% If the interior point is on the edge of the box, or even outside, this
% nudges the interior point to the true interior.

nudge = 1e-8;

if interior_point(1,1)<=box(1,1)
    interior_point(1,1) = nudge;
elseif interior_point(1,1)>=box(1,3)
    interior_point(1,1) = 1 - nudge;
end
if interior_point(1,2)<=box(1,2)
    interior_point(1,2) = nudge;
elseif interior_point(1,2)>=box(1,4)
    interior_point(1,2) = 1 - nudge;
end
end % Ends INTERNAL_fcn_nudgeInteriorPointInward



