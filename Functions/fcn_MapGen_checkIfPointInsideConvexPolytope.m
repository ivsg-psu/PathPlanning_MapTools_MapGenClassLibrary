function [ in_polytope ] = ...
    fcn_MapGen_checkIfPointInsideConvexPolytope( ...
    test_point, ...
    vertices, ...
    varargin...
    )
% fcn_MapGen_checkIfPointInsideConvexPolytope
% Given an convex polytope, checks if point is within the polytope
%
%
%
% FORMAT:
%
%    [ ...
%    snap_point ...
%    ] = ...
%    fcn_MapGen_checkIfPointInsideConvexPolytope( ...
%    test_point, ...
%    vertices, ...
%    (fig_num) ...
%    )
%
% INPUTS:
%
%     test_point: the test point, in format [x y]
%
%     vertices: the list of vertex points defining the polytope, in [x y]
%     format, where x and y are columns
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
%     in_polytope: 1 if the point is within the polytope, 0 otherwise
%
%
% DEPENDENCIES:
%
%     (none)
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_checkIfPointInsideConvexPolytope
% for a full test suite.
%
% This function was written on 2021_07_14 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

%
% REVISION HISTORY:
%
% 2021_07_14 by Sean Brennan
% -- first write of function
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
if (nargin==3 && isequal(varargin{end},-1))
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
        if nargin < 2 || nargin > 3
            error('Incorrect number of input arguments')
        end

        % Check the test_point input, make sure it is '2column_of_numbers' type
        fcn_DebugTools_checkInputsToFunctions(...
            test_point, '2column_of_numbers',1);

        % Check the vertices input, make sure it is '2column_of_numbers' type
        fcn_DebugTools_checkInputsToFunctions(...
            vertices, '2column_of_numbers');
    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  3 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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
mean_point = mean(vertices(vertices(1:end-1,1)~=inf,:),1);

[distance,~,~] = ...
    INTERNAL_fcn_geometry_findIntersectionOfSegments(...
    vertices(1:end-1,:), ... % wall_start,...
    vertices(2:end,:), ... % wall_end
    test_point,... % sensor_vector_start,...
    mean_point, ... % sensor_vector_end,...
    2); % Return all hits

% If the above test passes on or through a vertex, then there will be
% multiple hits with the same distance value. We eliminate them by keeping
% only the unique values.
unique_distances = unique(distance);
try
    if isempty(unique_distances) || all(unique_distances==0)
        in_polytope = true;
    else
        in_polytope = false;
    end
catch
    disp('Stop here');
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
    % clf;
    hold on;
    axis equal
    grid on;
    
    % Plot the test point
    plot(test_point(:,1),test_point(:,2),'x');
    
    % Plot the vertices
    plot(vertices(:,1),vertices(:,2),'b-');
    
    % Plot the mean point
    plot(test_point(:,1),test_point(:,2),'mx');
    
    % Plot the result
    if in_polytope
        plot(test_point(:,1),test_point(:,2),'go');
    end
    
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



function [distance,location,wall_that_was_hit] = ...
    INTERNAL_fcn_geometry_findIntersectionOfSegments(...
    wall_start,...
    wall_end,...
    sensor_vector_start,...
    sensor_vector_end,...
    varargin)
% fcn_geometry_findIntersectionOfSegments calculates hits between a sensor
% projection and a set of walls, both specified by start and end points,
% returning the distance, and location of the hit, and the wall number of
% each hit.
%
% FORMAT:
%
%      [distance,location,wall_that_was_hit] = ...
%         fcn_geometry_findIntersectionOfSegments(...
%         wall_start,...
%         wall_end,...
%         sensor_vector_start,...
%         sensor_vector_end,...
%         (flag_search_type),(fig_num))
%
% INPUTS:
%
%      wall_start: an N x 2 vector containing the X,Y points of the
%      starting points of each "wall".
%
%      wall_end: an N x 2 vector containing the X,Y points of the
%      ending points of each "wall".
%
%      sensor_vector_start: a 1 x 2 vector containing the X,Y points of the
%      sensor's start location
%
%      sensor_vector_end: a 1 x 2 vector containing the X,Y points of the
%      sensor's end location
%
%      (OPTIONAL INPUTS)
%      flag_search_type: an integer specifying the type of search.
%
%            0: return distance and location of the first point of overlap,
%            only if the given sensor_vector overlaps the wall (this is the
%            default)
%
%            1: return distane and location if any projection of the sensor
%            vector, in any direction, hits the wall (in other words, if
%            there is any intersection). Note that distance returned will
%            be negative if the nearest intersection is in the opposite
%            direction of the given sensor vector.
%
%            2: returns distances and locations of all hits (not just the
%            first one as in the options above) as M x 1 and M x 2 vectors
%            respectively, where the M rows represent all the detected
%            intersections.
%
%      fig_num: a figure number to plot results. Turns debugging on.
%
% OUTPUTS:
%
%      distance: a N x 1 scalar representing the distance to the closest
%      intersection of the sensor with the path. NaN is returned if not
%      detected.
%
%      location: a N x 2 vector of the X,Y location of intersection point
%
%      wall_that_was_hit: the segment number of the wall that was hit (1 is
%      the first segment, 2 is the second, etc)
%
% EXAMPLES:
%
%       See the script: script_test_fcn_geometry_findIntersectionOfSegments.m
%       for a full test suite.
%
% Adopted from https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
% This function was written on 2021_06_05 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2021_06_05
%      - wrote the code, templated from fcn_geometry_findIntersectionOfSegments



%% Set up for debugging
flag_do_debug = 0; % Flag to plot the results for debugging
flag_do_plot = 0; % Flag to plot the results for debugging
flag_check_inputs = 0; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'Starting function: %s, in file: %s\n',st(1).name,st(1).file);
end

%% check input arguments
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
% Set default values
flag_search_type = 0;

% check input arguments
if flag_check_inputs == 1
    if nargin < 4 || nargin > 6
        error('Incorrect number of input arguments.')
    end
    
    % Check wall_start input
    fcn_geometry_checkInputsToFunctions(wall_start, '2column_of_numbers');
    
    Nwalls = length(wall_start(:,1));
    
    % Check wall_end input
    fcn_geometry_checkInputsToFunctions(wall_end, '2column_of_numbers',Nwalls);
    
    % Check sensor_vector_start input
    fcn_geometry_checkInputsToFunctions(sensor_vector_start, '2column_of_numbers',1);
    
    % Check sensor_vector_end input
    fcn_geometry_checkInputsToFunctions(sensor_vector_end, '2column_of_numbers',1);
end

% Does user wish to specify search type?
if 5 <= nargin
    flag_search_type = varargin{1};
end


% Does user want to show the plots?
if 6 == nargin
    fig_num = varargin{end};
    figure(fig_num);
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plot = 1;
    end
end


%% Calculations begin here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wall_numbers = (1:length(wall_start(:,1)))';

% Define p, q, r and s vectors
p = wall_start;
q = sensor_vector_start;
r = wall_end - wall_start;
s = sensor_vector_end - sensor_vector_start;

% Define useful intermediate terms
r_cross_s = INTERNAL_crossProduct(r,s);
q_minus_p =  q - p;
q_minus_p_cross_s = INTERNAL_crossProduct(q_minus_p,s);
q_minus_p_cross_r = INTERNAL_crossProduct(q_minus_p,r);

% Are any of these parallel?
parallel_non_intersecting_indices = find((0==r_cross_s).*(0~=q_minus_p_cross_r));
if any(parallel_non_intersecting_indices)
    r_cross_s(parallel_non_intersecting_indices) = 1; % They are colinear or parallel, so make dummy length
end

% Are any of these colinear?
colinear_indices = find((0==r_cross_s).*(0==q_minus_p_cross_r));
if any(colinear_indices)
    r_cross_s(colinear_indices) = 1; % They are colinear or parallel, so make dummy length
    r_dot_r = sum(r.*r,2);
    q_minus_p_dot_r = sum(q_minus_p.*r,2);
    s_dot_r = sum(s.*r,2);
    t0 = q_minus_p_dot_r./r_dot_r;
    t1 = t0 + s_dot_r./r_dot_r;
    
    % Keep only the good indices
    %     % For debugging:
    %     conditions = [-0.5 -0.4; -0.5 0; -0.5 .2; 0 0.2; 0.2 0.4; 0.2 1; 0.2 1.2; 1 1.2; 1.2 1.3; -0.5 1.2]
    %     t0 = conditions(:,1);
    %     t1 = conditions(:,2);
    %     t0_inside = (t0>=0)&(t0<=1);
    %     t1_inside = (t1>=0)&(t1<=1);
    %     t0_t1_surround = (t0<0)&(t1>1) | (t1<0)&(t0>1);
    %     any_within = t0_inside | t1_inside | t0_t1_surround;
    %     fprintf(1,'t0 inside flag:\n');
    %     [conditions t0_inside]
    %     fprintf(1,'t1 inside flag:\n');
    %     [conditions t1_inside]
    %     fprintf(1,'surround flag:\n');
    %     [conditions t0_t1_surround]
    %     fprintf(1,'any_wighin flag:\n');
    %     [conditions any_within]
    
    % Check whether there is overlap by seeing of the t0 and t1 values are
    % within the interval of [0 1], endpoint inclusive
    t0_inside = (t0>=0)&(t0<=1);
    t1_inside = (t1>=0)&(t1<=1);
    t0_t1_surround = (t0<0)&(t1>1) | (t1<0)&(t0>1);
    any_within = t0_inside | t1_inside | t0_t1_surround;
    good_indices = find(any_within);
    good_colinear_indices = intersect(colinear_indices,good_indices);
    
    % Fix the ranges to be within 0 and 1
    t0(good_colinear_indices) = max(0,t0(good_colinear_indices));
    t0(good_colinear_indices) = min(1,t0(good_colinear_indices));
    
    t1(good_colinear_indices) = max(0,t1(good_colinear_indices));
    t1(good_colinear_indices) = min(1,t1(good_colinear_indices));
    
end

% Calculate t and u, where t is distance along the path, and u is distance
% along the sensor.

t = q_minus_p_cross_s./r_cross_s; % Distance along the path
u = q_minus_p_cross_r./r_cross_s; % Distance along the sensor

% Fix any situations that are parallel and non-intersecting, as these will
% give wrong calculation results from the t and u calculations above
t(parallel_non_intersecting_indices) = inf;
u(parallel_non_intersecting_indices) = inf;

% Fix any situations that are colinear. For these, we save the start point
% as the point where overlap starts, and the end point where overlap ends.
% Note that this creates NEW intersections beyond the number of segements
if any(colinear_indices)
    % Shut off colinear ones to start
    t(colinear_indices) = inf;
    u(colinear_indices) = inf;
    
    % Correct the t values
    u(good_colinear_indices) = 1;
    t(good_colinear_indices) = t0(good_colinear_indices);
    
    % Do we need to add more hit points?
    indices_hit_different_point = find(t0~=t1);
    more_indices = intersect(indices_hit_different_point,good_colinear_indices);
    
    % Make p and r, t and u longer so that additional hit points are
    % calculated in special case of overlaps
    p = [p; p(more_indices,:)];
    r = [r; r(more_indices,:)];
    u = [u; u(more_indices)];
    t = [t; t1(more_indices)];
    wall_numbers = [wall_numbers; wall_numbers(more_indices)];
    
end

% Initialize all intersections to infinity
intersections = NaN*ones(length(p(:,1)),2);

% Note: could speed this up with nested if logical statements, but only if
% are doing checks on single segments at a time. Since doing many segments
% at once, need to use vector form.
if 0 == flag_search_type % Overlap of wall and vector, keep only first hit
    good_vector = ((0<=t).*(1>=t).*(0<=u).*(1>=u));
elseif 1 == flag_search_type % Overlap of wall, any direction of vector, keep first hit
    good_vector = ((0<=t).*(1>=t));
elseif 2 == flag_search_type % Same as zero, but finds all the hits
    good_vector = ((0<=t).*(1>=t).*(0<=u).*(1>=u));
elseif 3 == flag_search_type % First wall in direction of vector, keep first hit
    good_vector = ((0<=t).*(1>=t).*(0<=u));
else
    error('Incorrect flag_search_type entered');
end

% Keep only the indices that work
good_indices = find(good_vector>0);

if ~isempty(good_indices)
    result = p + t.*r;
    intersections(good_indices,:) = result(good_indices,:);
end

% Find the distances via Euclidian distance to the sensor's origin
% note: a faster way to do this might be to just
% calculate t*r as a length
distances_squared = sum((intersections - sensor_vector_start).^2,2);

if flag_search_type ~=2
    % Keep just the minimum distance
    [closest_distance_squared,closest_index] = min(distances_squared);
    
    distance = closest_distance_squared^0.5*sign(u(closest_index));
    location = intersections(closest_index,:);
    wall_that_was_hit = wall_numbers(closest_index);
else
    % Return all the results, not just minimum
    good_indices = find(~isnan(distances_squared));
    distance = distances_squared(good_indices).^0.5.*sign(u(good_indices));
    location = intersections(good_indices,:);
    wall_that_was_hit = wall_numbers(good_indices);
end

%% Any debugging?
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
    
    % Set up the figure
    figure(fig_num);
    clf;
    hold on;
    axis equal;
    grid on; grid minor;
    
    % Calculate the midpoints
    midpoints= (wall_start+wall_end)/2;
    
    % Plot the walls
    for ith_wall = 1:Nwalls
        handle_plot = plot(...
            [wall_start(ith_wall,1) wall_end(ith_wall,1)],...
            [wall_start(ith_wall,2) wall_end(ith_wall,2)],...
            '.-','Linewidth',5);
        line_color = get(handle_plot,'Color');
        handle_text = text(midpoints(ith_wall,1),midpoints(ith_wall,2),sprintf('%.0d',wall_numbers(ith_wall,1)));
        set(handle_text,'Color',line_color);
        set(handle_text,'Fontsize',20);
    end
    
    % Plot the sensor vector
    quiver(q(:,1),q(:,2),s(:,1),s(:,2),'r','Linewidth',3);
    plot(sensor_vector_end(:,1),sensor_vector_end(:,2),'r.','Markersize',10);
    
    handle_text = text(q(:,1),q(:,2),'Sensor');
    set(handle_text,'Color',[1 0 0]);
    
    axis_size = axis;
    y_range = axis_size(4)-axis_size(3);
    
    % Plot any hits in blue
    for i_result = 1:length(distance)
        plot(location(i_result,1),location(i_result,2),'bo','Markersize',30);
        handle_text = text(location(i_result,1),location(i_result,2)-0.05*y_range,sprintf('Hit %.0d at distance: %.2f',wall_that_was_hit(i_result),distance(i_result)));
        set(handle_text,'Color',[0 0 1]);
    end
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end
end



%% Calculate cross products
function result = INTERNAL_crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end

