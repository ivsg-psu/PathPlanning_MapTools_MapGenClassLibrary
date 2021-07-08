function [cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies, interior_point)

flag_do_debug = 0;

% TO-DO:
% -- allow user to enter the allowable range (hard-coded now to 0 to 1)
% -- check that inrerior point is inside vertices

% snap prior to closest wall
box = [0 0 1 1]; % [xmin ymin xmax ymax]

% Convert axis-aligned bounding box to wall format
walls = [box(1,1) box(1,2); box(1,3) box(1,2); box(1,3) box(1,4); box(1,1) box(1,4); box(1,1) box(1,2)];

% % FOR DEBUGGING
% if flag_do_debug
%     if(interior_point(1,1)<0.002)&&(interior_point(1,2)>0.965)
%         disp('Stop here');
%     end
% end

% Open the figure if doing debugging
if flag_do_debug
    figure(1);
    clf;
    hold on;
    
    % Plot the vertices
    plot(...
        [verticies(:,1); verticies(1,1)],...
        [verticies(:,2); verticies(1,2)],...
        '.-');
    
    % Plot the walls
    plot(walls(:,1),walls(:,2),'k-');
    
    % Plot the interior point
    plot(interior_point(:,1),interior_point(:,2),'ro');

end

% Are any verticies infinite? If so, we need to check that the adjacent
% verticies will create a reasonable polytope. 
verticies_no_infinite = INTERNAL_fcn_removeInfiniteVerticies(verticies,box);

% Open the figure if doing debugging
if flag_do_debug
    % Plot the vertices
    plot(...
        [verticies_no_infinite(:,1); verticies_no_infinite(1,1)],...
        [verticies_no_infinite(:,2); verticies_no_infinite(1,2)],...
        '.-');
end

% Nudge the interior point inward, if it is on a border
interior_point = INTERNAL_fcn_nudgeInteriorPointInward(interior_point,box);

if flag_do_debug
    % Plot the new interior point
    plot(interior_point(:,1),interior_point(:,2),'ro');

end

% Sometimes the polytopes intersect the box boundaries. We can artificially
% add these border crossings as extra points so that we can project the
% polytope correctly back onto walls (in a later step).
[all_points, flag_was_intersection] = INTERNAL_fcn_findAllPoints(verticies_no_infinite,walls);

if flag_do_debug
    % Plot the all_points locations
    plot(all_points(:,1),all_points(:,2),'kx');
end

% Check for the enclosing case where the polytope goes completely around
% the bounding box (e.g. bounding box is INSIDE the polytope!?!). In this
% case, there will be no projection, and so we should just exit.
flag_vertices_outside = ((verticies_no_infinite(:,1)>=1) + (verticies_no_infinite(:,1)<=0)).*((verticies_no_infinite(:,2)>=1) + (verticies_no_infinite(:,2)<=0));
if all(flag_vertices_outside) && (flag_was_intersection==0)
    cropped_vertices = walls;
    return;
end

% From the interior point, project all_points back onto the wall to create
% a polytope limited by the bounding box.
projected_points = INTERNAL_fcn_projectAllPointsOntoWalls(interior_point, all_points,walls);

if flag_do_debug
    % Plot the projected_points locations
    plot(projected_points(:,1),projected_points(:,2),'go-');
end

% Use the cross-product to eliminate co-linear points, as sometimes the
% above process generates multiple points in a line, which is technically
% not a polytope.
cropped_vertices = INTERNAL_fcn_cropRepeatedPoints(projected_points);

% Sometimes the cross-product step above removes the repeated last vertex.
% So we may have to fix this
if ~isempty(cropped_vertices)
    if ~isequal(cropped_vertices(1,:),cropped_vertices(end,:))
        cropped_vertices = [cropped_vertices; cropped_vertices(1,:)];
    end
else
    disp('Stop here');
end

if flag_do_debug
    % Plot the projected_points locations
    plot(cropped_vertices(:,1),cropped_vertices(:,2),'mo-');
end

end

function cropped_vertices = INTERNAL_fcn_cropRepeatedPoints(projected_points)
% Use the cross-product to eliminate co-linear points

Npoints = length(projected_points(:,1));
good_indices = zeros(Npoints,1);

for ith_point = 1:Npoints
    
    if ith_point>1
        previous_vector = [projected_points(ith_point,:)-projected_points(ith_point-1,:), 0];
    else % must wrap around
        previous_vector = [projected_points(end,:)-projected_points(1,:), 0];
    end
    
    if ith_point<Npoints
        subsequent_vector = [projected_points(ith_point+1,:)-projected_points(ith_point,:), 0];
    else % must wrap around
        subsequent_vector = [projected_points(1,:)-projected_points(Npoints,:), 0];
    end
    
    cross_result = cross(previous_vector,subsequent_vector);
    
    if cross_result(1,3)~=0
        good_indices(ith_point,:) = 1;
    end
end

% Final polytope
cropped_vertices = projected_points(good_indices>0,:);

end % Ends INTERNAL_fcn_cropRepeatedPoints



function projected_points = INTERNAL_fcn_projectAllPointsOntoWalls(interior_point, all_points,walls)
% From the interior point, project all_points back onto the wall.

projected_points = 0*all_points;
for ith_point = 1:length(all_points(:,1))
    
    % Define start and end of the sensor
    sensor_vector_start = interior_point;
    sensor_vector_end = all_points(ith_point,:);
    
    % Define the start and end of the walls
    wall_start = walls(1:end-1,:);
    wall_end = walls(2:end,:);
    
    [distance,location,~] = ...
        INTERNAL_fcn_geometry_findIntersectionOfSegments(...
        wall_start,...
        wall_end,...
        sensor_vector_start,...
        sensor_vector_end);
    
    if ~isnan(distance)
        projected_points(ith_point,:) = location;
    else
        projected_points(ith_point,:) = sensor_vector_end;
    end
    
end
end


function [all_points, flag_was_intersection] = INTERNAL_fcn_findAllPoints(verticies,walls)

% Pad the vertices to wrap around, so we don't miss the last wall
verticies = [verticies; verticies(1,:)];
start_verticies = verticies(1:end-1,:);
end_verticies = verticies(2:end,:);


% Fill in blank all_points as starter
all_points = [];

flag_was_intersection = 0;
% Find all intersection points
for ith_point = 1:length(start_verticies(:,1))
    all_points = [all_points; start_verticies(ith_point,:)]; %#ok<AGROW>
    
    % Define start and end of the sensor
    sensor_vector_start = start_verticies(ith_point,:);
    sensor_vector_end = end_verticies(ith_point,:);
    
    % Define the start and end of the walls
    wall_start = walls(1:end-1,:);
    wall_end = walls(2:end,:);

    % Call a function that determines where and if the sensor crosses the
    % walls
    [distance,location,~] = ...
        INTERNAL_fcn_geometry_findIntersectionOfSegments(...
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


function verticies = INTERNAL_fcn_removeInfiniteVerticies(verticies,box)
% Goes through the verticies and removes infinite values by inserting
% points prior, and after the infinite one that "close" the polytope.

if any(isinf(verticies),'all')
    bad_indices = find(any(isinf(verticies),2));
    
    if length(bad_indices)>1
        warning('More than 2 infinities found in one vector. Code may break');
    end
    
    for ith_index = 1:length(bad_indices)
        bad_index = bad_indices(ith_index);
        
        % Find the prior and next points relative to the bad index point
        [prior_point, next_point, start_data, end_data] =...
            INTERNAL_fcn_findPriorNextPoints(bad_index,verticies);
     
        new_prior = fcn_MapGen_snapToAABB(box,prior_point);
        new_next  = fcn_MapGen_snapToAABB(box,next_point);
                      
        % Substitute data in, removing the infinite value
        verticies = [start_data; new_prior; new_next; end_data];
        
    end
end 

% Remove any vertices that are infinite, if any still remain (there should
% be none, so this is just in case)
flag_is_finite = ~isinf(verticies(:,1)).*~isinf(verticies(:,2));
verticies = verticies(flag_is_finite>0,:);

end % Ends function INTERNAL_fcn_removeInfiniteVerticies




function [prior_point, next_point, start_data, end_data] =...
    INTERNAL_fcn_findPriorNextPoints(index,verticies)
% Grab prior and next indices before and after an infinite value, being
% careful to check for situations where the infinite index is at start or
% end. Also creates vectors of points "start_data" and "end_data" which are
% the points before and after the start/end indices (inclusive)
Npoints = length(verticies(:,1));

prior_index = index-1;
next_index  = index+1;

% Is the index at the start?
if index == 1
    prior_index = Npoints;
    start_data = [];
else
    start_data = verticies(1:prior_index,:);
end

% Is the index at the end?
if index == Npoints
    next_index = 1;
    end_data = [];
else
    end_data = verticies(next_index:end,:);
end

% Fill the next points
prior_point = verticies(prior_index,:);
next_point  = verticies(next_index,:);
end % Ends INTERNAL_fcn_findPriorNextPoints


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
if 0 == flag_search_type
    good_vector = ((0<=t).*(1>=t).*(0<=u).*(1>=u));
elseif 1 == flag_search_type
    good_vector = ((0<=t).*(1>=t));
elseif 2 == flag_search_type
    good_vector = ((0<=t).*(1>=t).*(0<=u).*(1>=u));
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

