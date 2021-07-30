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
flag_check_inputs = 1; % Set equal to 1 to check the input arguments 
flag_do_plot = 0;      % % Set equal to 1 for plotting 
flag_do_debug = 0;     % Set equal to 1 for debugging 

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

tolerance = 0.001;
location = [0.9935 0.1354];
if (...
        (interior_point(1,1)<location(1)+tolerance) && ...
        (interior_point(1,1)>location(1)-tolerance) && ...
        (interior_point(1,2)<location(2)+tolerance) && ...
        (interior_point(1,2)>location(2)-tolerance))
    disp('stop here');
end

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
% % Are any vertices infinite? If so, we need to check that the adjacent
% % vertices will create a reasonable polytope. 
% [vertices_no_infinite,~] = ...
%     INTERNAL_fcn_removeInfiniteVertices(vertices,AABB,walls);
% 
% % Open the figure if doing debugging
% if flag_do_debug
%     % Plot the vertices
%     figure(1);
%     plot(...
%         [vertices_no_infinite(:,1); vertices_no_infinite(1,1)],...
%         [vertices_no_infinite(:,2); vertices_no_infinite(1,2)],...
%         '.-');
% end

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
    projected_points = ...
        INTERNAL_fcn_projectAllPointsOntoWalls(interior_point, all_points,walls);
    
    if flag_do_debug
        figure(1);
        % Plot the projected_points locations
        plot(projected_points(:,1),projected_points(:,2),'go-');
    end
    
    %     % Check to see if corners need to be added
    %     projected_with_corner_points = ...
    %         INTERNAL_fcn_addCorners(projected_points,flag_infinite_was_found,walls);
    %
    %     if flag_do_debug
    %         figure(1);
    %         % Plot the projected_points locations
    %         plot(projected_with_corner_points(:,1),projected_with_corner_points(:,2),'co-');
    %     end
    
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
% 
% function projected_with_corner_points = ...
%     INTERNAL_fcn_addCorners(vertices,flag_infinite_was_found,walls)
% 
% %
% % % Check if any of the corners, e.g. where the walls start, are inside the
% % % polytope. If so, need to add these.
% % for ith_wall = 1:length(walls(1:end-1,1))
% %     test_point = walls(ith_wall,:);
% %     [ in_polytope ] = ...
% %         fcn_MapGen_checkIfPointInsideConvexPolytope( ...
% %         test_point, ...
% %         vertices);
% %     if in_polytope
% %         all_points = [all_points; test_point]; %#ok<AGROW>
% %     end
% %
% % end
% 
% if flag_infinite_was_found
%     
%     % Is the point is in a corner? If so, add an extra point for the
%     % corner
%     extra_point = [];
%     if all(isinf(vertices(bad_index,:)))
%         % Calculate the angles covered by the vertices
%         end_points_offset = [prior_point; next_point] - mean([box(1:2);box(3:4)],1);
%         [angles, ~] = cart2pol(end_points_offset(:,1),end_points_offset(:,2));
%         min_angle = min(angles);
%         max_angle = max(angles);
%         
%         % Check to see if angle crosses over -180 degrees
%         if (max_angle-min_angle)>pi
%             angles = mod(angles,2*pi);
%             min_angle = min(angles);
%             max_angle = max(angles);
%         end
%         
%         box_angle_right = abs(atan2(box(4)-box(2),box(3)-box(1)));
%         box_angle_left  = pi-box_angle_right;
%         
%         % Find location of the corners, and add them
%         if min_angle<-box_angle_left && max_angle>= -box_angle_left
%             % Bottom left corner
%             extra_point = [box(1) box(2)];
%         elseif min_angle<-box_angle_right && max_angle>= -box_angle_right
%             % Bottom right corner
%             extra_point = [box(3) box(2)];
%         elseif min_angle<box_angle_left && max_angle>= box_angle_left
%             % Top left corner
%             extra_point = [box(1) box(4)];
%         elseif min_angle<box_angle_right && max_angle>= box_angle_right
%             % Top right corner
%             extra_point = [box(3) box(4)];
%         end
%     end % Ends if for all vertices bad
%     projected_with_corner_points = [start_data; new_prior; extra_point; new_next; end_data];
% else
%     projected_with_corner_points = vertices;
% end % Ends if for flag_infinite_was_found
% 
% 
% end

function cropped_vertices = INTERNAL_fcn_cropRepeatedPoints(projected_points_with_repeats)
% Remove repeats
[projected_points,~,~] = unique(projected_points_with_repeats,'rows','stable');

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



function projected_points = INTERNAL_fcn_projectAllPointsOntoWalls(interior_point, all_points_unsorted,walls)
% From the interior point, project all_points back onto the wall.

% Convert the input xy data into polar coordinates, so that we can sort by
% theta. If we don't do this, then sometimes (rarely) the walls of the
% polytope will jump around producing very weird results.
offset_points = all_points_unsorted - interior_point;
[theta,~]=cart2pol(offset_points(:,1),offset_points(:,2));
[~,index_sorted] = sort(theta,'descend');
all_points = all_points_unsorted(index_sorted,:);

projected_points = 0*all_points;
for ith_point = 1:length(all_points(:,1))
    
    % Define start and end of the sensor
    sensor_vector_start = interior_point;
    sensor_vector_end = all_points(ith_point,:);
    
    % Define the start and end of the walls
    wall_start = walls(1:end-1,:);
    wall_end = walls(2:end,:);
    
    [distance,location,~] = ...
        fcn_MapGen_findIntersectionOfSegments(...
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


% function [resulting_vertices,flag_infinite_was_found] = ...
%     INTERNAL_fcn_removeInfiniteVertices(vertices,AABB,walls)
% % Goes through the vertices and removes infinite values by inserting
% % points prior, and after the infinite one that "close" the polytope.
% 
% flag_do_debug = 0;
% if flag_do_debug
%     figure(77474);
%     clf;
%     hold on; 
% 
%     % Plot the vertices
%     plot(vertices(:,1),vertices(:,2),'r.-');
%     
%     % Plot the walls
%     plot(walls(:,1),walls(:,2),'k-');
% end
% 
% flag_infinite_was_found = 0;
% if any(isinf(vertices),'all') % Are there any infinite vertices
%     warning('Infinite vertices found ... this should not have happened as a prior function should have cleaned these!');
%     %     flag_infinite_was_found = 1;
%     %
%     %     bad_indices = find(any(isinf(vertices),2));
%     %
%     %     % Warn user if 2 infinite values found. This may cause the code to fail
%     %     % because it searches for the points before and after infinity,
%     %     % assuming these points are NOT infinite.
%     %     if length(bad_indices)>1
%     %         % Check to see if the infinity is at start and end, artificially
%     %         % repeated
%     %         if isequal(bad_indices,[1; length(vertices(:,1))])
%     %             vertices_no_repeats = vertices(1:end-1,:);
%     %         else
%     %             warning('More than 2 infinities found in one vector. Code may break');
%     %         end
%     %     else
%     %         vertices_no_repeats = vertices(1:end-1,:);
%     %     end
%     %
%     %     bad_index = bad_indices(1);
%     %     % Rearrange the points so that the infinite index is the first one.
%     %     % Makes things easier in later steps since we don't have to carry
%     %     % around two snips of data, just one
%     %     vertex_string = ...
%     %         [vertices_no_repeats(bad_index+1:end,:); ...
%     %         vertices_no_repeats(1:bad_index-1,:)];
%     %
%     %     if flag_do_debug
%     %         % Plot the vertex_string
%     %         plot(vertex_string(:,1),vertex_string(:,2),'b.-');
%     %     end
%     % else
%     %     vertex_string = vertices;
%     % end
%     %
%     % % Now we crop the vertices
%     % [cropped_vertices,~] = ...
%     %     fcn_MapGen_cropVerticesByWallIntersections(vertex_string,walls);
%     %
%     % if flag_do_debug
%     %     % Plot the cropped_vertices
%     %     plot(cropped_vertices(:,1),cropped_vertices(:,2),'g.-');
%     % end
%     %
%     %     % Find the prior and next points relative to the bad index point
%     %     % The prior_point and next_point are not used as vertices themselves,
%     %     % but rather to create new vertices by snapping to the bounding polygon
%     %     % or ABB. Thus, the data before and after these points, including these
%     %     % points, must be kept.
%     %     prior_point = cropped_vertices(end,:);
%     %     next_point = cropped_vertices(1,:);
%     %     start_data = [];
%     %     end_data = cropped_vertices(1:end,:);
%     %
%     %     [new_prior, ~] = fcn_MapGen_snapToAABB(AABB,prior_point);
%     %     [new_next, ~]  = fcn_MapGen_snapToAABB(AABB,next_point);
%     %
%     %
%     %     % Substitute data in, removing the infinite value
%     %     resulting_vertices = [start_data; new_prior; new_next; end_data];
%     %
%     %     if flag_do_debug
%     %         % Plot the resulting_vertices
%     %         plot(resulting_vertices(:,1),resulting_vertices(:,2),'g.-');
%     %     end
%     %
%     %     %         % Check if any of the corners, e.g. where the walls start, are inside the
%     %     %         % polytope. If so, need to add these.
%     %     %         for ith_wall = 1:length(walls(1:end-1,1))
%     %     %             test_point = walls(ith_wall,:);
%     %     %             [ in_polytope ] = ...
%     %     %                 fcn_MapGen_checkIfPointInsideConvexPolytope( ...
%     %     %                 test_point, ...
%     %     %                 vertices);
%     %     %             if in_polytope
%     %     %                 all_points = [all_points; test_point]; %#ok<AGROW>
%     %     %             end
%     %     %
%     %     %         end
%     %
%     %     % Is the point is in a corner? If so, add an extra point for the
%     %     % corner
%     %     extra_point = [];
%     %     if all(isinf(vertices(bad_index,:)))
%     %         % Calculate the angles covered by the vertices
%     %         end_points_offset = [prior_point; next_point] - mean([AABB(1:2);AABB(3:4)],1);
%     %         [angles, ~] = cart2pol(end_points_offset(:,1),end_points_offset(:,2));
%     %         min_angle = min(angles);
%     %         max_angle = max(angles);
%     %
%     %         % Check to see if angle crosses over -180 degrees
%     %         if (max_angle-min_angle)>pi
%     %             angles = mod(angles,2*pi);
%     %             min_angle = min(angles);
%     %             max_angle = max(angles);
%     %         end
%     %
%     %         box_angle_right = abs(atan2(AABB(4)-AABB(2),AABB(3)-AABB(1)));
%     %         box_angle_left  = pi-box_angle_right;
%     %
%     %         % Find location of the corners, and add them
%     %         if min_angle<-box_angle_left && max_angle>= -box_angle_left
%     %             % Bottom left corner
%     %             extra_point = [AABB(1) AABB(2)];
%     %         elseif min_angle<-box_angle_right && max_angle>= -box_angle_right
%     %             % Bottom right corner
%     %             extra_point = [AABB(3) AABB(2)];
%     %         elseif min_angle<box_angle_left && max_angle>= box_angle_left
%     %             % Top left corner
%     %             extra_point = [AABB(1) AABB(4)];
%     %         elseif min_angle<box_angle_right && max_angle>= box_angle_right
%     %             % Top right corner
%     %             extra_point = [AABB(3) AABB(4)];
%     %         end
%     %     end
%     %
%     %     non_repeating_resulting_vertices = [start_data; new_prior; extra_point; new_next; end_data];
%     %     resulting_vertices = [non_repeating_resulting_vertices; non_repeating_resulting_vertices(1,:)];
% else
%     resulting_vertices = vertices;
% end 
% 
% end % Ends function INTERNAL_fcn_removeInfiniteVertices



% function [prior_point, next_point, start_data, end_data] =...
%     INTERNAL_fcn_findPriorNextPoints(index,vertices)
% % Grab prior and next indices before and after an infinite value, being
% % careful to check for situations where the infinite index is at start or
% % end. Also creates vectors of points "start_data" and "end_data" which are
% % the points before and after the start/end indices (inclusive)
% Npoints = length(vertices(:,1));
% 
% prior_index = index-1;
% next_index  = index+1;
% 
% % Is the index at the start?
% if index == 1
%     prior_index = Npoints;
%     start_data = [];
% else
%     start_data = vertices(1:prior_index,:);
% end
% 
% % Is the index at the end?
% if index == Npoints
%     next_index = 1;
%     end_data = [];
% else
%     end_data = vertices(next_index:end,:);
% end
% 
% % Fill the next points
% prior_point = vertices(prior_index,:);
% next_point  = vertices(next_index,:);
% end % Ends INTERNAL_fcn_findPriorNextPoints

