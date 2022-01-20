function [shrunk_polytope] = ...
    fcn_MapGen_polytopeShrinkFromEdges(...
    shrinker,...
    edge_cut,...
    varargin)
% fcn_MapGen_polytopeShrinkFromEdges cuts edges off the polytopes
% Each edge is cut so that the entire polytope is trimmed exactly the same
% amount from each edge.
%
% FORMAT:
%
% [shrunk_polytopes,mu_final,sigma_final] = ...
%     fcn_MapGen_polytopeShrinkFromEdges(...
%     shrinker,...
%     edge_cut,...
%     (fig_num))
%
% INPUTS:
%
%     shrinker: original polytope with same fields as shrunk_polytopes
%     below
%
%     edge_cut: desired cut distance from each edge
%
%    (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%     SHRUNK_POLYTOPES: a 1-by-n seven field structure of shrunken polytopes,
%     where n <= number of polytopes with fields:
%       vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%         the number of the individual polytope vertices
%       xv: a 1-by-m vector of vertice x-coordinates
%       yv: a 1-by-m vector of vertice y-coordinates
%       distances: a 1-by-m vector of perimeter distances from one point to the
%         next point, distances(i) = distance from vertices(i) to vertices(i+1)
%       mean: centroid xy coordinate of the polytope
%       area: area of the polytope
%       max_radius: distance from the mean to the farthest vertex
%
% DEPENDENCIES:
%
%     fcn_MapGen_checkInputsToFunctions
%     fcn_MapGen_polytopeFindVertexAngles
%     fcn_MapGen_fillPolytopeFieldsFromVertices
%
% % EXAMPLES:
%
%
% For additional examples, see: script_test_fcn_MapGen_polytopeShrinkFromEdges
%
% This function was written on 2021_08_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%

% Revision History:
% 2021_08_02 - S. Brennan
% -- first write of code

% TO DO
% -- none

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 1;     % Set equal to 1 for debugging

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
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end

    % Check the shrinker input
    fcn_MapGen_checkInputsToFunctions(...
        shrinker, 'polytopes');

    % Check the edge_cut input
    fcn_MapGen_checkInputsToFunctions(...
        edge_cut, 'positive_1column_of_numbers',1);

end


% Does user want to show the plots?
if  3 == nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig_for_debug = 1584;
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

% find verticies, centroid, midpoints
vertices = shrinker.vertices;
centroid = shrinker.mean;
midpoints = (vertices(2:end,:)+vertices(1:end-1,:))/2;

% Plot results?
if flag_do_debug
    figure(fig_for_debug);
    clf;
    axis equal;
    hold on;

    % Plot the polytope in red
    plot(vertices(:,1),vertices(:,2),'r-','Linewidth',2);

    % Find size of vertices
    size = max(max(vertices)) - min(min(vertices));
    nudge = size*0.003;

    % Number the midpoints of vertices as labels
    for ith_midpoint = 1:length(midpoints(:,1))
        text(midpoints(ith_midpoint,1),midpoints(ith_midpoint,2),...
            sprintf('%.0d',ith_midpoint));
    end
end

% % Calculate the Voronoi diagram
% midpoints = (vertices(2:end,:)+vertices(1:end-1,:))/2;
% [V,C] = voronoin(midpoints);
% voronoi_vertices = V;
% if flag_do_debug
%
%     % Plot the voronoi_vertices
%     plot(voronoi_vertices(:,1),voronoi_vertices(:,2),'r.','Markersize',10);
%
% end

%% Find vertex skeleton
% This is the set of vertices created with different edge cuts
[new_vertices, new_direction, new_scale_factor] = ...
    INTERNAL_fcn_findVertexSkeleton(vertices,23123);


% Calculate the directions needed for projection
[unit_direction, scale_factor, half_angles] = INTERNAL_fcn_findUnitDirectionVectors(vertices,fig_for_debug);
distances = edge_cut.*scale_factor; % distance needed to move each point
projection_vectors = unit_direction.*distances;

%% Find new vertices based on projection
short_new_vertices = vertices(1:end-1,:)+projection_vectors;
new_vertices = [short_new_vertices; short_new_vertices(1,:)];

%% Check if cannot shrink
cannot_shrink = 0;
% If projection greater than the vertex, then can't shrink
% Find the vectors pointing from vertices to mean
vectors_vertices_to_mean = centroid - vertices(1:end-1,:);
lengths_vertices_to_mean = sum(vectors_vertices_to_mean.^2,2).^0.5;
unit_vectors_vertices_to_mean = vectors_vertices_to_mean./lengths_vertices_to_mean;

% Find the projection vector length in direction of mean, using dot product
projection_toward_mean = dot(projection_vectors,unit_vectors_vertices_to_mean,2);

if any(projection_toward_mean>lengths_vertices_to_mean)
    cannot_shrink = 1;
end

% Check the angles
new_angles = fcn_MapGen_polytopeFindVertexAngles(...
    new_vertices);



% Check if any angles are negative. If so, this is a self-intersection case
% and we can fix it.
if any(new_angles<0)

    if flag_do_debug
        figure(fig_for_debug);
        hold on;
        plot(new_vertices(:,1),new_vertices(:,2),'g-','Linewidth',2);

        % Find size of vertices
        size = max(max(new_vertices)) - min(min(new_vertices));
        nudge = size*0.01;

        % Label the vertices
        for ith_angle = 1:length(new_angles(:,1))
            ith_vertex = ith_angle;
            text(new_vertices(ith_vertex,1)+nudge,new_vertices(ith_vertex,2),...
                sprintf('%.0d = %.0f deg',ith_vertex, new_angles(ith_angle,1)*180/pi));
        end
    end

    % Find any self intersections
    [vertices_with_self_intersects] = ...
        fcn_MapGen_polytopeFindSelfIntersections(...
        new_vertices);

    if flag_do_debug
        figure(fig_for_debug);
        hold on;
        plot(vertices_with_self_intersects(:,1),vertices_with_self_intersects(:,2),'kx','Linewidth',2);
    end

    % Keep only the vertices that are positive. The hard part is that the
    % center of the polytope is hard to figure out since there will be
    % multiple polytopes. To find the new center, just match up the angles
    % from the new one to the old


    % Find the angles for the self-intersection points.
    self_intersection_angles = fcn_MapGen_polytopeFindVertexAngles(...
        vertices_with_self_intersects);
    angles = new_angles;
    angles*180/pi
    self_intersection_angles*180/pi


    % Which self-intersection angles match the old angles?
    Nvertices = length(vertices_with_self_intersects(:,1))-1;
    good_vertices = find(ismember(round(self_intersection_angles*1000),round(angles*1000)));

    % For the solution to be valid, at least 2 of the good vertices must be
    % consecutive. Must remember to check the end rollover condition as well
    consecutive_vertices = diff([good_vertices; good_vertices(1)+Nvertices]);
    is_good = find(consecutive_vertices==1);

    if isempty(is_good)
        cannot_shrink = 1;
    else

        good_consecutive_vertices = good_vertices(is_good);
        next_after = mod(good_consecutive_vertices,Nvertices)+1;
        good_consecutive_vertices = unique([good_consecutive_vertices;next_after]);
    end

    if 0==cannot_shrink

        % The points to keep will be these angles, plus one of the other
        % points leading in.
        next_point_after_good = good_consecutive_vertices+1;
        prev_point_before_good = mod(good_consecutive_vertices-2,Nvertices)+1;
        points_after = vertices_with_self_intersects(next_point_after_good,:);
        points_before = vertices_with_self_intersects(prev_point_before_good,:);
        all_points = [points_after; points_before];

        % Now we need to find which of these points is repeated.
        indices = ones(length(all_points(:,1)),1);
        [~,IA,~] = unique(all_points,'rows');
        indices(IA) = 0;
        repeated_index = indices>0;
        repeated_points = all_points(repeated_index,:);

        % Fill in the vertices to average
        vertices_to_average = [...
            vertices_with_self_intersects(good_consecutive_vertices,:);
            repeated_points];

        average_point = mean(vertices_to_average,1);

        % Check to see if the average point is further than distance of
        % projection
        vectors_vertices_to_average = average_point - vertices(1:end-1,:);
        lengths_vertices_to_average = sum(vectors_vertices_to_average.^2,2).^0.5;

        if flag_do_debug
            figure(fig_for_debug);
            hold on;
            plot(vertices_to_average(:,1),vertices_to_average(:,2),'mo','Markersize',10);
        end

        % Check to see if we can shrink. If not, no sense in continuing
        % calculations!
        tol = 0.0001;
        if any(distances+tol>=lengths_vertices_to_average)
            cannot_shrink = 1;
        end
    end

    if 0==cannot_shrink
        % Push points onto nearest wall
        [projected_points] = ...
            fcn_MapGen_polytopeProjectVerticesOntoWalls(...,
            average_point,...
            vertices_with_self_intersects,...
            vertices_with_self_intersects(1:end-1,:),...
            vertices_with_self_intersects(2:end,:));

        % Crop the vertices
        try
            [cropped_vertices] = ...
                fcn_MapGen_polytopeRemoveColinearVertices(...,
                projected_points);
            if ~isempty(cropped_vertices)
                cropped_vertices = [cropped_vertices; cropped_vertices(1,:)];
            else
                % Degenerate case where all the points are the same
                cropped_vertices = projected_points(1:3,:);
            end
        catch
            disp('Stop here');
        end
        % Check the angles again
        new_angles = fcn_MapGen_polytopeFindVertexAngles(...
            cropped_vertices);
        new_vertices = cropped_vertices;

        % Plot the results?
        if flag_do_debug
            figure(fig_for_debug);
            hold on;
            plot(cropped_vertices(:,1),cropped_vertices(:,2),'b-','Linewidth',2);

        end % Ends debug plotting
    end% Ends cannot shrink
end

% Check if overlap occurs
if abs(2*pi-sum(new_angles))>1E-5
    cannot_shrink = 1;
end

if 1==cannot_shrink
    new_vertices = ones(length(vertices(:,1)),1)*centroid;
end

%% Fill in the results
shrunk_polytope.vertices = new_vertices;

% fill in other fields from the vertices field
shrunk_polytope = fcn_MapGen_fillPolytopeFieldsFromVertices(shrunk_polytope);



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

    % Plot the cetroid in black
    plot(centroid(:,1),centroid(:,2),'ko','Markersize',10);

    % Plot the input shrinker in red
    fcn_MapGen_plotPolytopes(shrinker,fig_num,'r',2);

    % plot the output polytope in blue
    fcn_MapGen_plotPolytopes(shrunk_polytope,fig_num,'b',2);

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function

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
    fcn_MapGen_checkInputsToFunctions(...
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
    % TODO(@sjharnett) understand how this could give interior vertex normal vectors
    % Draw the unit vectors in the vertex direction
    quiver(vertices(1:end-1,1),vertices(1:end-1,2),unit_vectors_vertex_to_vertex(:,1),unit_vectors_vertex_to_vertex(:,2),0);
    % TODO(@sjharnett) understand how this could give interior vertex normal vectors
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends INTERNAL_fcn_findUnitDirectionVectors


function [new_vertices, new_projection_vectors, cut_distance] = ...
    INTERNAL_fcn_findVertexSkeleton(vertices, varargin)
% Calculates the VertexSkeleton for a polytope


% Revision History:
% 2021_08_06 - S. Brennan
% -- first write of code

% TO DO
% -- none

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 1;     % Set equal to 1 for debugging

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
    fcn_MapGen_checkInputsToFunctions(...
        vertices, '2column_of_numbers');

end


% Does user want to show the plots?
if  2 == nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig_for_debug = 333;
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
iteration = 1;
flag_stop_loop = 0;
total_cut = 0;
working_vertices = vertices;

while 0 == flag_stop_loop
    new_vertices{iteration} = working_vertices; %#ok<AGROW>
    cut_distance{iteration} = total_cut; %#ok<AGROW>

    Nvertices = length(working_vertices(:,1));
    if 2==Nvertices
        new_projection_vectors{iteration} = [0 0; 0 0]; %#ok<AGROW>
        flag_stop_loop = 1;
    else
        [vector_direction_of_unit_cut, ...
            half_angles,...
            distances_vertex_to_vertex,...
            unit_vectors_vertex_to_vertex] = ...
            INTERNAL_fcn_findUnitDirectionVectors(working_vertices,333);
        new_projection_vectors{iteration} = vector_direction_of_unit_cut; %#ok<AGROW>

        % Find the projection point
        looped_half_angles = [half_angles; half_angles(1)];
        theta2s = looped_half_angles(2:end);
        theta1s = looped_half_angles(1:end-1);

        A = tan(theta2s)./tan(theta1s);
        L1 = distances_vertex_to_vertex.*(A./(A+1));
        projection_points = working_vertices(1:end-1,:) + L1.*unit_vectors_vertex_to_vertex;

        if flag_do_debug
            figure(fig_for_debug);
            grid on
            grid minor
            hold on
            axis equal

            figure(333);
            plot(projection_points(:,1),projection_points(:,2),'b.','Markersize',20);
        end

        % Find the cut size
        Lcuts = distances_vertex_to_vertex.*tan(theta1s).*tan(theta2s)./(tan(theta1s)+tan(theta2s));
        projection_directions = unit_vectors_vertex_to_vertex*[0 1; -1 0]; % Rotate by 90 degrees
        intersection_points = projection_points + Lcuts.*projection_directions;

        if flag_do_debug
            plot(intersection_points(:,1),intersection_points(:,2),'g.','Markersize',20);
            % TODO(@sjharnett) understand how this could be used to find interior vertex normal vectors
        end

        % Find the tightest cut that is possible, and use this to update
        % the total cut
        min_cut = min(Lcuts);
        total_cut = total_cut + min_cut;

        % Sort by the cut size, finding the cuts that are nearly exactly
        % the same cut length
        indices_repeated = find(Lcuts<(min_cut+1E5*eps));

        % Tag the vertices that are merged
        all_indices = (1:Nvertices)';
        indices_following = indices_repeated+1;
        % indices_following = mod(indices_repeated,Nvertices-1)+1;

        vertices_merged = union(indices_repeated,indices_following);
        vertices_not_merged = ~ismember(all_indices(1:end-1),vertices_merged);

        % Associate the indices with the repeats
        for ith_repeat = 1:length(indices_following)
            current_index = indices_following(ith_repeat);
            all_indices(current_index) = all_indices(current_index-1);
        end

        % Check for the rollover condition, since the last point is first
        % point
        if all_indices(end)~=Nvertices
            first_indices = all_indices==1;
            all_indices(first_indices) = all_indices(end);
        end

        % Crop back all indices to avoid rollover end point
        all_indices = all_indices(1:(Nvertices-1),:);
        vertices_merged = vertices_merged(vertices_merged<Nvertices);


        % vertex_indices_to_merge = indices_repeated;
        %         points_to_merge = intersection_points(indices_repeated,:);
        %         for ith_index = 1:length(indices_repeated)
        %             current_index = indices_repeated(ith_index);
        %             diff_to_current = points_to_merge-points_to_merge(current_index,:);
        %             squared_distance_to_current = sum(diff_to_current.^2,2);
        %             same_as_current = squared_distance_to_current<1E-14;
        %             indices_repeated(same_as_current) = current_index;
        %         end
        %         intersection_points(vertex_indices_to_merge,:) = ...
        %             intersection_points(indices_repeated,:);
        % all_indices(vertex_indices_to_merge)=indices_repeated;

        % vertices_not_merged = ~ismember(all_indices,vertex_indices_to_merge);


        % Calculate the movements
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
    figure(fig_num);
    clf;
    grid on
    grid minor
    hold on
    axis equal

    % Plot the polytope in red
    plot(vertices(:,1),vertices(:,2),'r-','Linewidth',2);

    % Find size of vertices
    size = max(max(vertices)) - min(min(vertices));
    nudge = size*0.003;

    % Number the vertices with labels
    for ith_vertex = 1:length(vertices(:,1))
        text(vertices(ith_vertex,1)+nudge,vertices(ith_vertex,2),...
            sprintf('%.0d',ith_vertex));
    end

    % Plot each contraction
    for ith_contraction = 1:(length(new_vertices)-1)
        starting_points = new_vertices{ith_contraction};
        starting_points = starting_points(1:end-1,:);
        ending_points = starting_points + ...
            new_projection_vectors{ith_contraction}*...
            (cut_distance{ith_contraction+1}-cut_distance{ith_contraction});

        for jth_segment = 1:length(starting_points(:,1))
            plot(...
                [starting_points(jth_segment,1) ending_points(jth_segment,1)],...
                [starting_points(jth_segment,2) ending_points(jth_segment,2)],...
                '-');
        end
        plot(starting_points(:,1),starting_points(:,2),'.','Markersize',20);
        plot(ending_points(:,1),ending_points(:,2),'.','Markersize',20);

    end

end % Ends flag_do_plot if statement

end % ends fucntion INTERNAL_fcn_findVertexSkeleton

function [new_vertices, new_direction, new_scale_factor] = INTERNAL_fcn_findIntersectionOfVertexProjections(vertices,shrinker)
% creates skeleton of where vertecies shrink to the point of
% Finds intersection points between adjacent verticies
longest_distance = 2*shrinker.max_radius;
% find vertex projection and store in unit mean vector
longest_projections = unit_mean_vectors*longest_distance;
Nangles = length(angles);
distances_hit = zeros(Nangles,1);
first_wall_hit = zeros(Nangles,1);


for jth_angle = 1:Nangles
    prev_angle = mod(jth_angle-2,Nangles)+1;
    next_angle = mod(jth_angle,Nangles)+1;

    % Find hits on the adjacent walls
    wall_start = [vertices(prev_angle,:); vertices(next_angle,:)];
    wall_end = wall_start + ...
        [longest_projections(prev_angle,:);
        longest_projections(next_angle,:)];
    sensor_vector_start = vertices(jth_angle,:);
    sensor_vector_end = sensor_vector_start + ...
        longest_projections(jth_angle,:);

    [distances_hit(jth_angle),location_of_hit,wall_hit] = ...
        fcn_MapGen_findIntersectionOfSegments(...
        wall_start,...
        wall_end,...
        sensor_vector_start,...
        sensor_vector_end);

    vertex_hit = mod(jth_angle -1 + (wall_hit-1.5)*2,Nangles)+1;
    first_wall_hit(jth_angle) = vertex_hit;

    if flag_do_debug
        plot(location_of_hit(:,1),location_of_hit(:,2),'r.');
        text(location_of_hit(:,1),location_of_hit(:,2),sprintf('V %.0d hit %.0d',jth_angle,vertex_hit));
    end

end
[~,closest_index] = min(distances_hit);



end % Ends function INTERNAL_fcn_findIntersectionOfVertexProjections
