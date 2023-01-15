function flattened_polytopes = fcn_MapGen_flattenPolytopeMap(polytopes)
    % fcn_MapGen_flattenPolytopeMap
    % Given a structure array of polytopes that may overlap and/or be concave,
    % return a structure array of equivalent polytopes with overlapping regions
    % broken into new polytopes with traversal cost equal to the sum of the
    % traversal costs of the overlapping polytopes that created the region
    % additionally all polytopes will be broken up into triangles to enforce
    % convexity
    %
    %
    %
    % FORMAT:
    % flattened_polytopes = fcn_MapGen_flattenPolytopeMap(polytopes)
    %
    % INPUTS:
    %     polytopes - the initial polytope field, potentially containing concave
    %       and/or overlapping polytopes
    %
    %
    % OUTPUTS:
    %
    %
    %     flattened_polytopes - a polytope field with neither overlapping nor
    %       concave polytopes, that is equivalent to the oroginal polytope field
    %
    % DEPENDENCIES:
    %
    %     fcn_MapGen_plotPolytopes
    %     fcn_MapGen_fillPolytopeFieldsFromVertices
    %
    % EXAMPLES:
    %
    % See the script: script_fcn_MapGen_flattenPolytopeMap.m
    % for a full test suite.
    %
    % Questions or comments? contact sjh6473@psu.edu

    % REVISION HISTORY:
    % 2021_10_07
    % -- first written by Steve Harnett


    %% Debugging and Input checks
    flag_check_inputs = 1; % Set equal to 1 to check the input arguments
    flag_do_plot = 1;      % Set equal to 1 for plotting
    flag_do_debug = 0;     % Set equal to 1 for debugging
    if flag_do_debug
        fig_for_debug = 747;
        fig_num = 2100;
        st = dbstack; %#ok<*UNRCH>
        fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
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

    % convert polytopes to polyshapes
    if flag_do_plot
        fcn_MapGen_plotPolytopes(polytopes,1000,'-',2,[0 0 0],[],'square',[1 0 0 0 0.5]);
    end
    polyshapes = [];
    if flag_do_plot
        figure(1)
        clf
    end
    for i = 1:length(polytopes)
        polyshapes = [polyshapes, polyshape(polytopes(i).xv,polytopes(i).yv)];
        if flag_do_plot
            figure(1)
            hold on
            title('Polyshapes before subtraction')
            plot(polyshapes(i))
        end
    end
    % check each other polytope for overlap
    overlap_truth_table = overlaps(polyshapes);
    % remove self-overlaps from the diagonal (they don't count)
    overlap_truth_table = overlap_truth_table - eye(size(overlap_truth_table,1));
    % if an overlap is detected, subtract them
    [r, c] = find(overlap_truth_table==1);
    % check that at least one overlap exists
    if ~isempty(r) && ~isempty(c)
        % just look at the first overlap since this will be a recursive function
        p1 = polyshapes(r(1));
        % also note their costs since we'll need to flatten the cost into the new shape
        cost_of_p1 = polytopes(r(1)).cost;
        p2 = polyshapes(c(1));
        cost_of_p2 = polytopes(c(1)).cost;
        % make new polytope from overlapping region
        p3 = intersect(p1,p2);
        % modify old polytopes by removing overlapping region
        p1_new = subtract(p1,p3);
        p2_new = subtract(p2,p3);
        % tidy them
        p1_new = simplify(p1_new);
        p2_new = simplify(p2_new);
        p1_new = rmslivers(p1_new,0.001);
        p2_new = rmslivers(p2_new,0.001);
        p3 = simplify(p3);
        p3 = rmslivers(p3,0.001);
        if flag_do_plot
            figure(2)
            clf
            hold on
            plot(p1_new)
            plot(p2_new)
            plot(p3)
            title('overlapping polyshapes')
        end
        % make the polygonal polyshapes into series of triangular polyshapes
        % and triangular polytopes
        if p1_new.area > eps
            [~, p1_new_tris] = INTERNAL_fcn_triangulatePolyshape(p1_new,flag_do_plot);
        end
        if p2_new.area > eps
            [~, p2_new_tris] = INTERNAL_fcn_triangulatePolyshape(p2_new,flag_do_plot);
        end
        [~, p3_new_tris] = INTERNAL_fcn_triangulatePolyshape(p3,flag_do_plot);

        % TODO when making them triangles, set an ID field that the planner can use,
        % to know that when planning between two tris from the same parent poly,
        % the traversal counts as being interior

        % set cost to cost of triangles making up p1 and p2 to their original costs
        % set cost of intersection (p3) triangles to the sum of p1's and p2's costs
        if p1_new.area > eps
            p1_new_tris = fcn_polytope_editing_set_all_costs(p1_new_tris,cost_of_p1);
        end
        if p2_new.area > eps
            p2_new_tris = fcn_polytope_editing_set_all_costs(p2_new_tris,cost_of_p2);
        end
        p3_new_tris = fcn_polytope_editing_set_all_costs(p3_new_tris,cost_of_p1+cost_of_p2)
        % remove p1 and p2 from polytope list
        polytopes(r(1)) = [];
        polytopes(c(1)) = [];
        % append p1_new_tris, p2_new_tris, and p3_new_tris to polytopes list
        if p1_new.area > eps
            polytopes = [polytopes, p1_new_tris];
        end
        if p2_new.area > eps
            polytopes = [polytopes, p2_new_tris];
        end
        polytopes = [polytopes, p3_new_tris];
        if flag_do_plot
            fcn_MapGen_plotPolytopes(polytopes,1000,'-',2,[0 0 0],[],'square',[1 0 0 0 0.5]);
            title('originaly polytopes')
            if p1_new.area > eps
                fcn_MapGen_plotPolytopes(polytopes,1000,'-',2,[0 0 0],[],'square',[1 0 0 0 0.5]);
                title('triangulated parent 1, less, intersection')
            end
            if p2_new.area > eps
                fcn_MapGen_plotPolytopes(polytopes,1000,'-',2,[0 0 0],[],'square',[1 0 0 0 0.5]);
                title('triangulated parent 2, less intersection')
            end
            fcn_MapGen_plotPolytopes(polytopes,1000,'-',2,[0 0 0],[],'square',[1 0 0 0 0.5]);
            title('triangulated intersection')
        end
        % remove r(1), c(1) from truth table
        overlap_truth_table(r(1),c(1)) = 0;
        overlap_truth_table(c(1),r(1)) = 0;
        % recursively call this function again becaus the initial polytope list is
        % effectively different now
        polytopes = fcn_MapGen_flattenPolytopeMap(polytopes);
    end
    % if there are no intersections the loop will exit
    flattened_polytopes = polytopes;
    % check for very small polytopes
    small_poly_logical = extractfield(flattened_polytopes,'area')<1e-10;
    [small_poly_indecies] = find(small_poly_logical==1);
    flattened_polytopes(small_poly_indecies) = [];
    if true
        fcn_MapGen_plotPolytopes(polytopes,1000,'-',2,[0 0 0],[],'square',[1 0 0 0 0.5])
        title('final obstacle field')
    end
end

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

function [p_tri_polyshapes, p_tri_polytopes] = INTERNAL_fcn_triangulatePolyshape(my_polyshape,flag_do_plot)
    % make polyshape into triangulation
    p_tri = triangulation(my_polyshape);
    if flag_do_plot
        figure
        hold on
        triplot(p_tri)
    end
    % initialize cell array to store triangulation converted to polyshapes
    p_tri_polyshapes = {};
    % initialize array to store triangulation converted to polytopes
    p_tri_polytopes = [];
    % go through each set of connected triangle vertecies
    for i=1:size(p_tri.ConnectivityList,1)
        x1 = p_tri.Points(p_tri.ConnectivityList(i,1),1);
        y1 = p_tri.Points(p_tri.ConnectivityList(i,1),2);
        x2 = p_tri.Points(p_tri.ConnectivityList(i,2),1);
        y2 = p_tri.Points(p_tri.ConnectivityList(i,2),2);
        x3 = p_tri.Points(p_tri.ConnectivityList(i,3),1);
        y3 = p_tri.Points(p_tri.ConnectivityList(i,3),2);
%         if norm(([x1, y1] - [x2, y2]) == [0, 0]) == 0 || norm(([x1, y1] - [x3, y3]) == [0, 0]) == 0 || norm(([x2, y2] - [x3, y3]) == [0, 0]) == 0
%             continue
%         end
        % turn each triangle into a polyshape
        p_tri_polyshapes{i} = polyshape([x1 x2 x3],[y1 y2 y3]);
        if flag_do_plot
            plot(p_tri_polyshapes{i})
        end
        % turn each triangle into polytope
        p_tri_polytopes(i).vertices = [x1 y1; x2 y2; x3 y3; x1 y1];
    end
    % fill out all polytope fields from vertices
    p_tri_polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(p_tri_polytopes);
end

% example code doing different things
%    % for each polytope
%    for i = 1:length(polytopes)
%        for j = 1:length(polytopes)
%            if i ~= j
%            end
%        end
%    end
%
%    %% test case simple intersection
%    p1 = polyshape([1 0 0 1],[1 1 0 0]);
%    p2 = polyshape([1.5 0.5 0.5 1.5],[1.5 1.5 0.5 0.5]);
%    figure
%    hold on
%    plot(p1)
%    plot(p2)
%    p3 = intersect(p1,p2)
%    plot(p3)
%    p1_new = subtract(p1,p3)
%    p2_new = subtract(p2,p3)
%    figure
%    hold on
%    plot(p1_new)
%    plot(p2_new)
%    plot(p3)
%    p1_new = simplify(p1_new)
%    p2_new = simplify(p2_new)
%    p3 = simplify(p3)
%
%    %% test case intersecting and passing through
%    p1 = polyshape([1.4 0.4 0.4 1.4],[1.4 1.4 0.4 0.4]);
%    p2 = polyshape([0.6 0.6 1 1],[0 2 2 0]);
%    figure
%    hold on
%    plot(p1)
%    plot(p2)
%    p3 = intersect(p1,p2)
%    plot(p3)
%    p1_new = subtract(p1,p3)
%    p2_new = subtract(p2,p3)
%    figure
%    hold on
%    plot(p1_new)
%    plot(p2_new)
%    plot(p3)
%    p1_new = simplify(p1_new)
%    p2_new = simplify(p2_new)
%    p1_new = rmslivers(p1_new,0.001)
%    p2_new = rmslivers(p2_new,0.001)
%    p3 = simplify(p3)
%    % TODO p3 has to have cost of p1.cost + p2.cost
%    % triangulation code
%
%
%    %% test case enclave
%    p1 = polyshape([1.4 0.4 0.4 1.4],[1.4 1.4 0.4 0.4]);
%    p2 = polyshape([1.2 0.6 0.6 1.2],[1.2 1.2 0.6 0.6]);
%    figure
%    hold on
%    plot(p1)
%    plot(p2)
%    p3 = intersect(p1,p2)
%    plot(p3)
%    p1_new = subtract(p1,p3)
%    p2_new = subtract(p2,p3)
%    figure
%    hold on
%    plot(p1_new)
%    plot(p2_new)
%    plot(p3)
%    p1_new = simplify(p1_new)
%    p2_new = simplify(p2_new)
%    p3 = simplify(p3)
%    p1_tri = triangulation(p1_new)
%    figure
%    hold on
%    triplot(p1_tri)
%    p1_tri_polyshapes = {};
%    for i=1:size(p1_tri.ConnectivityList,1)
%        x1 = p1_tri.Points(p1_tri.ConnectivityList(i,1),1)
%        y1 = p1_tri.Points(p1_tri.ConnectivityList(i,1),2)
%        x2 = p1_tri.Points(p1_tri.ConnectivityList(i,2),1)
%        y2 = p1_tri.Points(p1_tri.ConnectivityList(i,2),2)
%        x3 = p1_tri.Points(p1_tri.ConnectivityList(i,3),1)
%        y3 = p1_tri.Points(p1_tri.ConnectivityList(i,3),2)
%        p1_tri_polyshapes{i} = ...
%            polyshape([x1 x2 x3],[y1 y2 y3]);
%        % make polyshape object into polytope by calling fill polytope from vertices
%        plot(p1_tri_polyshapes{i})
%    end
%
