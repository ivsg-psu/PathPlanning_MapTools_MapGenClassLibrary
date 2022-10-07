clear all; close all; clc
function flattened_polytopes = fcn_MapGen_flattenPolytopeMap(polytopes)
    % convert polytopes to polyshapes
    polyshapes = NaN(1,length(polytopes));
    for i = 1:length(polytopes)
        polyshapes(i) = polyshape(polytopes(i).xv,polytopes(i).yv);
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
        p2 = polyshapse(c(1));
        cost_of_p2 = polytopes(c(1)).cost;
        % make new polytope from overlapping region
        p_3 = intersect(p1,p2);
        % modify old polytopes by removing overlapping region
        p1_new = subtract(p1,p3);
        p2_new = subtract(p2,p3);
        % tidy them
        p1_new = simplify(p1_new)
        p2_new = simplify(p2_new)
        p1_new = rmslivers(p1_new,0.001)
        p2_new = rmslivers(p2_new,0.001)
        p3 = simplify(p3)
        p2_new = rmslivers(p2_new,0.001)
        % make them triangles
        p1_new_tris = INTERNAL_fcn_triangulatePolyshape(p1_new);
        p2_new_tris = INTERNAL_fcn_triangulatePolyshape(p2_new);
        p3_new_tris = INTERNAL_fcn_triangulatePolyshape(p3);
        % make p1_new_tris into list of polytopes (do this below in internal function)
        % set cost to cost of p1
        % repeat for p2_new_tris and p3_new_tris
        % remove p1 and p2 from polytope list
        % append p1_new_tris, p2_new_tris, and p3_new_tris to polytopes list
        % recursively call this function again becaus the initial polytope list is
        % effectively different now
    % if there are no intersections the loop will exit
end
function p_tri_polyshapes = INTERNAL_fcn_triangulatePolyshape(polyshape)
    p_tri = triangulation(polyshape)
    if flag_do_plot
        figure
        hold on
        triplot(p_tri)
    end
    p_tri_polyshapes = {};
    for i=1:size(p_tri.ConnectivityList,1)
        x1 = p_tri.Points(p_tri.ConnectivityList(i,1),1)
        y1 = p_tri.Points(p_tri.ConnectivityList(i,1),2)
        x2 = p_tri.Points(p_tri.ConnectivityList(i,2),1)
        y2 = p_tri.Points(p_tri.ConnectivityList(i,2),2)
        x3 = p_tri.Points(p_tri.ConnectivityList(i,3),1)
        y3 = p_tri.Points(p_tri.ConnectivityList(i,3),2)
        p_tri_polyshapes{i} = ...
            polyshape([x1 x2 x3],[y1 y2 y3]);
        if flag_do_plot
            plot(p_tri_polyshapes{i})
        end
    end
end

% example code doing different things
    % for each polytope
    for i = 1:length(polytopes)
        for j = 1:length(polytopes)
            if i ~= j
            end
        end
    end

    %% test case simple intersection
    p1 = polyshape([1 0 0 1],[1 1 0 0]);
    p2 = polyshape([1.5 0.5 0.5 1.5],[1.5 1.5 0.5 0.5]);
    figure
    hold on
    plot(p1)
    plot(p2)
    p3 = intersect(p1,p2)
    plot(p3)
    p1_new = subtract(p1,p3)
    p2_new = subtract(p2,p3)
    figure
    hold on
    plot(p1_new)
    plot(p2_new)
    plot(p3)
    p1_new = simplify(p1_new)
    p2_new = simplify(p2_new)
    p3 = simplify(p3)

    %% test case intersecting and passing through
    p1 = polyshape([1.4 0.4 0.4 1.4],[1.4 1.4 0.4 0.4]);
    p2 = polyshape([0.6 0.6 1 1],[0 2 2 0]);
    figure
    hold on
    plot(p1)
    plot(p2)
    p3 = intersect(p1,p2)
    plot(p3)
    p1_new = subtract(p1,p3)
    p2_new = subtract(p2,p3)
    figure
    hold on
    plot(p1_new)
    plot(p2_new)
    plot(p3)
    p1_new = simplify(p1_new)
    p2_new = simplify(p2_new)
    p1_new = rmslivers(p1_new,0.001)
    p2_new = rmslivers(p2_new,0.001)
    p3 = simplify(p3)
    % TODO p3 has to have cost of p1.cost + p2.cost
    % triangulation code


    %% test case enclave
    p1 = polyshape([1.4 0.4 0.4 1.4],[1.4 1.4 0.4 0.4]);
    p2 = polyshape([1.2 0.6 0.6 1.2],[1.2 1.2 0.6 0.6]);
    figure
    hold on
    plot(p1)
    plot(p2)
    p3 = intersect(p1,p2)
    plot(p3)
    p1_new = subtract(p1,p3)
    p2_new = subtract(p2,p3)
    figure
    hold on
    plot(p1_new)
    plot(p2_new)
    plot(p3)
    p1_new = simplify(p1_new)
    p2_new = simplify(p2_new)
    p3 = simplify(p3)
    p1_tri = triangulation(p1_new)
    figure
    hold on
    triplot(p1_tri)
    p1_tri_polyshapes = {};
    for i=1:size(p1_tri.ConnectivityList,1)
        x1 = p1_tri.Points(p1_tri.ConnectivityList(i,1),1)
        y1 = p1_tri.Points(p1_tri.ConnectivityList(i,1),2)
        x2 = p1_tri.Points(p1_tri.ConnectivityList(i,2),1)
        y2 = p1_tri.Points(p1_tri.ConnectivityList(i,2),2)
        x3 = p1_tri.Points(p1_tri.ConnectivityList(i,3),1)
        y3 = p1_tri.Points(p1_tri.ConnectivityList(i,3),2)
        p1_tri_polyshapes{i} = ...
            polyshape([x1 x2 x3],[y1 y2 y3]);
        % make polyshape object into polytope by calling fill polytope from vertices
        plot(p1_tri_polyshapes{i})
    end
