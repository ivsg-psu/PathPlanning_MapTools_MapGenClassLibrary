clear all; close all; clc
function flattened_polytopes = fcn_MapGen_flattenPolytopeMap(polytopes)
    % convert polytopes to polyshapes
    polyshapes = [];
    for i = 1:length(polytopes)
        polytopes(i).xv,
    % for each polytope
    % check each other polytope for overlap
    % if an overlap is detected, subtract them
    % tidy them (make them triangles)
    % add p3 to polytopes list
    % add new p1 and new p2 to list
    % delete original p1 and p2
    % recursively call this function again becaus the initial polytope list is
    % effectively different now
    % if there are no intersections the loop will exit
    % then simply

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
        plot(p1_tri_polyshapes{i})
    end


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
        plot(p1_tri_polyshapes{i})
    end
end
