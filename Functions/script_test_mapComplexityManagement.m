close all; clear all; clc

addpath([pwd '\..\..\Errata_Tutorials_DebugTools\Functions'])

flag_do_plot = 1;
n = 25;
m = 25;
sz = [n m]; % size of board
surface_grid = rand(n,m);
cells = [];
figure(1); hold on
xlabel('x [km]')
ylabel('y [km]')
for ind = 1:1:n
    yvec = [(ind-1)/n (ind-1)/n ind/n ind/n];
    for jind = 1:1:m
        xvec = [(jind-1)/m jind/m jind/m (jind-1)/m];
        face_alpha = surface_grid(ind,jind);
        fill(xvec,yvec,'black','FaceAlpha',face_alpha)
        cells(ind*jind).xv = xvec;
        cells(ind*jind).yv = yvec;
        cells(ind*jind).cost = face_alpha;
    end
end
figure(2); hold on
for ind = 1:1:n
    yvec = [(ind-1)/n (ind-1)/n ind/n ind/n];
    for jind = 1:1:m
        xvec = [(jind-1)/m jind/m jind/m (jind-1)/m];
        face_alpha = surface_grid(ind,jind);
        fill(xvec,yvec,'black','FaceAlpha',face_alpha)
        cells(ind*jind).xv = xvec;
        cells(ind*jind).yv = yvec;
        cells(ind*jind).cost = face_alpha;
    end
end
figure(5); hold on
for ind = 1:1:n
    yvec = [(ind-1)/n (ind-1)/n ind/n ind/n];
    for jind = 1:1:m
        xvec = [(jind-1)/m jind/m jind/m (jind-1)/m];
        face_alpha = surface_grid(ind,jind);
        fill(xvec,yvec,'black','FaceAlpha',face_alpha)
        cells(ind*jind).xv = xvec;
        cells(ind*jind).yv = yvec;
        cells(ind*jind).cost = face_alpha;
    end
end
% generate map
Halton_seed = 10;
low_pt = 1+Halton_seed; high_pt = 30+Halton_seed; % range of Halton points to use to generate the tiling
trim_polytopes = fcn_MapGen_haltonVoronoiTiling([low_pt,high_pt],[1 1]);
% shink the polytopes so that they are no longer tiled
gap_size = 0.08; % desired average maximum radius
polytopes = fcn_MapGen_polytopesShrinkFromEdges(trim_polytopes,gap_size);
% plot the map

    fig = 2; % figure to plot on
    line_spec = 'b-'; % edge line plotting
    line_width = 2; % linewidth of the edge
    axes_limits = [0 1 0 1]; % x and y axes limits
    axis_style = 'square'; % plot axes style
    fcn_MapGen_plotPolytopes(polytopes,fig,line_spec,line_width,axes_limits,axis_style);
    hold on
    box on
    xlabel('x [km]')
    ylabel('y [km]')
    poly_map_stats = fcn_MapGen_polytopesStatistics(polytopes);
    my_title = sprintf('obstacle departure ratio: %2f',poly_map_stats.avg_r_D);
    title(my_title);
%% create polyshapes array from polytopes array
    polyshapes = [];
    polyshape_background = polyshape([0 0 1 1],[0 1 1 0]);
    if flag_do_plot
        figure(3)
        clf
    end
    for i = 1:length(polytopes)
        new_polyshape = polyshape(polytopes(i).xv,polytopes(i).yv);
        polyshapes = [polyshapes, new_polyshape];
        if flag_do_plot
            figure(3)
            hold on
            title('Polyshapes before subtraction')
            plot(polyshapes(i))
        end
        polyshape_background = subtract(polyshape_background,new_polyshape);
    end
figure(3)
plot(polyshape_background)
figure(4)
plot(polyshape_background)

background_verts = polyshape_background.Vertices;
background_verts(sum(isnan(background_verts), 2) >= 1, :) = [];
DT = delaunayTriangulation(background_verts);
figure(6)
triplot(DT)
hold on; box on;
plot(polyshape_background)

[p_tri_polyshapes, p_tri_polytopes] = INTERNAL_fcn_triangulatePolyshape(polyshape_background,flag_do_plot)
figure(5)
for i = 1:length(p_tri_polyshapes)
    plot(p_tri_polyshapes{i})
end
% make a grid Thrusday
% convert some grid cells to polytopes Friday-Monday
    % this means find departure ratio
    % find cricical traversal cost
    % note all cell vertices above critical cost
    % some polys are inside, some outside, some on edge
    % find overlapping polytopes
    % remove them from the list replacing them with their union()
    % non overlapping squares just get added to the list

% subtract polytopes from background and triangulate Monday-Tuesday
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
