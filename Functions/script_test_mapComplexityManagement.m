close all; clear all; clc

addpath([pwd '\..\..\Errata_Tutorials_DebugTools\Functions'])


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
