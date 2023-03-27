% script_test_fcn_MapGen_polytopesStatistics
% Tests: fcn_MapGen_polytopesStatistics

%
% REVISION HISTORY:
%
% 2021_07_12 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§

close all;

% Generate a set of polytopes from the Halton set
fig_num = 12;
Halton_range = [200 301]; % range of Halton points to use to generate the tiling
polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1]);

pre_shrink_stats = fcn_MapGen_polytopesStatistics(...
    polytopes,...
    fig_num);


shrinkage = 0.001;
%% Shrink all polytopes by a gap using radial shrinking
shrunk_polytopes_radial = polytopes;
for ith_poly = 1:length(polytopes)
    orig_radius = polytopes(ith_poly).max_radius;
    des_rad = orig_radius - shrinkage;

    tolerance = 1e-5; % This is the edge distance below which vertices are merged together in the polytope
    shrunk_polytopes_radial(ith_poly) =...
        fcn_MapGen_polytopeShrinkToRadius(...
        polytopes(ith_poly),des_rad,tolerance);
end

% Check the statistics on these
fig_num = 13;
fcn_MapGen_polytopesStatistics(...
    shrunk_polytopes_radial,...
    fig_num);

%% Shrink all polytopes by a gap using edge shrinking
shrunk_polytopes_edge = polytopes;
for ith_poly = 1:length(polytopes)
    shrunk_polytopes_edge(ith_poly) = ...
        fcn_MapGen_polytopeShrinkFromEdges(...
        polytopes(ith_poly),shrinkage);
end

% Check the statistics on these
fig_num = 14;
fcn_MapGen_polytopesStatistics(...
    shrunk_polytopes_edge,...
    fig_num);

% Show the original polytopes
h_plot = subplot(2,3,1);
hold on;
line_width = 2;
fcn_MapGen_plotPolytopes(polytopes,h_plot,'r',line_width);
