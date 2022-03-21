% script_fcn_MapGen_polytopeShrinkFromEdges
% Tests function: fcn_MapGen_polytopeShrinkFromEdges

% REVISION HISTORY:
% 2022_01_17
% -- first written by S. Harnett using
% script_test_fcn_MapGen_polytopesShrinkToRadius as a template

%% Set up variables
close all;
polytopes = fcn_MapGen_haltonVoronoiTiling([1 100]);

fig_num = 1;
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);

%% Basic example of uniform shrinking
fig_num = 2;
des_gap_size = 0.05;
shrunk_polytopes1=...
    fcn_MapGen_polytopesShrinkFromEdges(...
    trim_polytopes,des_gap_size,fig_num);
field_stats = fcn_MapGen_polytopesStatistics(shrunk_polytopes1);
r_occ_meas = field_stats.occupancy_ratio; % calculated occupancy ratio
r_unocc_meas = field_stats.unoccupancy_ratio;
G_bar = field_stats.average_gap_size_G_bar;
rho = field_stats.linear_density_mean;
estimated_unoccupancy_ratio = G_bar^2*rho; % theoretial occupancy ratio from gap size
assert(isequal(round(G_bar,4),round(des_gap_size,4)));
assert(isequal(round(r_unocc_meas,4),round(estimated_unoccupancy_ratio,4)));
