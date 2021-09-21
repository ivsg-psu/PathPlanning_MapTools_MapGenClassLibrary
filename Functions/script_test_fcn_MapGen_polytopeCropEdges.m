% script_test_fcn_MapGen_polytopeCropEdges
% Tests function: fcn_MapGen_polytopeCropEdges

% REVISION HISTORY:
% 2021_06_06
% -- first written by S. Brennan.

close all;
polytopes = fcn_MapGen_haltonVoronoiTiling([1 1000]);

fig_num = 2;
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);

assert(true);
