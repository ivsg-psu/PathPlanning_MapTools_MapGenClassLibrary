% script_fcn_MapGen_polytopeShrinkFromEdges
% Tests function: fcn_MapGen_polytopeShrinkFromEdges

% REVISION HISTORY:
% 2025_04_28
% -- first written by S. Harnett using single case from data generation script:
% demo_fcn_MapGen_polytopeShrinkFromEdges


% Set up variables
polytopes = fcn_MapGen_haltonVoronoiTiling([1 20]);
fig_num = 1;
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);

des_gap_size = 0.1;
shrunk_polytopes1=...
    fcn_MapGen_polytopesShrinkFromEdges(...
    trim_polytopes,des_gap_size,fig_num);
