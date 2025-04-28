% script_test_fcn_MapGen_polytopesPredictLengthCostStraight
% Test function: fcn_MapGen_polytopesPredictLengthCostRatioStraightPath

% REVISION HISTORY:
% 2025_04_28
% -- first written by S. Harnett
rng(1)
low_pt = 1; high_pt = 20; % range of Halton points to use to generate the tiling
trim_polytopes = fcn_MapGen_haltonVoronoiTiling([low_pt,high_pt],[1 1]);
gap_size = 0.05; % desired average maximum radius
shrunk_polytopes = fcn_MapGen_polytopesShrinkFromEdges(trim_polytopes,gap_size);
A.x = 0; A.y = 0.5; B.x = 1; B.y = 0.5;
r_lc_straight_through = ...
    fcn_MapGen_polytopesPredictLengthCostRatioStraightPath(trim_polytopes,shrunk_polytopes,gap_size,A.x,A.y,B.x,B.y);
assert(round(r_lc_straight_through,4) == 1.3135)
