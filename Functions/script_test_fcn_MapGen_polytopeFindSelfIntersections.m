% script_test_fcn_MapGen_polytopeFindSelfIntersections
% Tests function: fcn_MapGen_polytopeFindSelfIntersections

% REVISION HISTORY:
% 2021_08_02
% -- first written by S. Brennan


%% Basic example of self-intersection
fig_num = 1;
vertices = [0 0; 1 0; 0.5 1.5; 1 1; 0 1; 0 0];
fcn_MapGen_polytopeFindSelfIntersections(...
    vertices,fig_num);

