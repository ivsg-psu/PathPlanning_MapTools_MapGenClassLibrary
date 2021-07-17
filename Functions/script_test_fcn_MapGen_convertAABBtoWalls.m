% script_test_fcn_MapGen_convertAABBtoWalls
% Tests: fcn_MapGen_convertAABBtoWalls

% 
% REVISION HISTORY:
% 
% 2021_07_17 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§


%% simple example
fig_num = 1;
AABB = [0 0 1 1];
walls = fcn_MapGen_convertAABBtoWalls(AABB,fig_num);

