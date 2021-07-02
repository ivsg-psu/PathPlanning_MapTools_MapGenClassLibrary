% script_test_fcn_MapGen_polytopeMapGen
% Tests function: fcn_MapGen_polytopeMapGen

% REVISION HISTORY:
% 2021_06_06 
% -- first written by S. Brennan. 

close all;
fig_num = 1;
stretch = [1 1];
polytopes = fcn_MapGen_haltonVoronoiTiling([1 1000],stretch,fig_num);

%% Show that the stretch works
fig_num = 200;
stretch = [100 200];
polytopes = fcn_MapGen_haltonVoronoiTiling([1 1000],stretch,fig_num);

