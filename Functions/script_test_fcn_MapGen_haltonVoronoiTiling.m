% script_test_fcn_MapGen_polytopeMapGen
% Tests function: fcn_MapGen_polytopeMapGen

% REVISION HISTORY:
% 2021_06_06
% -- first written by S. Brennan.

close all;

%%
fig_num = 1;
figure(fig_num);
clf;

stretch = [1 1];
polytopes = fcn_MapGen_haltonVoronoiTiling([1 1000],stretch,fig_num);

assert(isstruct(polytopes));

%% Show that the stretch works
fig_num = 200;
figure(fig_num);
clf;

stretch = [100 200];
polytopes = fcn_MapGen_haltonVoronoiTiling([1 1000],stretch,fig_num);

assert(isstruct(polytopes));
