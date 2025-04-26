% script_test_fcn_MapGen_randomNormalVoronoiTiling
% Tests function: fcn_MapGen_randomNormalVoronoiTiling

% REVISION HISTORY:
% 2021_06_06
% -- first written by S. Brennan.

close all;

%% Basic call example
fig_num = 1;
figure(fig_num);
clf;

polytopes = fcn_MapGen_randomNormalVoronoiTiling([1 200],[1 1],fig_num);

assert(isstruct(polytopes));
