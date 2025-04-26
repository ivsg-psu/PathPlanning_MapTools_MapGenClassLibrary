% script_test_fcn_MapGen_randVoronoiTiling
% Tests function: fcn_MapGen_randVoronoiTiling

% REVISION HISTORY:
% 2021_06_06
% -- first written by S. Brennan.
close all

%% Example of basic call
fig_num = 1;
figure(fig_num);
clf;

polytopes = fcn_MapGen_randVoronoiTiling([1 1000],[1 1],fig_num);

assert(isstruct(polytopes));
