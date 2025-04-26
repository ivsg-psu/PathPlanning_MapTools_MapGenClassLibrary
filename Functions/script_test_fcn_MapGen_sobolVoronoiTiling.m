% script_test_fcn_MapGen_sobolVoronoiTiling
% Tests function: fcn_MapGen_sobolVoronoiTiling

% REVISION HISTORY:
% 2021_06_06 
% -- first written by S. Brennan. 

close all;

%% Basic call example
fig_num = 1;
figure(fig_num);
clf;


polytopes = fcn_MapGen_sobolVoronoiTiling([1 1000],[1 1],fig_num);

assert(isstruct(polytopes));
