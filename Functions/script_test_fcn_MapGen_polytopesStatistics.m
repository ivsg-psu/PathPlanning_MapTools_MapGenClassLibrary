% script_test_fcn_MapGen_polytopesStatistics
% Tests: fcn_MapGen_polytopesStatistics

% 
% REVISION HISTORY:
% 
% 2021_07_12 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§

close all;

% Generate a set of polytopes from the Halton set
fig_num = 12;
Halton_range = [1 1000]; % range of Halton points to use to generate the tiling
polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1]);

fcn_MapGen_polytopesStatistics(...
    polytopes,...
    fig_num)



