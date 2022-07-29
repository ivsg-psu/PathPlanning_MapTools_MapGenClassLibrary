% script_test_fcn_MapGen_polytopesRadiusDistributions
% Tests: fcn_MapGen_polytopesRadiusDistributions

%
% REVISION HISTORY:
%
% 2022_07_28
% -- first write of script
%%%%%%%%%%%%%%§

close all;

% Generate a set of polytopes from the Halton set
Halton_range = [200 1001]; % range of Halton points to use to generate the tiling
polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1]);

poly_size_stats = fcn_MapGen_polytopesRadiusDistributions(polytopes);