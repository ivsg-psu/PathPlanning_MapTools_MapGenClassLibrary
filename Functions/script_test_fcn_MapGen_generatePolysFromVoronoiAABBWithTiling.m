% script_test_fcn_MapGen_generatePolysFromVoronoiAABBWithTiling
% Tests: fcn_MapGen_generatePolysFromVoronoiAABBWithTiling

%
% REVISION HISTORY:
%
% 2023_03_13 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§

close all;


%% fill polytopes from tiling
fig_num = 1;

% pull halton set
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values

% pick values from halton set
Halton_range = [1801 1901];
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
seed_points = points_scrambled(low_pt:high_pt,:);


AABB = [0 0 1 1];
stretch = [1 1];
polytopes = fcn_MapGen_generatePolysFromVoronoiAABBWithTiling(seed_points,AABB, stretch,fig_num);
assert(true);


%% check empty figure number (should not plot)
fig_num = [];

% pull halton set
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values

% pick values from halton set
Halton_range = [1801 1901];
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
seed_points = points_scrambled(low_pt:high_pt,:);

AABB = [0 0 1 1];
stretch = [1 1];
polytopes = fcn_MapGen_generatePolysFromVoronoiAABBWithTiling(seed_points,AABB, stretch,fig_num);
