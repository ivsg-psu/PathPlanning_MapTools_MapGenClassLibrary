% script_test_fcn_MapGen_generatePolysFromTiling
% Tests: fcn_MapGen_generatePolysFromTiling

% 
% REVISION HISTORY:
% 
% 2021_07_02 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§


%% pull halton set
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values

%% pick values from halton set
Halton_range = [3901        4001];
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
seed_points = points_scrambled(low_pt:high_pt,:);
[V,C] = voronoin(seed_points);
% V = V.*stretch;

%% fill polytopes from tiling
fig_num = 1;
AABB = [0 0 1 1];
stretch = [1 1];
polytopes = fcn_MapGen_generatePolysFromTiling(seed_points,V,C,AABB, stretch,fig_num);

