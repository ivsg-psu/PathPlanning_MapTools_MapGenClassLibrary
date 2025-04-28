% script_test_fcn_MapGen_polytopesPredictLengthCostRatio
% Test function: fcn_MapGen_polytopesPredictLengthCostRatio

% REVISION HISTORY:
% 2025_04_28
% -- first written by S. Harnett


Halton_range = [1 10]; % range of Halton points to use to generate the tiling
tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1]);
des_rad = 0.1; sigma_radius = 0.02; min_rad = 0.001;
[shrunk_field,mu_final,sigma_final] = fcn_MapGen_polytopesShrinkToRadius(tiled_polytopes,des_rad,sigma_radius,min_rad);
field_stats = fcn_MapGen_polytopesStatistics(shrunk_field);
field_avg_r_D = field_stats.avg_r_D;
field_stats_pre_shrink = fcn_MapGen_polytopesStatistics(tiled_polytopes);
field_avg_r_D_pre_shrink = field_stats_pre_shrink.avg_r_D;
R_bar_initial = field_stats_pre_shrink.average_max_radius;
avg_max_rad = field_stats.average_max_radius;
shrunk_distance = R_bar_initial - avg_max_rad;
gap_size = 2*shrunk_distance;
% travel direction is assumed to be left ot right, 0 degrees from horizontal
travel_direction = [1,0];
% traversal crosses the entire field from left to right, going unit length
L_E = 1;
% prediction may fail but should just result in a discarded data point
[field_small_choice_angles,field_big_choice_angles,r_lc_estimates] = ...
    fcn_MapGen_polytopesPredictLengthCostRatio(tiled_polytopes,shrunk_field,gap_size,travel_direction,L_E)
assert(true);
% expected output:
%              r_lc_sparse_worst: 1.0258
%            r_lc_sparse_average: 1.0098
%                r_lc_sparse_std: 0.0133
%          r_lc_sparse_worst_new: 0.8211
%        r_lc_sparse_average_new: 1.0085
%            r_lc_sparse_std_new: 0.0116
%       r_lc_sparse_worst_actual: 0.7725
%     r_lc_sparse_average_actual: 1.0010
%         r_lc_sparse_std_actual: 0.0014
%       r_lc_sparse_worst_linear: 0.7421
%     r_lc_sparse_average_linear: 1.0011
%         r_lc_sparse_std_linear: 0.0016
