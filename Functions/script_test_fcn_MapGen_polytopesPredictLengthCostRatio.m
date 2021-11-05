% script_test_fcn_MapGen_polytopesPredictLengthCostRatio
% Tests function: fcn_MapGen_polytopesPredictLengthCostRatio

% REVISION HISTORY:
% 2021_10_22
% -- first written by S. Harnett

clear all; close all;clc;
fig_num = 1;
% begin r_D range generation
r_D = [];
r_lc_max_all = [];
r_lc_avg_all = [];
r_lc_iterative_all = [];
% shrink_distance = [];
for tiles=50:25:75%10:80:500
    Halton_range = [1 tiles]; % range of Halton points to use to generate the tiling
    tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);
    title('Halton set');
    fig_num = fig_num+1;
    % find r_D for this field
    field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
    field_avg_r_D = field_stats.avg_r_D;
    r_D = [r_D, field_avg_r_D];
    gap_size = 0;
    [r_lc_max,r_lc_avg,r_lc_iterative] = fcn_MapGen_polytopesPredictLengthCostRatio(tiled_polytopes,gap_size)
    r_lc_max_all = [r_lc_max_all, r_lc_max];
    r_lc_avg_all = [r_lc_avg_all, r_lc_avg];
    r_lc_iterative_all = [r_lc_iterative_all, r_lc_iterative];
    % shrink_distance = [shrink_distance, 0];
    for radii_goals=0.02:0.02:0.1%0.001:0.010:0.1
        try
            des_rad = radii_goals; sigma_radius = 0; min_rad = 0.001;
            % TODO switch this to side shrinking to get gap distance as an output so it can be given to predictor as input
            [shrunk_field,mu_final,sigma_final] = fcn_MapGen_polytopesShrinkToRadius(tiled_polytopes,des_rad,sigma_radius,min_rad,fig_num);
            field_stats = fcn_MapGen_polytopesStatistics(shrunk_field);
            gap_size = field_stats.average_gap_size_G_bar;
            field_avg_r_D = field_stats.avg_r_D;
            r_D = [r_D, field_avg_r_D];
            % avg_max_rad = field_stats.average_max_radius;
            % shrink_distance = [shrink_distance, avg_max_rad-des_rad];
            [r_lc_max,r_lc_avg,r_lc_iterative] = fcn_MapGen_polytopesPredictLengthCostRatio(shrunk_field,gap_size)
            r_lc_max_all = [r_lc_max_all, r_lc_max];
            r_lc_avg_all = [r_lc_avg_all, r_lc_avg];
            r_lc_iterative_all = [r_lc_iterative_all, r_lc_iterative];
        end
    end
end
close all
figure
plot(r_D,'ko')
title('r_D distribution')
figure(609)
title('r_D versus r_{LC}')
xlabel('r_D from average radius')
ylabel('r_{LC} from side length and vertex angle')
hold on
plot(r_D,r_lc_max_all,'ro')
plot(r_D,r_lc_avg_all,'bo')
plot(r_D,r_lc_iterative_all,'go')
legend('theoretical max from half of average angle','average from average side length and average angle','measured from subset of polytopes')
% figure
% plot(shrink_distance,r_D,'go')
% end r_D range generation

