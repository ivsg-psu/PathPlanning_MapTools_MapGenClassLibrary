close all; clear all; clc;
N_int_actual_all = [];
rd_all= [];
tiles_all = [];
rad_all = [];
std_all = [];
sharpness_all = [];
for tiles=25:25:1000
    Halton_range = [1 tiles]; % range of Halton points to use to generate the tiling
    tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1]);
    field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
    N_int_actual = field_stats.linear_density_mean;
    rd = field_stats.avg_r_D;
    rad = field_stats.average_max_radius;
    std = field_stats.std_max_radius;
    sharpness = field_stats.average_sharpness;
    N_int_actual_all = [N_int_actual_all, N_int_actual];
    tiles_all = [tiles_all, tiles];
    rd_all = [rd_all, rd];
    rad_all = [rad_all, rad];
    std_all = [std_all, std];
    sharpness_all = [sharpness_all, sharpness];
end
% repeat with random set instead of Halton set
N_int_actual_all_rand = [];
rd_all_rand= [];
tiles_all_rand = [];
rad_all_rand = [];
std_all_rand = [];
sharpness_all_rand = [];
for tiles=25:25:1000
    tile_range = [1 tiles];
    tiled_polytopes = fcn_MapGen_randomNormalVoronoiTiling(tile_range,[1 1]);
    field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
    N_int_actual = field_stats.linear_density_mean;
    rd = field_stats.avg_r_D;
    rad = field_stats.average_max_radius;
    std = field_stats.std_max_radius;
    sharpness = field_stats.average_sharpness;
    N_int_actual_all_rand = [N_int_actual_all_rand, N_int_actual];
    tiles_all_rand = [tiles_all_rand, tiles];
    rd_all_rand = [rd_all_rand, rd];
    rad_all_rand = [rad_all_rand, rad];
    std_all_rand = [std_all_rand, std];
    sharpness_all_rand = [sharpness_all_rand, sharpness];
end
% plot
figure(1)
subplot(2,2,1)
hold on
plot(tiles_all,N_int_actual_all,'go');
plot(tiles_all_rand,N_int_actual_all_rand,'mo');
legend('Halton set','Random normal set')
title('Estimated Number of Encountered Obstalces')
xlabel('Number of tiles')
subplot(2,2,2)
hold on
plot(rd_all,N_int_actual_all,'go');
plot(rd_all_rand,N_int_actual_all_rand,'mo');
title('Estimated Number of Encountered Obstalces')
xlabel('Mapped Departure Ratio [r_D]')
subplot(2,2,3)
hold on
plot(rad_all,N_int_actual_all,'go');
plot(rad_all_rand,N_int_actual_all_rand,'mo');
title('Estimated Number of Encountered Obstalces')
xlabel('Average max radius [km]')
subplot(2,2,4)
hold on
plot(std_all,N_int_actual_all./tiles_all,'go');
plot(std_all_rand,N_int_actual_all_rand./tiles_all_rand,'mo');
title('Estimated Number of Encountered Obstalces normalized by tile count')
xlabel('Std deviation of max radius')

figure(2)
hold on
plot(sharpness_all,N_int_actual_all,'go');
plot(sharpness_all_rand,N_int_actual_all_rand,'mo');
legend('Halton set','Random normal set')
title('Estimated Number of Encountered Obstalces')
xlabel('Average Sharpness Ratio for Field')

figure(3)
hold on
plot(std_all,sharpness_all,'go');
plot(std_all_rand,sharpness_all_rand,'mo');
legend('Halton set','Random normal set')
xlabel('Std deviation of max radius')
ylabel('Average Sharpness Ratio for Field')

figure(4)
hold on
plot(tiles_all,sharpness_all,'go');
plot(tiles_all_rand,sharpness_all_rand,'mo');
legend('Halton set','Random normal set')
ylabel('Average Sharpness Ratio for Field')
xlabel('Number of tiles')
