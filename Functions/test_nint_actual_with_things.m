close all; clear all; clc;
N_int_actual_all = [];
rd_all= [];
tiles_all = [];
rad_all = [];
std_all = [];
for tiles=25:25:500
    Halton_range = [1 tiles]; % range of Halton points to use to generate the tiling
    tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1]);%,fig_num);
    field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
    N_int_actual = field_stats.linear_density_mean;
    rd = field_stats.avg_r_D;
    rad = field_stats.average_max_radius;
    std = field_stats.std_max_radius;
    N_int_actual_all = [N_int_actual_all, N_int_actual];
    tiles_all = [tiles_all, tiles];
    rd_all = [rd_all, rd];
    rad_all = [rad_all, rad];
    std_all = [std_all, std];
end
figure
hold on
subplot(2,2,1)
plot(tiles_all,N_int_actual_all,'go');
title('Estimated Number of Encountered Obstalces')
xlabel('Number of tiles')
subplot(2,2,2)
plot(rd_all,N_int_actual_all,'ro');
title('Estimated Number of Encountered Obstalces')
xlabel('Mapped Departure Ratio [r_D]')
subplot(2,2,3)
plot(rad_all,N_int_actual_all,'bo');
title('Estimated Number of Encountered Obstalces')
xlabel('Average max radius [km]')
subplot(2,2,4)
plot(std_all,N_int_actual_all./tiles_all,'ko');
title('Estimated Number of Encountered Obstalces normalized by tile count')
xlabel('Std deviation of max radius')
