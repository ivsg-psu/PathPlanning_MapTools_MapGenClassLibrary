% script_test_fcn_MapGen_polytopesPredictLengthCostRatio
% Tests function: fcn_MapGen_polytopesPredictLengthCostRatio

% REVISION HISTORY:
% 2021_10_22
% -- first written by S. Harnett

% clear all; close all;clc;
fig_num = 1;

do_single_test = false;
if do_single_test
    Halton_range = [1 20];
    tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1]);%,fig_num);
    % field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
    % gap_size = 0;
    % [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse] = fcn_MapGen_polytopesPredictLengthCostRatio(tiled_polytopes,gap_size)
    field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
    radii_goals = field_stats.average_max_radius*.40
    des_rad = radii_goals; sigma_radius = 0; min_rad = 0.001;
    [shrunk_field,mu_final,sigma_final] = fcn_MapGen_polytopesShrinkToRadius(tiled_polytopes,des_rad,sigma_radius,min_rad);%,fig_num);
    field_stats = fcn_MapGen_polytopesStatistics(shrunk_field);
    gap_size = field_stats.average_gap_size_G_bar;
    r_D = field_stats.avg_r_D;
    [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse_worst,r_lc_sparse_average,r_lc_sparse_std] = fcn_MapGen_polytopesPredictLengthCostRatio(shrunk_field,gap_size)
end

do_range_test = true;
if do_range_test
    % begin r_D range generation
    r_D = [];
    R_bar_initials = [];
    R_bar_finals = [];
    r_D_theoretical = [];
    r_lc_max_all = [];
    r_lc_avg_all = [];
    r_lc_iterative_all = [];
    r_lc_max_effective_all = [];
    r_lc_avg_effective_all = [];
    r_lc_iterative_effective_all = [];
    r_lc_sparse_worst_all = [];
    r_lc_sparse_average_all = [];
    r_lc_sparse_std_all = [];
    r_lc_sparse_worst_all_new = [];
    r_lc_sparse_average_all_new = [];
    r_lc_sparse_std_all_new = [];
    r_lc_sparse_worst_all_actual = [];
    r_lc_sparse_average_all_actual = [];
    r_lc_sparse_std_all_actual = [];
    shrink_distance = [];
    tiles_failed = [];
    N_int_from_density_all = [];
    N_int_from_shrink_dist_all = [];
    size_percent_failed = [];
    shrunk_distances = [];
    linear_unocc = [];
    debug_ang_all = [];
    debug_sine_all = [];
    N_int_actual_all = [];
    N_int_actual_std_all = [];
    for tiles=100%25:25:25%25:25:125%10:80:500
        Halton_range = [1 tiles]; % range of Halton points to use to generate the tiling
        % tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1]);%,fig_num);
        % title('Halton set');
        % fig_num = fig_num+1;
        % % find r_D for this field
        % field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
        % field_avg_r_D = field_stats.avg_r_D;
        % r_D = [r_D, field_avg_r_D];
        % gap_size = 0;
        % [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse_worst,r_lc_sparse_average,r_lc_sparse_std] = fcn_MapGen_polytopesPredictLengthCostRatio(tiled_polytopes,gap_size)
        % r_lc_max_all = [r_lc_max_all, r_lc_max];
        % r_lc_max_effective_all = [r_lc_max_effective_all, r_lc_max_effective];
        % r_lc_avg_all = [r_lc_avg_all, r_lc_avg];
        % r_lc_avg_effective_all = [r_lc_avg_effective_all, r_lc_avg_effective];
        % r_lc_iterative_all = [r_lc_iterative_all, r_lc_iterative];
        % r_lc_iterative_effective_all = [r_lc_iterative_effective_all, r_lc_iterative_effective];
        % r_lc_sparse_worst_all = [r_lc_sparse_worst_all, r_lc_sparse_worst];
        % r_lc_sparse_average_all = [r_lc_sparse_average_all, r_lc_sparse_average];
        % r_lc_sparse_std_all = [r_lc_sparse_std_all, r_lc_sparse_std];
        % shrink_distance = [shrink_distance, 0];
        % for radii_goals=0.25%0.02:0.02:0.1%0.001:0.010:0.1
        sd_radius_values = [0, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32];
        for sd_radius_index = 1:1:length(sd_radius_values)
            sd_radius = sd_radius_values(sd_radius_index);
            for radii_goals = 0.001:0.005:0.081
                % TODO(@sjharnett) functionalize code at each SD
                % TODO(@sjharnett) loop through, assigning function to SD distinct arrays for plotting
                % TODO(@sjharnett) return means and raw data
                % record multiple points at each departure ratio by generating
                % several maps for this halton tiling and radius
                r_lc_sparse_worst_this_map = [];
                r_lc_sparse_average_this_map = [];
                r_lc_sparse_std_this_map = [];
                r_lc_sparse_worst_this_map_new = [];
                r_lc_sparse_average_this_map_new = [];
                r_lc_sparse_std_this_map_new = [];
                r_lc_sparse_worst_this_map_actual = [];
                r_lc_sparse_average_this_map_actual = [];
                r_lc_sparse_std_this_map_actual = [];
                for i = 1:1:1
                    tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1]);%,fig_num);
                    des_rad = radii_goals; sigma_radius = sd_radius; min_rad = 0.001;
                    % TODO switch this to side shrinking to get gap distance as an output so it can be given to predictor as input
                    [shrunk_field,mu_final,sigma_final] = fcn_MapGen_polytopesShrinkToRadius(tiled_polytopes,des_rad,sigma_radius,min_rad);%,fig_num);
                    field_stats = fcn_MapGen_polytopesStatistics(shrunk_field);
                    gap_size = field_stats.average_gap_size_G_bar;
                    field_avg_r_D = field_stats.avg_r_D;
                    field_stats_pre_shrink = fcn_MapGen_polytopesStatistics(tiled_polytopes);
                    field_avg_r_D_pre_shrink = field_stats_pre_shrink.avg_r_D;
                    % avg_max_rad = field_stats.average_max_radius;
                    shrunk_distance = field_stats_pre_shrink.average_max_radius - field_stats.average_max_radius;
                    shrink_ang = field_stats_pre_shrink.average_vertex_angle;
                    [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse_worst,r_lc_sparse_average,r_lc_sparse_std,r_lc_sparse_worst_new,r_lc_sparse_average_new,r_lc_sparse_std_new,r_lc_sparse_worst_actual,r_lc_sparse_average_actual,r_lc_sparse_std_actual] = fcn_MapGen_polytopesPredictLengthCostRatio(shrunk_field,gap_size,shrunk_distance,shrink_ang)
                    % r_lc_max_all = [r_lc_max_all, r_lc_max];
                    % r_lc_max_effective_all = [r_lc_max_effective_all, r_lc_max_effective];
                    % r_lc_avg_all = [r_lc_avg_all, r_lc_avg];
                    % r_lc_avg_effective_all = [r_lc_avg_effective_all, r_lc_avg_effective];
                    % r_lc_iterative_all = [r_lc_iterative_all, r_lc_iterative];
                    % r_lc_iterative_effective_all = [r_lc_iterative_effective_all, r_lc_iterative_effective];
                    r_lc_sparse_worst_this_map = [r_lc_sparse_worst_this_map, r_lc_sparse_worst];
                    r_lc_sparse_average_this_map = [r_lc_sparse_average_this_map, r_lc_sparse_average];
                    r_lc_sparse_std_this_map = [r_lc_sparse_std_this_map, r_lc_sparse_std];
                    r_lc_sparse_worst_this_map_new = [r_lc_sparse_worst_this_map_new, r_lc_sparse_worst_new];
                    r_lc_sparse_average_this_map_new = [r_lc_sparse_average_this_map_new, r_lc_sparse_average_new];
                    r_lc_sparse_std_this_map_new = [r_lc_sparse_std_this_map_new, r_lc_sparse_std_new];
                    r_lc_sparse_worst_this_map_actual = [r_lc_sparse_worst_this_map_actual, r_lc_sparse_worst_actual];
                    r_lc_sparse_average_this_map_actual = [r_lc_sparse_average_this_map_actual, r_lc_sparse_average_actual];
                    r_lc_sparse_std_this_map_actual = [r_lc_sparse_std_this_map_actual, r_lc_sparse_std_actual];
                    r_D = [r_D, field_avg_r_D];
                    R_bar_initials = [R_bar_initials, field_stats_pre_shrink.average_max_radius];
                    R_bar_finals = [R_bar_finals, field_stats.average_max_radius];
                    N_int_from_density = field_stats.linear_density;
                    N_int_from_density_all = [N_int_from_density_all, N_int_from_density];
                    linear_unocc = [linear_unocc, sqrt(field_stats.unoccupancy_ratio)];
                    shrunk_distances = [shrunk_distances, shrunk_distance];
                    N_int_from_shrink_dist = (sqrt(field_stats.unoccupancy_ratio)*1)/(2*shrunk_distance*sind(field_stats_pre_shrink.average_vertex_angle/2));
                    debug_ang = field_stats_pre_shrink.average_vertex_angle;
                    debug_sine = sind(debug_ang/2);
                    debug_ang_all = [debug_ang_all, debug_ang];
                    debug_sine_all = [debug_sine_all, debug_sine];
                    N_int_from_shrink_dist_all = [N_int_from_shrink_dist_all, N_int_from_shrink_dist];
                    N_int_actual = field_stats.linear_density_mean;
                    N_int_actual_std = field_stats.linear_density_std;
                    N_int_actual_all = [N_int_actual_all, N_int_actual];
                    N_int_actual_std_all = [N_int_actual_std_all, N_int_actual_std];
                end
                r_lc_sparse_worst = mean(r_lc_sparse_worst_this_map);
                r_lc_sparse_average = mean(r_lc_sparse_average_this_map);
                r_lc_sparse_std = mean(r_lc_sparse_std_this_map);
                r_lc_sparse_worst_new = mean(r_lc_sparse_worst_this_map_new);
                r_lc_sparse_average_new = mean(r_lc_sparse_average_this_map_new);
                r_lc_sparse_std_new = mean(r_lc_sparse_std_this_map_new);
                r_lc_sparse_worst_all = [r_lc_sparse_worst_all, r_lc_sparse_worst];
                r_lc_sparse_average_all = [r_lc_sparse_average_all, r_lc_sparse_average];
                r_lc_sparse_std_all = [r_lc_sparse_std_all, r_lc_sparse_std];
                r_lc_sparse_worst_all_new = [r_lc_sparse_worst_all_new, r_lc_sparse_worst_new];
                r_lc_sparse_average_all_new = [r_lc_sparse_average_all_new, r_lc_sparse_average_new];
                r_lc_sparse_std_all_new = [r_lc_sparse_std_all_new, r_lc_sparse_std_new];
                r_lc_sparse_worst_all_actual = [r_lc_sparse_worst_all_actual, r_lc_sparse_worst_actual];
                r_lc_sparse_average_all_actual = [r_lc_sparse_average_all_actual, r_lc_sparse_average_actual];
                r_lc_sparse_std_all_actual = [r_lc_sparse_std_all_actual, r_lc_sparse_std_actual];
                r_D_theoretical = [r_D_theoretical, sqrt(tiles)*des_rad];
            end

%             r_lc_sparse_worst_all = [];
%             r_lc_sparse_average_all = [];
%             r_lc_sparse_std_all = [];
%             r_D = [];
        end
    end
   fprintf('Obstacle fields that could not be predicted:\n')
   fprintf('Num. tiles | Size percent\n')
   fprintf('%d         | %d\n',tiles_failed,size_percent_failed)
end
plot_flag = true; if plot_flag
    figure(607)
    hold on;
    plot(r_D, N_int_from_density_all, 'bo');
    plot(r_D, N_int_from_shrink_dist_all, 'ro');
    % errorbar(r_D,N_int_actual_all,2*N_int_actual_std_all,'go')
    plot(r_D,N_int_actual_all,'go');
    plot(r_D, N_int_from_density_all.*(1-shrunk_distances./R_bar_initials),'mo');
    title('Estimated Number of Encountered Obstalces')
    legend('estimated from linear point density','estimated from occupancy ratio','estimated from ray cast','estimated from percent shrunk')
    xlabel('Mapped Departure Ratio [r_D]')
    ylabel('Estiamted number of encountered obstacles [N_{int}]')
    figure(606)
    hold on;
    plot(r_D, linear_unocc, 'bo');
    plot(r_D, shrunk_distances, 'ro');
    legend('linear unocc','shrunk distances')
    figure(660)
    hold on;
    plot(r_D, shrunk_distances, 'ko');
    plot(r_D, R_bar_initials, 'co');
    plot(r_D, R_bar_finals, 'mo');
    legend('d_{shrink}','R_{bar, initial}','R_{bar, final}')
    figure(608)
    plot(r_D,'ko')
    title('r_D distribution')
    figure(1)
    title('Comparison of Length Cost Ratio versus Departure Ratio from Path Planning and from Geometric Estimate')
    % xlabel('r_D from average radius')
    % ylabel('r_{LC} from side length and vertex angle')
    hold on
    try
        % plot(r_D,r_lc_max_all,'ro')
        % plot(r_D,r_lc_avg_all,'bo')
        % plot(r_D,r_lc_iterative_all,'go')
        % plot(r_D,r_lc_max_effective_all,'rx')
        % plot(r_D,r_lc_avg_effective_all,'bx')
        % plot(r_D,r_lc_iterative_effective_all,'gx')
        plot(r_D_theoretical,r_lc_sparse_worst_all,'md')
        plot(r_D_theoretical,r_lc_sparse_worst_all_new,'gd')
        plot(r_D_theoretical,r_lc_sparse_worst_all_actual,'rd')
        % plot(r_D,r_lc_sparse_average_all,'cd')
        % positive and negative errorbars
        % errorbar(r_D,r_lc_sparse_average_all,2*r_lc_sparse_std_all,'cd')
        % positive only errorbars
        % errorbar(r_D_theoretical,r_lc_sparse_average_all,zeros(1,length(r_lc_sparse_average_all)),2*r_lc_sparse_std_all,zeros(1,length(r_lc_sparse_average_all)),zeros(1,length(r_lc_sparse_average_all)),'cd')
        % errorbar(r_D_theoretical,r_lc_sparse_average_all_new,zeros(1,length(r_lc_sparse_average_all_new)),2*r_lc_sparse_std_all_new,zeros(1,length(r_lc_sparse_average_all_new)),zeros(1,length(r_lc_sparse_average_all_new)),'bd')
        % errorbar(r_D_theoretical,r_lc_sparse_average_all_actual,zeros(1,length(r_lc_sparse_average_all_actual)),2*r_lc_sparse_std_all_actual,zeros(1,length(r_lc_sparse_average_all_actual)),zeros(1,length(r_lc_sparse_average_all_actual)),'yd')
        plot(r_D_theoretical,r_lc_sparse_average_all,'cd')
        plot(r_D_theoretical,r_lc_sparse_average_all_new,'bd')
        plot(r_D_theoretical,r_lc_sparse_average_all_actual,'yd')
        legend('0','0.01','0.02','0.04','0.08','0.16','0.32','worst case original','worst case new','worst case actual','average original','average new','average actual');
        legend('worst case original','worst case new','worst case actual','average original','average new','average actual');
    end
%     x1 = linspace(0,0.65,300);
%     x2=linspace(0.65,.78,100);
%     k1 = 0.4124*x1+41.91*x1.^2;
%     k2 = 0.4124*0.65+41.91*0.65^2-120.3*(x2-0.65)+17.47*(x2-0.65).^2;
%     t1 = 0.0048*x1-0.0016*x1.^2;
%     t2 = 0.0048*0.68-0.0016*0.65^2-0.0009/((0.8118-0.65)^1.25)+0.0009./((0.8118-x2).^1.25);
%     mean1 = k1.*t1 + 1;
%     mode1 = (k1-1).*t1+1;
%     mean2 = k2.*t2 + 1;
%     mode2 = (k2-1).*t2+1;
%     var1 = k1.*t1.^2;
%     var2 = k2.*t2.^2;
%     sd1 = (var1).^2;
%     sd2 = (var2).^2;
%     plot(x1,mean1,'k-')
%     plot(x2,mean2,'k-')
%     bot1 = mean1 - var1;
%     bot2 = mean2 - var2;
%     top1 = mean1 + var1;
%     top2 = mean2 + var2;
%     plot(x1,bot1,'g-');
%     plot(x2,bot2,'g-');
%     plot(x1,top1,'g-');
%     plot(x2,top2,'g-');
%     patch([x1 fliplr(x1)], [bot1 fliplr(top1)], 'g','FaceAlpha',0.2);
%     patch([x2 fliplr(x2)], [bot2 fliplr(top2)], 'g','FaceAlpha',0.2);
%    legend('theoretical max',...
%         'average from side and angle',...
%         'iterative from side and angle',...
%         'effective theoretical max',...
%         'average from effective side and angle',...
%         'iterative from effective side and angle',...
%         'sparse formula, worst case',...
%         'sparse formula, average case',...
%         'gamma distribution curve fit',...
%         'gamma distribution +/- one SD');
end
