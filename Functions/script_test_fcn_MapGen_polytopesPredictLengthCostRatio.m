% script_test_fcn_MapGen_polytopesPredictLengthCostRatio
% Tests function: fcn_MapGen_polytopesPredictLengthCostRatio

% REVISION HISTORY:
% 2021_10_22
% -- first written by S. Harnett

clear all; close all;clc;
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
    r_lc_max_all = [];
    r_lc_avg_all = [];
    r_lc_iterative_all = [];
    r_lc_max_effective_all = [];
    r_lc_avg_effective_all = [];
    r_lc_iterative_effective_all = [];
    r_lc_sparse_worst_all = [];
    r_lc_sparse_average_all = [];
    r_lc_sparse_std_all= [];
    shrink_distance = [];
    tiles_failed = [];
    size_percent_failed = [];
    for tiles=100%25:25:25%25:25:125%10:80:500
        Halton_range = [1 tiles]; % range of Halton points to use to generate the tiling
        tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1]);%,fig_num);
        title('Halton set');
        fig_num = fig_num+1;
        % find r_D for this field
        field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
        field_avg_r_D = field_stats.avg_r_D;
        r_D = [r_D, field_avg_r_D];
        gap_size = 0;
        [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse_worst,r_lc_sparse_average,r_lc_sparse_std] = fcn_MapGen_polytopesPredictLengthCostRatio(tiled_polytopes,gap_size)
        r_lc_max_all = [r_lc_max_all, r_lc_max];
        r_lc_max_effective_all = [r_lc_max_effective_all, r_lc_max_effective];
        r_lc_avg_all = [r_lc_avg_all, r_lc_avg];
        r_lc_avg_effective_all = [r_lc_avg_effective_all, r_lc_avg_effective];
        r_lc_iterative_all = [r_lc_iterative_all, r_lc_iterative];
        r_lc_iterative_effective_all = [r_lc_iterative_effective_all, r_lc_iterative_effective];
        r_lc_sparse_worst_all = [r_lc_sparse_worst_all, r_lc_sparse_worst];
        r_lc_sparse_average_all = [r_lc_sparse_average_all, r_lc_sparse_average];
        r_lc_sparse_std_all = [r_lc_sparse_std_all, r_lc_sparse_std];
        shrink_distance = [shrink_distance, 0];
        % for radii_goals=0.25%0.02:0.02:0.1%0.001:0.010:0.1
        for radii_goals = 0.001:0.005:0.081
            for sd_radius = 0:0.005:0.032
                des_rad = radii_goals; sigma_radius = sd_radius; min_rad = 0.001;
                % TODO switch this to side shrinking to get gap distance as an output so it can be given to predictor as input
                [shrunk_field,mu_final,sigma_final] = fcn_MapGen_polytopesShrinkFromEdges(tiled_polytopes,des_rad,sigma_radius,min_rad);%,fig_num);
                field_stats = fcn_MapGen_polytopesStatistics(shrunk_field);
                gap_size = field_stats.average_gap_size_G_bar;
                field_avg_r_D = field_stats.avg_r_D;
                % avg_max_rad = field_stats.average_max_radius;
                shrink_distance = [shrink_distance, gap_size];%avg_max_rad-des_rad];
                try
                   [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse_worst,r_lc_sparse_average,r_lc_sparse_std] = fcn_MapGen_polytopesPredictLengthCostRatio(shrunk_field,gap_size)
                catch
                    tiles_failed = [tiles_failed, tiles];
                    size_percent_failed = [size_percent_failed, size_percent];
                end
                r_lc_max_all = [r_lc_max_all, r_lc_max];
                r_lc_max_effective_all = [r_lc_max_effective_all, r_lc_max_effective];
                r_lc_avg_all = [r_lc_avg_all, r_lc_avg];
                r_lc_avg_effective_all = [r_lc_avg_effective_all, r_lc_avg_effective];
                r_lc_iterative_all = [r_lc_iterative_all, r_lc_iterative];
                r_lc_iterative_effective_all = [r_lc_iterative_effective_all, r_lc_iterative_effective];
                r_lc_sparse_worst_all = [r_lc_sparse_worst_all, r_lc_sparse_worst];
                r_lc_sparse_average_all = [r_lc_sparse_average_all, r_lc_sparse_average];
                r_lc_sparse_std_all = [r_lc_sparse_std_all, r_lc_sparse_std];
                r_D = [r_D, field_avg_r_D];
                close all;
           end
        end
    end%     close all
   fprintf('Obstacle fields that could not be predicted:\n')
   fprintf('Num. tiles | Size percent\n')
   fprintf('%d        | %d\n',tiles_failed,size_percent_failed)
end
plot_flag = true;
if plot_flag
    figure(608)
    plot(r_D,'ko')
    title('r_D distribution')
    figure(609)
    title('r_D versus r_{LC}')
    xlabel('r_D from average radius')
    ylabel('r_{LC} from side length and vertex angle')
    hold on
    try
        plot(r_D,r_lc_max_all,'ro')
        plot(r_D,r_lc_avg_all,'bo')
        plot(r_D,r_lc_iterative_all,'go')
        plot(r_D,r_lc_max_effective_all,'rx')
        plot(r_D,r_lc_avg_effective_all,'bx')
        plot(r_D,r_lc_iterative_effective_all,'gx')
        plot(r_D,r_lc_sparse_worst_all,'md')
        % plot(r_D,r_lc_sparse_average_all,'cd')
        % errorbar(r_D,r_lc_sparse_average_all,2*r_lc_sparse_std_all,'cd')
        % positive only errorbars
        errorbar(r_D,r_lc_sparse_average_all,zeros(1,length(r_lc_sparse_average_all)),2*r_lc_sparse_std_all,zeros(1,length(r_lc_sparse_average_all)),zeros(1,length(r_lc_sparse_average_all)),'cd')
    end
    x1 = linspace(0,0.65,300);
    x2=linspace(0.65,.78,100);
    k1 = 0.4124*x1+41.91*x1.^2;
    k2 = 0.4124*0.65+41.91*0.65^2-120.3*(x2-0.65)+17.47*(x2-0.65).^2;
    t1 = 0.0048*x1-0.0016*x1.^2;
    t2 = 0.0048*0.68-0.0016*0.65^2-0.0009/((0.8118-0.65)^1.25)+0.0009./((0.8118-x2).^1.25);
    mean1 = k1.*t1 + 1;
    mode1 = (k1-1).*t1+1;
    mean2 = k2.*t2 + 1;
    mode2 = (k2-1).*t2+1;
    var1 = k1.*t1.^2;
    var2 = k2.*t2.^2;
    sd1 = (var1).^2;
    sd2 = (var2).^2;
    plot(x1,mean1,'k-')
    plot(x2,mean2,'k-')
    bot1 = mean1 - var1;
    bot2 = mean2 - var2;
    top1 = mean1 + var1;
    top2 = mean2 + var2;
    plot(x1,bot1,'g-');
    plot(x2,bot2,'g-');
    plot(x1,top1,'g-');
    plot(x2,top2,'g-');
    patch([x1 fliplr(x1)], [bot1 fliplr(top1)], 'g','FaceAlpha',0.2);
    patch([x2 fliplr(x2)], [bot2 fliplr(top2)], 'g','FaceAlpha',0.2);
    legend('theoretical max',...
        'average from side and angle',...
        'iterative from side and angle',...
        'effective theoretical max',...
        'average from effective side and angle',...
        'iterative from effective side and angle',...
        'sparse formula, worst case',...
        'sparse formula, average case',...
        'gamma distribution curve fit',...
        'gamma distribution +/- one SD');
    figure
end
