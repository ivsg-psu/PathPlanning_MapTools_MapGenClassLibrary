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
    tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);
    % field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
    % gap_size = 0;
    % [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse] = fcn_MapGen_polytopesPredictLengthCostRatio(tiled_polytopes,gap_size)
    field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
    radii_goals = field_stats.average_max_radius*.40
    des_rad = radii_goals; sigma_radius = 0; min_rad = 0.001;
    [shrunk_field,mu_final,sigma_final] = fcn_MapGen_polytopesShrinkToRadius(tiled_polytopes,des_rad,sigma_radius,min_rad,fig_num);
    field_stats = fcn_MapGen_polytopesStatistics(shrunk_field);
    gap_size = field_stats.average_gap_size_G_bar;
    r_D = field_stats.avg_r_D;
    [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse_worst,r_lc_sparse_average] = fcn_MapGen_polytopesPredictLengthCostRatio(shrunk_field,gap_size)
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
    shrink_distance = [];
    for tiles=25:25:25%25:25:125%10:80:500
        Halton_range = [1 tiles]; % range of Halton points to use to generate the tiling
        tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);
        title('Halton set');
        fig_num = fig_num+1;
        % find r_D for this field
        field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
        field_avg_r_D = field_stats.avg_r_D;
        r_D = [r_D, field_avg_r_D];
        gap_size = 0;
        [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse_worst,r_lc_sparse_average] = fcn_MapGen_polytopesPredictLengthCostRatio(tiled_polytopes,gap_size)
        r_lc_max_all = [r_lc_max_all, r_lc_max];
        r_lc_max_effective_all = [r_lc_max_effective_all, r_lc_max_effective];
        r_lc_avg_all = [r_lc_avg_all, r_lc_avg];
        r_lc_avg_effective_all = [r_lc_avg_effective_all, r_lc_avg_effective];
        r_lc_iterative_all = [r_lc_iterative_all, r_lc_iterative];
        r_lc_iterative_effective_all = [r_lc_iterative_effective_all, r_lc_iterative_effective];
        r_lc_sparse_worst_all = [r_lc_sparse_worst_all, r_lc_sparse_worst];
        r_lc_sparse_average_all = [r_lc_sparse_average_all, r_lc_sparse_average];
        shrink_distance = [shrink_distance, 0];
        for size_percent = 0.2:0.2:0.8
            radii_goals = field_stats.average_max_radius*size_percent;
             des_rad = radii_goals; sigma_radius = 0; min_rad = 0.001;
             % TODO switch this to side shrinking to get gap distance as an output so it can be given to predictor as input
             [shrunk_field,mu_final,sigma_final] = fcn_MapGen_polytopesShrinkToRadius(tiled_polytopes,des_rad,sigma_radius,min_rad,fig_num);
             field_stats = fcn_MapGen_polytopesStatistics(shrunk_field);
             gap_size = field_stats.average_gap_size_G_bar;
             field_avg_r_D = field_stats.avg_r_D;
             % avg_max_rad = field_stats.average_max_radius;
            shrink_distance = [shrink_distance, gap_size];%avg_max_rad-des_rad];
            [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse_worst,r_lc_sparse_average] = fcn_MapGen_polytopesPredictLengthCostRatio(shrunk_field,gap_size)
            r_lc_max_all = [r_lc_max_all, r_lc_max];
            r_lc_max_effective_all = [r_lc_max_effective_all, r_lc_max_effective];
            r_lc_avg_all = [r_lc_avg_all, r_lc_avg];
            r_lc_avg_effective_all = [r_lc_avg_effective_all, r_lc_avg_effective];
            r_lc_iterative_all = [r_lc_iterative_all, r_lc_iterative];
            r_lc_iterative_effective_all = [r_lc_iterative_effective_all, r_lc_iterative_effective];
            r_lc_sparse_worst_all = [r_lc_sparse_worst_all, r_lc_sparse_worst];
            r_lc_sparse_average_all = [r_lc_sparse_average_all, r_lc_sparse_average];
            r_D = [r_D, field_avg_r_D];
        end
    end
    for tiles=25:25:100%25:25:125%10:80:500
            Halton_range = [1 tiles]; % range of Halton points to use to generate the tiling
            tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);
            title('Halton set');
            fig_num = fig_num+1;
            % find r_D for this field
            field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
             field_avg_r_D = field_stats.avg_r_D;
             r_D = [r_D, field_avg_r_D];
            gap_size = 0;
            [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse_worst,r_lc_sparse_average] = fcn_MapGen_polytopesPredictLengthCostRatio(tiled_polytopes,gap_size)
            r_lc_max_all = [r_lc_max_all, r_lc_max];
            r_lc_max_effective_all = [r_lc_max_effective_all, r_lc_max_effective];
            r_lc_avg_all = [r_lc_avg_all, r_lc_avg];
            r_lc_avg_effective_all = [r_lc_avg_effective_all, r_lc_avg_effective];
            r_lc_iterative_all = [r_lc_iterative_all, r_lc_iterative];
            r_lc_iterative_effective_all = [r_lc_iterative_effective_all, r_lc_iterative_effective];
            r_lc_sparse_worst_all = [r_lc_sparse_worst_all, r_lc_sparse_worst];
            r_lc_sparse_average_all = [r_lc_sparse_average_all, r_lc_sparse_average];
            shrink_distance = [shrink_distance, 0];
            for size_percent = 0.55:0.05:0.85
                radii_goals = field_stats.average_max_radius*size_percent;
                des_rad = radii_goals; sigma_radius = 0; min_rad = 0.001;
                % TODO switch this to side shrinking to get gap distance as an output so it can be given to predictor as input
                [shrunk_field,mu_final,sigma_final] = fcn_MapGen_polytopesShrinkToRadius(tiled_polytopes,des_rad,sigma_radius,min_rad,fig_num);
                field_stats = fcn_MapGen_polytopesStatistics(shrunk_field);
                gap_size = field_stats.average_gap_size_G_bar;
                field_avg_r_D = field_stats.avg_r_D;
                % avg_max_rad = field_stats.average_max_radius;
                shrink_distance = [shrink_distance, gap_size];%avg_max_rad-des_rad];
                [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse_worst,r_lc_sparse_average] = fcn_MapGen_polytopesPredictLengthCostRatio(shrunk_field,gap_size)
                r_lc_max_all = [r_lc_max_all, r_lc_max];
                r_lc_max_effective_all = [r_lc_max_effective_all, r_lc_max_effective];
                r_lc_avg_all = [r_lc_avg_all, r_lc_avg];
                r_lc_avg_effective_all = [r_lc_avg_effective_all, r_lc_avg_effective];
                r_lc_iterative_all = [r_lc_iterative_all, r_lc_iterative];
                r_lc_iterative_effective_all = [r_lc_iterative_effective_all, r_lc_iterative_effective];
                r_lc_sparse_worst_all = [r_lc_sparse_worst_all, r_lc_sparse_worst];
                r_lc_sparse_average_all = [r_lc_sparse_average_all, r_lc_sparse_average];
                r_D = [r_D, field_avg_r_D];
            end
    end
    tiles_failed = [];
    size_percent_failed = [];
    for tiles=25:75:500
            Halton_range = [1 tiles]; % range of Halton points to use to generate the tiling
            tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);
            title('Halton set');
            fig_num = fig_num+1;
            % find r_D for this field
            field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
            field_avg_r_D = field_stats.avg_r_D;
            r_D = [r_D, field_avg_r_D];
            gap_size = 0;
            [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse_worst,r_lc_sparse_average] = fcn_MapGen_polytopesPredictLengthCostRatio(tiled_polytopes,gap_size)
            r_lc_max_all = [r_lc_max_all, r_lc_max];
            r_lc_max_effective_all = [r_lc_max_effective_all, r_lc_max_effective];
            r_lc_avg_all = [r_lc_avg_all, r_lc_avg];
            r_lc_avg_effective_all = [r_lc_avg_effective_all, r_lc_avg_effective];
            r_lc_iterative_all = [r_lc_iterative_all, r_lc_iterative];
            r_lc_iterative_effective_all = [r_lc_iterative_effective_all, r_lc_iterative_effective];
            r_lc_sparse_worst_all = [r_lc_sparse_worst_all, r_lc_sparse_worst];
            r_lc_sparse_average_all = [r_lc_sparse_average_all, r_lc_sparse_average];
            shrink_distance = [shrink_distance, 0];
            % for radii_goals=0.25%0.02:0.02:0.1%0.001:0.010:0.1
            for size_percent = 0.85:0.0125:1
                radii_goals = field_stats.average_max_radius*size_percent;
    %             try
                des_rad = radii_goals; sigma_radius = 0; min_rad = 0.001;
                % TODO switch this to side shrinking to get gap distance as an output so it can be given to predictor as input
                [shrunk_field,mu_final,sigma_final] = fcn_MapGen_polytopesShrinkToRadius(tiled_polytopes,des_rad,sigma_radius,min_rad,fig_num);
                field_stats = fcn_MapGen_polytopesStatistics(shrunk_field);
                gap_size = field_stats.average_gap_size_G_bar;
                field_avg_r_D = field_stats.avg_r_D;
                % avg_max_rad = field_stats.average_max_radius;
                shrink_distance = [shrink_distance, gap_size];%avg_max_rad-des_rad];
                try
                   [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse_worst,r_lc_sparse_average] = fcn_MapGen_polytopesPredictLengthCostRatio(shrunk_field,gap_size)
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
                r_D = [r_D, field_avg_r_D];
           end
   %         end
   %         end
       end
   %     close all
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
        plot(r_D,r_lc_sparse_average_all,'cd')
    end
    legend('theoretical max',...
        'average from side and angle',...
        'iterative from side and angle',...
        'effective theoretical max',...
        'average from effective side and angle',...
        'iterative from effective side and angle',...
        'sparse formula, worst case',...
        'sparse formula, average case');
    figure
end