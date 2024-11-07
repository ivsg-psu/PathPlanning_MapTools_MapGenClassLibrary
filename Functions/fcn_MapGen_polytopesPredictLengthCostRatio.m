function [field_small_choice_angles,field_big_choice_angles,r_lc_estimates] =...
    fcn_MapGen_polytopesPredictLengthCostRatio(pre_shrink_polytopes,polytopes,gap_size,travel_direction,L_E)
    % fcn_MapGen_polytopesPredictLengthCostRatio
    % Given an polytope field, predict the length cost ratio from geometry
    %
    %
    %
    % FORMAT:
    % [field_small_choice_angles,field_big_choice_angles,r_lc_estimates] =...
    %  fcn_MapGen_polytopesPredictLengthCostRatio(pre_shrink_polytopes,polytopes,gap_size,travel_direction)
    %
    % INPUTS:
    %     pre_shrink_polytopes - the fully tiled field
    %     polytopes - the field after shrinking has been applied
    %     gap_size - the commanded gap size used when shrinking was applied
    %     travel_direction - the unit vectory pointing in the direciton of travel
    %         commonly [1,0] for a path from left to right
    %     L_E - Euclidean distance from start to goal.  The quantity in the denomonator of length cost ratio
    %
    %
    % OUTPUTS:
    %
    %
    %     field_small_choice_angles - array of all smaller (likely chosen) deflection angles for the field
    %     field_big_choice_angles - array of all larger (likely unchosen) deflection angles for the field
    %     r_lc_estimates - struct of different cost ratio estimation methods
    %     r_lc_estimates - struct of different cost ratio estimation methods
    %
    % DEPENDENCIES:
    %
    %     fcn_MapGen_polytopeFindVertexAngles
    %     fcn_MapGen_polytopesStatistics
    %     fcn_MapGen_polytopesPredictUnoccupancyRatio
    %
    % EXAMPLES:
    %
    % See the script: script_fcn_MapGen_polytopesPredictLengthCostRatio.m
    % for a full test suite.
    %
    % Questions or comments? contact sjh6473@psu.edu

    % REVISION HISTORY:
    % 2021_10_22
    % -- first written by Steve Harnett


    %% Debugging and Input checks
    flag_check_inputs = 1; % Set equal to 1 to check the input arguments
    flag_do_plot = 0;      % Set equal to 1 for plotting
    flag_do_debug = 0;     % Set equal to 1 for debugging
    if flag_do_debug
        fig_for_debug = 747;
        fig_num = 2100;
        st = dbstack; %#ok<*UNRCH>
        fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    end

    %% check input arguments?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   _____                   _
    %  |_   _|                 | |
    %    | |  _ __  _ __  _   _| |_ ___
    %    | | | '_ \| '_ \| | | | __/ __|
    %   _| |_| | | | |_) | |_| | |_\__ \
    %  |_____|_| |_| .__/ \__,_|\__|___/
    %              | |
    %              |_|
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if 1 == flag_check_inputs

        % Are there the right number of inputs?
        if nargin < 5 || nargin > 5
            error('Incorrect number of input arguments')
        end

    end

    %% Start of main code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   __  __       _
    %  |  \/  |     (_)
    %  | \  / | __ _ _ _ __
    %  | |\/| |/ _` | | '_ \
    %  | |  | | (_| | | | | |
    %  |_|  |_|\__,_|_|_| |_|
    %
    %See: http://patorjk.com/software/taag/#p=display&f=Big&t=Main
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

    fig_num = 12;
    field_stats_pre_shrink = fcn_MapGen_polytopesStatistics(pre_shrink_polytopes);
    R_bar_initial = field_stats_pre_shrink.average_max_radius;
    shrink_ang = field_stats_pre_shrink.average_vertex_angle;
    field_stats = fcn_MapGen_polytopesStatistics(polytopes);
    field_avg_r_D = field_stats.avg_r_D;
    field_away_normals = [];
    field_away_angles = [];
    field_away_angle_vertex_normal_to_travel_direction = [];
    field_small_choice_angles = [];
    field_big_choice_angles = [];
    field_chosen_side_length = [];
    for j=1:length(polytopes)
        shrinker = polytopes(j);
        vertices = shrinker.vertices;
        if flag_do_plot
            [angles, unit_in_vectors, unit_out_vectors] =...
            fcn_MapGen_polytopeFindVertexAngles(...
            vertices,fig_num);
        else
            [angles, unit_in_vectors, unit_out_vectors] =...
            fcn_MapGen_polytopeFindVertexAngles(...
            vertices);
        end
        mean_vectors = (unit_out_vectors-unit_in_vectors)/2;
        length_mean_vectors = sum(mean_vectors.^2,2).^0.5;
        unit_direction_of_cut = mean_vectors./length_mean_vectors;
        if flag_do_plot
            % plots interior vertex normal
            vertex_normal_plot = quiver(vertices(1:end-1,1),vertices(1:end-1,2),unit_direction_of_cut(:,1),unit_direction_of_cut(:,2),'g');
        end
        % interior vertex normals
        vertex_normal_vectors = [unit_direction_of_cut(:,1),unit_direction_of_cut(:,2)];
        % vectors in directions of polytope sides
        side_vectors = zeros(length(vertices)-1,2);
        % this will store the direction from the vertex normal to the side
        theta_normal_to_side = zeros(length(vertices)-1,1);
        % this will store the angle between the vertex normal and travel direction
        angle_vertex_normal_to_travel_direction = NaN(length(vertices)-1,1);
        % a boolean determining if vertex opens away from travel direction
        vertex_is_away = zeros(1, length(vertices)-1);
        % a boolean determining if this vertex "blocks" the direction of travel
        travel_direction_is_within_polytope = zeros(1, length(vertices)-1);
        % these will store normals, vertices, and angles that block travel and are valid decision points
        away_normals = NaN(length(vertices)-1,2);
        away_vertices = NaN(length(vertices)-1,2);
        away_angles = NaN(length(vertices)-1,1);
        % this will store side lengths associated with the chosen deflection direction (the small angle)
        chosen_side_lengths = NaN(length(vertices)-1,1);
        away_angle_vertex_normal_to_travel_direction = NaN(length(vertices)-1,1);
        % make distances array circular by putting the first entry at the end to support end+1 circular indexing
        shrinker.distances = [shrinker.distances;shrinker.distances(1)];

        % begin looping through vertices
        for i=1:length(vertices)-1
            angle_vertex_normal_to_travel_direction(i) = INTERNAL_angle_between_vectors(vertex_normal_vectors(i,:),travel_direction);
            vertex_is_away(i) = dot(vertex_normal_vectors(i,:),travel_direction);
            travel_direction_is_within_polytope(i) = angle_vertex_normal_to_travel_direction(i)<=(pi-angles(i))/2;
            % begin looping through away vertices that block travel direction
            if vertex_is_away(i)>0 && travel_direction_is_within_polytope(i)
                away_vertices(i,:) = vertices(i,:);
                away_normals(i,:) = vertex_normal_vectors(i,:);
                % angles need to be subtracted from 180 deg to get interior angles
                away_angles(i) = pi-angles(i);
                away_angle_vertex_normal_to_travel_direction(i) =  angle_vertex_normal_to_travel_direction(i);
                % for the chosen divergence angle, find the length of the side we turned towards
                if INTERNAL_is_left_turn_smaller(vertex_normal_vectors(i,:),travel_direction)
                    % the vertex we're at has the side length to its left at the next index
                    side_length = shrinker.distances(i+1);
                 else
                    % the vertex we're at has the side length associated with it on the right
                    side_length = shrinker.distances(i);
                end
                chosen_side_lengths(i) = side_length;

            % end looping through away vertices
            end
            side_vectors(i,1) = vertices(i+1,1)-vertices(i,1);
            side_vectors(i,2) = vertices(i+1,2)-vertices(i,2);
            theta_normal_to_side(i) = INTERNAL_angle_between_vectors(side_vectors(i,:),vertex_normal_vectors(i,:));
            if flag_do_plot
                % plot the side vectors on top of the original polytope to check for discrepancies
                side_vec_plot = quiver(vertices(i,1),vertices(i,2),vertices(i+1,1)-vertices(i,1),vertices(i+1,2)-vertices(i,2),'-b');
                travel_vec_plot = quiver(vertices(i,1),vertices(i,2),travel_direction(1)*0.05,travel_direction(2)*0.05,'-m');
            end
        % end looping through vertices
        end
        % remove data for angles pointing towards travel direction, these aren't validly blocking travel
        away_normals(any(isnan(away_angles),2),:)=[];
        away_angles(any(isnan(away_angles),2),:)=[];
        away_angle_vertex_normal_to_travel_direction(any(isnan(away_angle_vertex_normal_to_travel_direction),2),:)=[];
        chosen_side_lengths(any(isnan(chosen_side_lengths),2),:)=[];
        % find large and small deflection choice angles
        small_choice_angles = away_angles./2-away_angle_vertex_normal_to_travel_direction;
        big_choice_angles = away_angles./2+away_angle_vertex_normal_to_travel_direction;
        % log data for the entire field
        field_away_normals = [field_away_normals;away_normals];
        field_away_angles = [field_away_angles;away_angles];
        field_away_angle_vertex_normal_to_travel_direction = [field_away_angle_vertex_normal_to_travel_direction;away_angle_vertex_normal_to_travel_direction];
        field_small_choice_angles = [field_small_choice_angles;small_choice_angles];
        field_big_choice_angles = [field_big_choice_angles;big_choice_angles];
        field_chosen_side_length = [field_chosen_side_length;chosen_side_lengths];

        if flag_do_plot
            % code to plot travel directions on each relevant vertex with angle labels
            % with labels
            for i=1:length(away_angles)
                away_travel_vec_plot = ...
                quiver(away_vertices(i,1),away_vertices(i,2),travel_direction(1)*0.05,travel_direction(2)*0.05,'-m');
            end
            for i=1:length(vertices)-1
            side_vectors(i,1) = vertices(i+1,1)-vertices(i,1);
            side_vectors(i,2) = vertices(i+1,2)-vertices(i,2);

            theta_normal_to_side(i) = INTERNAL_angle_between_vectors(side_vectors(i,:),vertex_normal_vectors(i,:));
            % label large and small sizes
               size = max(max(vertices)) - min(min(vertices));
               nudge = size*0.01;
               % Label the vertices
               for ith_angle = 1:length(angles)
                   ith_vertex = ith_angle;
                   text(vertices(ith_vertex,1)+nudge,vertices(ith_vertex,2),...
                       sprintf('%.0f deg',180-angles(ith_angle,1)*180/pi));
               end
            end
        end
    end
    % begin results generation

    side_and_ang = [field_chosen_side_length, field_small_choice_angles];
    rounded_side_and_ang = round(side_and_ang,6);
    % remove zeros (usually from zero angle divergences)
    rounded_side_and_ang( any(rounded_side_and_ang==0,2), : ) = [];
    divergence_heights = rounded_side_and_ang(:,1).*sin(rounded_side_and_ang(:,2));
    linear_density = field_stats.linear_density;
    linear_density_int = round(linear_density,0);
    unocc_ests = fcn_MapGen_polytopesPredictUnoccupancyRatio(pre_shrink_polytopes,polytopes,gap_size);
    N_int_from_shrink_dist = (sqrt(unocc_ests.A_unocc_est_poly_fit)*1)/(gap_size*sind(shrink_ang/2));
    N_int_actual = field_stats.linear_density_mean;
    N_int_linear = linear_density*(1-(gap_size/2)/R_bar_initial);
    s_divergence_heights = sort(divergence_heights);
    worst_divergence_heights = s_divergence_heights(end-linear_density_int+1:end);
    worst_divergence_heights_new = s_divergence_heights(end-N_int_from_shrink_dist+1:end);
    worst_divergence_heights_actual = s_divergence_heights(end-N_int_actual+1:end);
    worst_divergence_heights_linear = s_divergence_heights(end-N_int_linear+1:end);

    r_lc_estimates.r_lc_sparse_worst = sum((worst_divergence_heights.^2+(1/linear_density_int)^2).^0.5)/L_E;
    r_lc_estimates.r_lc_sparse_average = linear_density_int*mean((s_divergence_heights.^2+(1/linear_density_int)^2).^0.5)/L_E;
    r_lc_estimates.r_lc_sparse_std = linear_density_int*std((s_divergence_heights.^2+(1/linear_density_int)^2).^0.5)/L_E;
    r_lc_estimates.r_lc_sparse_worst_new = sum((worst_divergence_heights_new.^2+(1/N_int_from_shrink_dist)^2).^0.5)/L_E;
    r_lc_estimates.r_lc_sparse_average_new = N_int_from_shrink_dist*mean((s_divergence_heights.^2+(1/N_int_from_shrink_dist)^2).^0.5)/L_E;
    r_lc_estimates.r_lc_sparse_std_new = N_int_from_shrink_dist*std((s_divergence_heights.^2+(1/N_int_from_shrink_dist)^2).^0.5)/L_E;
    r_lc_estimates.r_lc_sparse_worst_actual = sum((worst_divergence_heights_actual.^2+(1/N_int_actual)^2).^0.5)/L_E;
    r_lc_estimates.r_lc_sparse_average_actual = N_int_actual*mean((s_divergence_heights.^2+(1/N_int_actual)^2).^0.5)/L_E;
    r_lc_estimates.r_lc_sparse_std_actual = N_int_actual*std((s_divergence_heights.^2+(1/N_int_actual)^2).^0.5)/L_E;
    r_lc_estimates.r_lc_sparse_worst_linear = sum((worst_divergence_heights_linear.^2+(1/N_int_linear)^2).^0.5)/L_E;
    r_lc_estimates.r_lc_sparse_average_linear = N_int_linear*mean((s_divergence_heights.^2+(1/N_int_linear)^2).^0.5)/L_E;
    r_lc_estimates.r_lc_sparse_std_linear = N_int_linear*std((s_divergence_heights.^2+(1/N_int_linear)^2).^0.5)/L_E;

    if flag_do_plot
        figure;
        hold on;
        histogram(field_big_choice_angles*180/pi,'BinWidth',2,'FaceColor','r','FaceAlpha',0.4)
        histogram(field_small_choice_angles*180/pi,'BinWidth',2,'FaceColor','g','FaceAlpha',0.4)
        legend('large, unchosen divergence angles','small, chosen divergence angles')
        xlabel('interior angle [deg]')
        ylabel('count')
        title('Histogram of Divergence Angles')
        figure;
        hold on;
        histogram(field_big_choice_angles*180/pi,'BinWidth',2,'FaceColor','r','FaceAlpha',0.4)
        histogram(field_small_choice_angles*180/pi,'BinWidth',2,'FaceColor','g','FaceAlpha',0.4)
        histogram(field_away_angles/2*180/pi,'BinWidth',2,'FaceColor','b','FaceAlpha',0.4)
        legend('large, unchosen divergence angles','small, chosen divergence angles','all away facing polytope angles, halved')
        xlabel('interior angle [deg]')
        ylabel('count')
        title('Histogram of Divergence Angles')
        figure;
        hold on;
        histogram(field_chosen_side_length,'BinWidth',2,'FaceColor','g','FaceAlpha',0.4)
        legendstr_side = sprintf('effective chosen side length for gap size %.1f',gap_size);
        legend('chosen side length', legendstr_side)
        xlabel('side length [m]')
        ylabel('count')
        title('Histogram of Traveled Side Lengths')
    end
end

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

function ang = INTERNAL_angle_between_vectors(a,b)
    % note this function returns [0,180] so is not directional (i.e. a
    % vector that is rotated 182 deg CCW will return a value of 178 CW
    % this doesn't matter as we expect all vertex normal vectors we care
    % about to be in the direction of the path
    % theta = arccos((x dot y)/(mag x * mag y))
    ang = acos(dot(a,b)/(norm(a)*norm(b)));
    if a(2)<0
        ang = ang;
    end
end

function left_turn_is_smaller = INTERNAL_is_left_turn_smaller(vertex_normal,travel_direction)
    % leverage the right hand rule (RHR) to determine the position of the
    % smaller deflection angle relative to
    % an actor looking in the direciton of travel
    cross_prod = cross([vertex_normal,0],[travel_direction,0]);
    left_turn_is_smaller = cross_prod(3)>0;
end
