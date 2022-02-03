function [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,r_lc_sparse_worst,r_lc_sparse_average,r_lc_sparse_std,r_lc_sparse_worst_new,r_lc_sparse_average_new,r_lc_sparse_std_new,r_lc_sparse_worst_actual,r_lc_sparse_average_actual,r_lc_sparse_std_actual,r_lc_sparse_worst_linear,r_lc_sparse_average_linear,r_lc_sparse_std_linear] = fcn_MapGen_polytopesPredictLengthCostRatio(tiled_polytopes,gap_size,shrunk_distance,shrink_ang,R_bar_initial) % TODO(@sjharnett) put outputs into struct
    % fcn_MapGen_polytopesPredictLengthCostRatio
    % Given an polytope field, predict the length cost ratio from geometry
    %
    %
    %
    % FORMAT:
    %
    % [r_lc_max,r_lc_avg,r_lc_iterative,...
    % r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective,...
    % r_lc_sparse_worst,r_lc_sparse_average,r_lc_sparse_std] = ...
    % fcn_MapGen_polytopesPredictLengthCostRatio(tiled_polytopes,gap_size)
    %
    % INPUTS:
    %
    %     tiled_polytopes
    %     gap_size
    %
    %
    % OUTPUTS:
    %
    %     r_lc_max
    %     r_lc_avg
    %     r_lc_iterative
    %     r_lc_max_effective
    %     r_lc_avg_effective
    %     r_lc_iterative_effective
    %     r_lc_sparse_worst
    %     r_lc_sparse_average
    %     r_lc_sparse_std
    %
    % DEPENDENCIES:
    %
    %     fcn_MapGen_polytopeFindVertexAngles
    %     fcn_MapGen_polytopesStatistics
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

    % TODO add outputs for chosen angle, chosen side length etc

    %% Debugging and Input checks

    % set an environment variable on your machine with the getenv function...
    % in the Matlab console.  Char array of '1' will be true and '0' will be false.
    flag_check_inputs = getenv('ENV_FLAG_CHECK_INPUTS');  % '1' will check input args
    flag_do_plot = getenv('ENV_FLAG_DO_PLOT'); % '1' will make plots
    flag_do_debug = getenv('ENV_FLAG_DO_DEBUG'); % '1' will enable debugging

    % if the char array has length 0, assume the env var isn't set and default to...
    % dipslaying more information rather than potentially hiding an issue
    if length(flag_check_inputs) == 0
        flag_check_inputs = '1';
    end
    if length(flag_do_plot) == 0
        flag_do_plot = '1';
    end
    if length(flag_do_debug) == 0
        flag_do_debug = '1';
    end

    % convert flag from char string to logical
    flag_check_inputs = flag_check_inputs == '1';
    flag_do_plot = flag_do_plot == '1';
    flag_do_debug = flag_do_debug == '1';

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
        if nargin < 2 || nargin > 2
            error('Incorrect number of input arguments')
        end

        % Check the test_point input, make sure it is 'positive_1column_of_numbers' type, size 1
        fcn_MapGen_checkInputsToFunctions(column_of_numbers_test, 'positive_1column_of_numbers',1);

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
    field_stats = fcn_MapGen_polytopesStatistics(tiled_polytopes);
    field_avg_r_D = field_stats.avg_r_D;
    field_away_normals = [];
    field_away_angles = [];
    field_away_angle_vertex_normal_to_travel_direction = [];
    field_small_choice_angles = [];
    field_big_choice_angles = [];
    field_chosen_side_length = [];
    % initialize path for iterative solution
    path = NaN(1,2);
    iterative_chosen_side_lengths = [];
    iterative_small_choice_angles = [];
    start_not_found = true;
    % begin looping through polytopes in a field
    for j=1:length(tiled_polytopes)
        shrinker = tiled_polytopes(j);
        if flag_do_plot
            [angles, unit_in_vectors, unit_out_vectors] =...
            fcn_MapGen_polytopeFindVertexAngles(...
            shrinker.vertices,fig_num);
        else
            [angles, unit_in_vectors, unit_out_vectors] =...
            fcn_MapGen_polytopeFindVertexAngles(...
            shrinker.vertices);
        end
        assert(1000*eps>abs(360-sum(angles)*180/pi));
        vertices = shrinker.vertices;
        mean_vectors = (unit_out_vectors-unit_in_vectors)/2;
        length_mean_vectors = sum(mean_vectors.^2,2).^0.5;
        unit_direction_of_cut = mean_vectors./length_mean_vectors;
        if flag_do_plot
            vertex_normal_plot = quiver(vertices(1:end-1,1),vertices(1:end-1,2),unit_direction_of_cut(:,1),unit_direction_of_cut(:,2),'g');
        end
        vertex_normal_vectors = [unit_direction_of_cut(:,1),unit_direction_of_cut(:,2)];
        % vertex_normal_vectors = circshift(vertex_normal_vectors,1)
        side_vectors = zeros(length(vertices)-1,2);
        theta_normal_to_side = zeros(length(vertices)-1,1);
        angle_vertex_normal_to_travel_direction = NaN(length(vertices)-1,1);
        vertex_is_away = zeros(1, length(vertices)-1);
        travel_direction_is_within_polytope = zeros(1, length(vertices)-1);
        travel_direction = [1,0];
        away_normals = NaN(length(vertices)-1,2);
        away_vertices = NaN(length(vertices)-1,2);
        away_angles = NaN(length(vertices)-1,1);
        chosen_side_lengths = NaN(length(vertices)-1,1);
        away_angle_vertex_normal_to_travel_direction = NaN(length(vertices)-1,1);
        % make distances array circular by putting the first entry at the end to support end+1 circular indexing
        shrinker.distances = [shrinker.distances;shrinker.distances(1)];
        % end initialize path for iterative solution
        % begin looping through vertices
        for i=1:length(vertices)-1
            angle_vertex_normal_to_travel_direction(i) = INTERNAL_angle_between_vectors(vertex_normal_vectors(i,:),travel_direction);
            vertex_is_away(i) = dot(vertex_normal_vectors(i,:),travel_direction);
            travel_direction_is_within_polytope(i) = angle_vertex_normal_to_travel_direction(i)<=(pi-angles(i))/2;
            % begin looping through away vertices
            if vertex_is_away(i)>0 && travel_direction_is_within_polytope(i)
                % pick a point at x=0 y=/= 0 to start iterative solution
                if vertices(i,1) == 0 && vertices(i,2) ~= 0 && start_not_found
                    start_location = vertices(i,:);
                    path = [path; start_location];
                    path(any(isnan(path),2),:)=[];
                    start_not_found = false;
                end
                % end pick a point at x=0 y=/= 0 to start iterative solution
                away_vertices(i,:) = vertices(i,:);
                away_normals(i,:) = vertex_normal_vectors(i,:);
                away_angles(i) = pi-angles(i);
                away_angle_vertex_normal_to_travel_direction(i) =  angle_vertex_normal_to_travel_direction(i);
            % for the chosen divergence angle, find the length of the side we turned towards
                if INTERNAL_is_left_turn_smaller(vertex_normal_vectors(i,:),travel_direction)
                     side_length = shrinker.distances(i+1); % the vertex we're at has the side length to its left at the next index
                 else
                     side_length = shrinker.distances(i); % the vertex we're at has the side length associated with it on the right
                end
                chosen_side_lengths(i) = side_length;
                % log chosen_side_length and small angle separately for iterative solution
                if vertices(i,:) == path(end,:)
                    iterative_chosen_side_lengths = [iterative_chosen_side_lengths, side_length];
                    iterative_small_choice_angles = [iterative_small_choice_angles, away_angles(i)./2-away_angle_vertex_normal_to_travel_direction(i)];
                    if INTERNAL_is_left_turn_smaller(vertex_normal_vectors(i,:),travel_direction)
                        try
                            next_location = vertices(i-1,:);
                        catch
                            next_location = vertices(end,:);
                        end
                    else
                        next_location = vertices(i+1,:);
                    end
                    path = [path; next_location];
                end
                % end log chosen_side_lengths separately for iterative solution
            % end looping through away vertices
            end
            side_vectors(i,1) = vertices(i+1,1)-vertices(i,1);
            side_vectors(i,2) = vertices(i+1,2)-vertices(i,2);
            theta_normal_to_side(i) = INTERNAL_angle_between_vectors(side_vectors(i,:),vertex_normal_vectors(i,:));
            if flag_do_plot
                side_vec_plot = quiver(vertices(i,1),vertices(i,2),vertices(i+1,1)-vertices(i,1),vertices(i+1,2)-vertices(i,2),'-b');
                travel_vec_plot = quiver(vertices(i,1),vertices(i,2),travel_direction(1)*0.05,travel_direction(2)*0.05,'-m');
            end
        % end looping through vertices
        end
        % remove data for angles pointing towards travel direction
        away_normals(any(isnan(away_angles),2),:)=[];
        away_angles(any(isnan(away_angles),2),:)=[];
        away_angle_vertex_normal_to_travel_direction(any(isnan(away_angle_vertex_normal_to_travel_direction),2),:)=[];
        chosen_side_lengths(any(isnan(chosen_side_lengths),2),:)=[];
        % find large and small choice angles
        small_choice_angles = away_angles./2-away_angle_vertex_normal_to_travel_direction;
        big_choice_angles = away_angles./2+away_angle_vertex_normal_to_travel_direction;
        % log things for the entire field
        field_away_normals = [field_away_normals;away_normals];
        field_away_angles = [field_away_angles;away_angles];
        field_away_angle_vertex_normal_to_travel_direction = [field_away_angle_vertex_normal_to_travel_direction;away_angle_vertex_normal_to_travel_direction];
        field_small_choice_angles = [field_small_choice_angles;small_choice_angles];
        field_big_choice_angles = [field_big_choice_angles;big_choice_angles];
        field_chosen_side_length = [field_chosen_side_length;chosen_side_lengths];
        % change chosen angle to effective angle based on gap size
        field_away_angles_effective = field_away_angles - atan2(gap_size,field_chosen_side_length);
        field_small_choice_angles_effective = field_small_choice_angles - atan2(gap_size,field_chosen_side_length);
        % change chosen side length to effective angle side length based on gap size
        field_chosen_side_length_effective = (field_chosen_side_length.^2+gap_size.^2).^0.5;

        % TODO put all plotting code behind debug flag
        % plot
        % for i=1:length(away_angles)
        %     away_travel_vec_plot = quiver(away_vertices(i,1),away_vertices(i,2),travel_direction(1)*0.05,travel_direction(2)*0.05,'-m');
        % end
        % for i=1:length(vertices)-1
        % side_vectors(i,1) = vertices(i+1,1)-vertices(i,1);
        % side_vectors(i,2) = vertices(i+1,2)-vertices(i,2);

        % theta = arccos((x dot y)/(mag x * mag y))
        % theta_normal_to_side(i) = INTERNAL_angle_between_vectors(side_vectors(i,:),vertex_normal_vectors(i,:));
        % label large and small sizes
%         size = max(max(vertices)) - min(min(vertices));
%         nudge = size*0.01;
%         % Label the vertices
%         for ith_angle = 1:length(angles)
%             ith_vertex = ith_angle;
%             text(vertices(ith_vertex,1)+nudge,vertices(ith_vertex,2),...
%                 sprintf('%.0f deg',180-angles(ith_angle,1)*180/pi));
%         end
    end
    % begin results generation
    r_lc_max = 1/cos(mean(field_away_angles)/2);
    r_lc_max_effective = 1/cos(mean(field_away_angles_effective)/2);
    r_lc_avg = 1/cos(mean(field_small_choice_angles));
    r_lc_avg_effective = 1/cos(mean(field_small_choice_angles_effective));
    field_traveled_distance_L = cos(field_small_choice_angles).*field_chosen_side_length;
    field_traveled_distance_L_effective = cos(field_small_choice_angles_effective).*field_chosen_side_length_effective;
    field_path_distance_H = field_chosen_side_length;
    field_path_distance_H_effective = field_chosen_side_length_effective;
    % TODO debug here
    L_E = 0;
    L_P = 0;
    L_E_effective = 0;
    L_P_effective = 0;
    i = 1;
    % this loops through polytopes in left to right order
    % TODO use vertices to loop through polytopes in path order
        % iterative solution
        % pick a point at x=0 y=/= 0
        % log the short side and short angle
        % note the vertex you wind up at
        % go to next polytope, if the general conditions aren't met AND that vertex is not present, skip
        % repeat until x=1
    try
        while L_E < 1
            L_E = L_E + field_traveled_distance_L(i);
            L_E_effective = L_E_effective + field_traveled_distance_L_effective(i);
            L_P = L_P + field_path_distance_H(i);
            L_P_effective = L_P_effective + field_path_distance_H_effective(i);
            i = i + 1;
        end
    catch
        L_E_effective = 1;
        L_P_effective = 1;
    end
    r_lc_iterative = L_P/L_E;
    r_lc_iterative_effective = L_P_effective/L_E_effective;
    side_and_ang = [field_chosen_side_length, field_small_choice_angles];
    rounded_side_and_ang = round(side_and_ang,6);
    % remove zeros (usually from zero angle divergences)
    rounded_side_and_ang( any(rounded_side_and_ang==0,2), : ) = [];
    divergence_heights = rounded_side_and_ang(:,1).*sin(rounded_side_and_ang(:,2));
    linear_density = field_stats.linear_density;
    linear_density_int = round(linear_density,0);
    N_int_from_shrink_dist = (sqrt(field_stats.unoccupancy_ratio)*1)/(2*shrunk_distance*sind(shrink_ang/2));
    N_int_actual = field_stats.linear_density_mean;
    N_int_linear = 1.3*linear_density*(1-shrunk_distance/R_bar_initial);
    s_divergence_heights = sort(divergence_heights);
    worst_divergence_heights = s_divergence_heights(end-linear_density_int+1:end);
    worst_divergence_heights_new = s_divergence_heights(end-N_int_from_shrink_dist+1:end);
    worst_divergence_heights_actual = s_divergence_heights(end-N_int_actual+1:end);
    worst_divergence_heights_linear = s_divergence_heights(end-N_int_linear+1:end);
    % TODO @sjharnett change /1 to be over L_P
    r_lc_sparse_worst = sum((worst_divergence_heights.^2+(1/linear_density_int)^2).^0.5)/1;
    r_lc_sparse_average = linear_density_int*mean((s_divergence_heights.^2+(1/linear_density_int)^2).^0.5)/1;
    r_lc_sparse_std = linear_density_int*std((s_divergence_heights.^2+(1/linear_density_int)^2).^0.5)/1;
    r_lc_sparse_worst_new = sum((worst_divergence_heights_new.^2+(1/N_int_from_shrink_dist)^2).^0.5)/1;
    r_lc_sparse_average_new = N_int_from_shrink_dist*mean((s_divergence_heights.^2+(1/N_int_from_shrink_dist)^2).^0.5)/1;
    r_lc_sparse_std_new = N_int_from_shrink_dist*std((s_divergence_heights.^2+(1/N_int_from_shrink_dist)^2).^0.5)/1;
    r_lc_sparse_worst_actual = sum((worst_divergence_heights_actual.^2+(1/N_int_actual)^2).^0.5)/1;
    r_lc_sparse_average_actual = N_int_actual*mean((s_divergence_heights.^2+(1/N_int_actual)^2).^0.5)/1;
    r_lc_sparse_std_actual = N_int_actual*std((s_divergence_heights.^2+(1/N_int_actual)^2).^0.5)/1;
    r_lc_sparse_worst_linear = sum((worst_divergence_heights_linear.^2+(1/N_int_linear)^2).^0.5)/1;
    r_lc_sparse_average_linear = N_int_linear*mean((s_divergence_heights.^2+(1/N_int_linear)^2).^0.5)/1;
    r_lc_sparse_std_linear = N_int_linear*std((s_divergence_heights.^2+(1/N_int_linear)^2).^0.5)/1;

    if flag_do_plot
        % end results generation
        figure;
        hold on;
        histogram(field_big_choice_angles*180/pi,'BinWidth',2,'FaceColor','g','FaceAlpha',0.4)
        histogram(field_small_choice_angles*180/pi,'BinWidth',2,'FaceColor','b','FaceAlpha',0.4)
        legend('large, unchosen divergence angles','small, chosen divergence angles')
        xlabel('interior angle [deg]')
        ylabel('count')
        title('Histogram of Divergence Angles')
        figure;
        hold on;
        histogram(field_big_choice_angles*180/pi,'BinWidth',2,'FaceColor','g','FaceAlpha',0.4)
        histogram(field_small_choice_angles*180/pi,'BinWidth',2,'FaceColor','b','FaceAlpha',0.4)
        histogram(field_away_angles/2*180/pi,'BinWidth',2,'FaceColor','r','FaceAlpha',0.4)
        legend('large, unchosen divergence angles','small, chosen divergence angles','all away facing polytope angles, halved')
        xlabel('interior angle [deg]')
        ylabel('count')
        title('Histogram of Divergence Angles')
        figure;
        hold on;
        histogram(field_big_choice_angles*180/pi,'BinWidth',2,'FaceColor','g','FaceAlpha',0.4)
        histogram(field_small_choice_angles*180/pi,'BinWidth',2,'FaceColor','b','FaceAlpha',0.4)
        histogram(field_away_angles/2*180/pi,'BinWidth',2,'FaceColor','r','FaceAlpha',0.4)
        histogram(field_small_choice_angles_effective/2*180/pi,'BinWidth',2,'FaceColor','c','FaceAlpha',0.4)
        legendstr = sprintf('effective chosen angle for gap size %.1f',gap_size);
        legend('large, unchosen divergence angles','small, chosen divergence angles','all away facing polytope angles, halved',legendstr)
        xlabel('interior angle [deg]')
        ylabel('count')
        title('Histogram of Divergence Angles')
        figure;
        hold on;
        histogram(field_chosen_side_length,'BinWidth',2,'FaceColor','b','FaceAlpha',0.4)
        histogram(field_chosen_side_length_effective,'BinWidth',2,'FaceColor','c','FaceAlpha',0.4)
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
    % ang = atan2(norm(cross([a,0],[b,0])), dot(a,b));
    ang = acos(dot(a,b)/(norm(a)*norm(b)));
    if a(2)<0
        ang = ang;
    end
end

function left_turn_is_smaller = INTERNAL_is_left_turn_smaller(vertex_normal,travel_direction)
    cross_prod = cross([vertex_normal,0],[travel_direction,0]);
    left_turn_is_smaller = cross_prod(3)>0;
end
