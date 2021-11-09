% fcn_MapGen_polytopesPredictLengthCostRatio

% REVISION HISTORY:
% 2021_10_22
% -- first written by S. Harnett
% TODO add outputs for chosen angle, chosen side length etc
function [r_lc_max,r_lc_avg,r_lc_iterative,r_lc_max_effective,r_lc_avg_effective,r_lc_iterative_effective] = fcn_MapGen_polytopesPredictLengthCostRatio(tiled_polytopes,gap_size)
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
        [angles, unit_in_vectors, unit_out_vectors] =...
            fcn_MapGen_polytopeFindVertexAngles(...
            shrinker.vertices,fig_num);
        assert(1000*eps>abs(360-sum(angles)*180/pi));
        vertices = shrinker.vertices;
        mean_vectors = (unit_out_vectors-unit_in_vectors)/2;
        length_mean_vectors = sum(mean_vectors.^2,2).^0.5;
        unit_direction_of_cut = mean_vectors./length_mean_vectors;
        vertex_normal_plot = quiver(vertices(1:end-1,1),vertices(1:end-1,2),unit_direction_of_cut(:,1),unit_direction_of_cut(:,2),'g')
        vertex_normal_vectors = [unit_direction_of_cut(:,1),unit_direction_of_cut(:,2)]
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
            angle_vertex_normal_to_travel_direction(i) = angle_between_vectors(vertex_normal_vectors(i,:),travel_direction);
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
                if is_left_turn_smaller(vertex_normal_vectors(i,:),travel_direction)
                     side_length = shrinker.distances(i+1); % the vertex we're at has the side length to its left at the next index
                 else
                     side_length = shrinker.distances(i); % the vertex we're at has the side length associated with it on the right
                end
                chosen_side_lengths(i) = side_length;
                % log chosen_side_length and small angle separately for iterative solution
                if vertices(i,:) == path(end,:)
                    iterative_chosen_side_lengths = [iterative_chosen_side_lengths, side_length];
                    iterative_small_choice_angles = [iterative_small_choice_angles, away_angles(i)./2-away_angle_vertex_normal_to_travel_direction(i)];
                    if is_left_turn_smaller(vertex_normal_vectors(i,:),travel_direction)
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
            theta_normal_to_side(i) = angle_between_vectors(side_vectors(i,:),vertex_normal_vectors(i,:));
            side_vec_plot = quiver(vertices(i,1),vertices(i,2),vertices(i+1,1)-vertices(i,1),vertices(i+1,2)-vertices(i,2),'-b');
            travel_vec_plot = quiver(vertices(i,1),vertices(i,2),travel_direction(1)*0.05,travel_direction(2)*0.05,'-m');
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
        % theta_normal_to_side(i) = angle_between_vectors(side_vectors(i,:),vertex_normal_vectors(i,:));
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
    while L_E < 1
        L_E = L_E + field_traveled_distance_L(i);
        L_E_effective = L_E_effective + field_traveled_distance_L_effective(i);
        L_P = L_P + field_path_distance_H(i);
        L_P_effective = L_P_effective + field_path_distance_H_effective(i);
        i = i + 1;
    end
    r_lc_iterative = L_P/L_E;
    r_lc_iterative_effective = L_P_effective/L_E_effective;
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
% TODO make internal function
function ang = angle_between_vectors(a,b)
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
% TODO make internal function
function left_turn_is_smaller = is_left_turn_smaller(vertex_normal,travel_direction)
    cross_prod = cross([vertex_normal,0],[travel_direction,0]);
    left_turn_is_smaller = cross_prod(3)>0;
end
