% script_test_fcn_MapGen_polytopeFindVertexAngles
% Tests function: fcn_MapGen_polytopeFindVertexAngles

% REVISION HISTORY:
% 2021_08_01
% -- first written by S. Brennan

function main()
    %% Basic example of vertex calculation - a square
    fig_num = 1;
    vertices = [0 0; 1 0; 1 1; 0 1; 0 0];
    [angles, unit_in_vectors, unit_out_vectors] =...
        fcn_MapGen_polytopeFindVertexAngles(...
        vertices,fig_num);
    assert(1000*eps>abs(360-sum(angles)*180/pi));
    mean_vectors = (unit_out_vectors-unit_in_vectors)/2;
    length_mean_vectors = sum(mean_vectors.^2,2).^0.5;
    unit_direction_of_cut = mean_vectors./length_mean_vectors;
    quiver(vertices(1:end-1,1),vertices(1:end-1,2),unit_direction_of_cut(:,1),unit_direction_of_cut(:,2),'g')
    vertex_normal_vectors = [unit_direction_of_cut(:,1)-vertices(1:end-1,1),unit_direction_of_cut(:,2)-vertices(1:end-1,2)];
    side_vectors = zeros(length(vertices)-1,2);
    for i=1:length(vertices)-1
        side_vectors(i,1) = vertices(i+1,1)-vertices(i,1);
        side_vectors(i,2) = vertices(i+1,2)-vertices(i,2);
        quiver(vertices(i,1),vertices(i,2),vertices(i+1,1)-vertices(i,1),vertices(i+1,2)-vertices(i,2),'-b');
    end
    %% Basic example of vertex calculation - a triangle
    fig_num = 2;
    vertices = [0 0; 1 1; 0 1; 0 0];
    [angles, unit_in_vectors, unit_out_vectors] =...
        fcn_MapGen_polytopeFindVertexAngles(...
        vertices,fig_num);
    assert(1000*eps>abs(360-sum(angles)*180/pi));
    mean_vectors = (unit_out_vectors-unit_in_vectors)/2;
    length_mean_vectors = sum(mean_vectors.^2,2).^0.5;
    unit_direction_of_cut = mean_vectors./length_mean_vectors;
    quiver(vertices(1:end-1,1),vertices(1:end-1,2),unit_direction_of_cut(:,1),unit_direction_of_cut(:,2),'g')
    vertex_normal_vectors = [unit_direction_of_cut(:,1)-vertices(1:end-1,1),unit_direction_of_cut(:,2)-vertices(1:end-1,2)];
    side_vectors = zeros(length(vertices)-1,2);
    for i=1:length(vertices)-1
        side_vectors(i,1) = vertices(i+1,1)-vertices(i,1);
        side_vectors(i,2) = vertices(i+1,2)-vertices(i,2);
        quiver(vertices(i,1),vertices(i,2),vertices(i+1,1)-vertices(i,1),vertices(i+1,2)-vertices(i,2),'-b');
    end
    %% Random polytope calculation
    % Set up polytopes
    polytopes = fcn_MapGen_haltonVoronoiTiling([1 100],[1 1]);


    bounding_box = [0,0; 1,1];
    trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box);

    % Pick a random polytope
    Npolys = length(trim_polytopes);
    rand_poly = 1+floor(rand*Npolys);
    shrinker = trim_polytopes(rand_poly);

    % Basic example of vertex calculation
    fig_num = 11;
    fig_num = 12;
    Halton_range = [1 500]; % range of Halton points to use to generate the tiling
    tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);
    title('Halton set');
    field_away_normals = [];
    field_away_angles = [];
    field_away_angle_vertex_normal_to_travel_direction = [];
    field_small_choice_angles = [];
    field_big_choice_angles = [];
    field_chosen_side_length = [];
    for i=1:length(tiled_polytopes)
        shrinker = tiled_polytopes(i);
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
        angle_vertex_normal_to_travel_direction = zeros(length(vertices)-1,1);
        is_away = zeros(1, length(vertices)-1);
        travel_direction = [1,0];
        away_normals = NaN(length(vertices)-1,2);
        away_vertices = NaN(length(vertices)-1,2);
        away_angles = NaN(length(vertices)-1,1);
        chosen_side_lengths = NaN(length(vertices)-1,1);
        away_angle_vertex_normal_to_travel_direction = NaN(length(vertices)-1,1);
        % make distances array circular by putting the first entry at the end to support end+1 circular indexing
        shrinker.distances = [shrinker.distances;shrinker.distances(1)];
        for i=1:length(vertices)-1
            is_away(i) = dot(vertex_normal_vectors(i,:),travel_direction);
            angle_vertex_normal_to_travel_direction(i) = angle_between_vectors(vertex_normal_vectors(i,:),travel_direction);
            if is_away(i)>0 && angle_vertex_normal_to_travel_direction(i)<(vertices(i)/2)
                away_vertices(i,:) = vertices(i,:);
                away_normals(i,:) = vertex_normal_vectors(i,:);
                away_angles(i) = pi-angles(i);
                away_angle_vertex_normal_to_travel_direction(i) =  angle_vertex_normal_to_travel_direction(i);
            end
            side_vectors(i,1) = vertices(i+1,1)-vertices(i,1);
            side_vectors(i,2) = vertices(i+1,2)-vertices(i,2);
            theta_normal_to_side(i) = angle_between_vectors(side_vectors(i,:),vertex_normal_vectors(i,:));
            side_vec_plot = quiver(vertices(i,1),vertices(i,2),vertices(i+1,1)-vertices(i,1),vertices(i+1,2)-vertices(i,2),'-b');
            travel_vec_plot = quiver(vertices(i,1),vertices(i,2),travel_direction(1)*0.05,travel_direction(2)*0.05,'-m');
            % for the chosen divergence angle, find the length of the side we turned towards
            if is_left_turn_smaller(vertex_normal_vectors(i,:),travel_direction)
                 side_length = shrinker.distances(i); % the vertex we're at has the side length associated with it on the left
             else
                 side_length = shrinker.distances(i+1); % the vertex we're at has the side length to its right at the next index
            end
            chosen_side_length(i) = side_length;
        end
        % remove data for angles pointing towards travel direction
        away_normals(any(isnan(away_angles),2),:)=[];
        away_angles(any(isnan(away_angles),2),:)=[];
        away_angle_vertex_normal_to_travel_direction(any(isnan(away_angle_vertex_normal_to_travel_direction),2),:)=[];
        chosen_side_length(any(isnan(away_angles),2),:)=[];
        % find large and small choice angles
        small_choice_angles = away_angles./2-away_angle_vertex_normal_to_travel_direction;
        big_choice_angles = away_angles./2+away_angle_vertex_normal_to_travel_direction;
        % log things for the entire field
        field_away_normals = [field_away_normals;away_normals];
        field_away_angles = [field_away_angles;away_angles];
        field_away_angle_vertex_normal_to_travel_direction = [field_away_angle_vertex_normal_to_travel_direction;away_angle_vertex_normal_to_travel_direction];
        field_small_choice_angles = [field_small_choice_angles;small_choice_angles];
        field_big_choice_angles = [field_big_choice_angles;big_choice_angles];
        field_chosen_side_length = [field_chosen_side_length;chosen_side_length];
        %plot
        for i=1:length(away_angles)
            away_travel_vec_plot = quiver(away_vertices(i,1),away_vertices(i,2),travel_direction(1)*0.05,travel_direction(2)*0.05,'-m');
        end
%         for i=1:length(vertices)-1
%             side_vectors(i,1) = vertices(i+1,1)-vertices(i,1);
%             side_vectors(i,2) = vertices(i+1,2)-vertices(i,2);

        % theta = arccos((x dot y)/(mag x * mag y))
%             theta_normal_to_side(i) = angle_between_vectors(side_vectors(i,:),vertex_normal_vectors(i,:));
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
end
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
function left_turn_is_smaller = is_left_turn_smaller(vertex_normal,travel_direction)
    cross_prod = cross([vertex_normal,0],[travel_direction,0]);
    left_turn_is_smaller = cross_prod(3)>0;
end