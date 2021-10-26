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
    Halton_range = [1 1000]; % range of Halton points to use to generate the tiling
    tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);
    title('Halton set');
    field_away_normals = [];
    field_away_angles = [];
    field_away_angle_vertex_normal_to_travel_direction = [];
    field_small_choice_angles = [];
    field_big_choice_angles = [];
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
        quiver(vertices(1:end-1,1),vertices(1:end-1,2),unit_direction_of_cut(:,1),unit_direction_of_cut(:,2),'g')
        vertex_normal_vectors = [unit_direction_of_cut(:,1),unit_direction_of_cut(:,2)]
        % vertex_normal_vectors = circshift(vertex_normal_vectors,1)
        side_vectors = zeros(length(vertices)-1,2);
        theta_normal_to_side = zeros(length(vertices)-1,1);
        angle_vertex_normal_to_travel_direction = zeros(length(vertices)-1,1);
        is_away = zeros(1, length(vertices)-1);
        travel_direction = [1,0];
        away_normals = NaN(length(vertices)-1,2);
        away_angles = NaN(length(vertices)-1,1);
        away_angle_vertex_normal_to_travel_direction = NaN(length(vertices)-1,1);
        for i=1:length(vertices)-1
            is_away(i) = dot(vertex_normal_vectors(i,:),travel_direction);
            angle_vertex_normal_to_travel_direction(i) = angle_between_vectors(vertex_normal_vectors(i,:),travel_direction);
            if is_away(i)>0
                away_normals(i,:) = vertex_normal_vectors(i,:);
                away_angles(i) = angles(i);
                away_angle_vertex_normal_to_travel_direction(i) =  angle_vertex_normal_to_travel_direction(i);
            end
        end
        away_normals(any(isnan(away_angles),2),:)=[];
        away_angles(any(isnan(away_angles),2),:)=[];
        away_angle_vertex_normal_to_travel_direction(any(isnan(away_angle_vertex_normal_to_travel_direction),2),:)=[];
        small_choice_angles = away_angles./2-away_angle_vertex_normal_to_travel_direction;
        big_choice_angles = away_angles./2+away_angle_vertex_normal_to_travel_direction;
        field_away_normals = [field_away_normals;away_normals];
        field_away_angles = [field_away_angles;away_angles];
        field_away_angle_vertex_normal_to_travel_direction = [field_away_angle_vertex_normal_to_travel_direction;away_angle_vertex_normal_to_travel_direction];
        field_small_choice_angles = [field_small_choice_angles;small_choice_angles];
        field_big_choice_angles = [field_big_choice_angles;big_choice_angles];
        for i=1:length(vertices)-1
            side_vectors(i,1) = vertices(i+1,1)-vertices(i,1);
            side_vectors(i,2) = vertices(i+1,2)-vertices(i,2);
            quiver(vertices(i,1),vertices(i,2),vertices(i+1,1)-vertices(i,1),vertices(i+1,2)-vertices(i,2),'-b');
            % theta = arccos((x dot y)/(mag x * mag y))
            theta_normal_to_side(i) = angle_between_vectors(side_vectors(i,:),vertex_normal_vectors(i,:));
        end
    end
end
function ang = angle_between_vectors(a,b)
    % note this function returns [0,180] so is not directional (i.e. a
    % vector that is rotated 182 deg CCW will return a value of 178 CW
    % this doesn't matter as we expect all vertex normal vectors we care
    % about to be in the direction of the path
    % ang = atan2(norm(cross([a,0],[b,0])), dot(a,b));
    ang = acos(dot(a,b)/(norm(a)*norm(b)));
end
