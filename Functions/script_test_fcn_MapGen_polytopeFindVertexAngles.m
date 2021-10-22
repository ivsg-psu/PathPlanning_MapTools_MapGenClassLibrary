% script_test_fcn_MapGen_polytopeFindVertexAngles
% Tests function: fcn_MapGen_polytopeFindVertexAngles

% REVISION HISTORY:
% 2021_08_01
% -- first written by S. Brennan


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
vertex_normal_vectors = [unit_direction_of_cut(:,1)-vertices(1:end-1,1),unit_direction_of_cut(:,2)-vertices(1:end-1,2)]
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
vertex_normal_vectors = [unit_direction_of_cut(:,1)-vertices(1:end-1,1),unit_direction_of_cut(:,2)-vertices(1:end-1,2)]
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
vertex_normal_vectors = circshift(vertex_normal_vectors,1)
side_vectors = zeros(length(vertices)-1,2);
theta_normal_to_side = zeros(length(vertices)-1);
for i=1:length(vertices)-1
    side_vectors(i,1) = vertices(i+1,1)-vertices(i,1);
    side_vectors(i,2) = vertices(i+1,2)-vertices(i,2);
    quiver(vertices(i,1),vertices(i,2),vertices(i+1,1)-vertices(i,1),vertices(i+1,2)-vertices(i,2),'-b');
    % theta = arccos((x dot y)/(mag x * mag y))
    theta_normal_to_side(i) = rad2deg(acos((dot([side_vectors(i,1),side_vectors(i,2)],[vertex_normal_vectors(i,1),vertex_normal_vectors(i,2)]))/(norm([side_vectors(i,1),side_vectors(i,2)])*norm([vertex_normal_vectors(i,1),vertex_normal_vectors(i,2)]))));
end

