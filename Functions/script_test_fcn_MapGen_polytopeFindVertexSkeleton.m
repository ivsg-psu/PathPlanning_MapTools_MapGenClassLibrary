% script_test_fcn_MapGen_polytopeFindVertexSkeleton
% Tests function: fcn_MapGen_polytopeFindVertexSkeleton

% REVISION HISTORY:
% 2022_02_15
% -- first written by S. Brennan

fig_num = 675;
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num);

figure(4747);
grid on
grid minor
hold on
axis equal

for cut = 0:0.5:cut_distance(end)
    
    % Find the shape that is less than or equal to the cut
    shape_index = find(cut_distance<=cut,1,'last');
    
    % Grab vertices to start from, cut to start from
    starting_vertices = new_vertices{shape_index};
    starting_cut = cut_distance(shape_index);
    
    % Calculate projection distance
    projection_distance = cut - starting_cut;
    
    % Determine final vertices
    final_vertices = starting_vertices + new_projection_vectors{shape_index}*projection_distance;
    
    % Plot results
    plot(final_vertices(:,1),final_vertices(:,2),'r.-','Linewidth',2,'Markersize',20);
end



%% Example case 1: non-normal wall shrinking
fig_num = 675;
% this polytope has a vertical wall
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num); %#ok<*ASGLU>


%% Example case 2: square
fig_num = 234343;
vertices = [0 0; 1 0; 1 1; 0 1; 0 0]*10;
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num);


% assert that final vertices are within 5% error of having the same x
% position
final_vertices = new_vertices{end}(1,:);
expected_final_vertices = [5 5];
error_tolerance = 0.05;
vertical_error = sum(sum((final_vertices - expected_final_vertices).^2,1).^0.5,2);
% note this error is currently about 33% as of the commit this line was
% added
assert(vertical_error <= error_tolerance,['Expected and actual final points are: ',...
    '\n\t(%d,%d)  (expected) \n\t(%d,%d) (actual) \nyielding a vertical error of %d. Error tolerance was %d.\n'],...
    expected_final_vertices(1,1), expected_final_vertices(1,2),...
    final_vertices(1,1),final_vertices(1,2),vertical_error,error_tolerance);

%% Example case 3: wide rectangle
fig_num = 2464;
vertices = [0 0; 1 0; 1 0.5; 0 0.5; 0 0]*10;
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num);


% assert that final vertices are within 5% error of having the same x
% position
final_vertices = new_vertices{end}(1,:);
expected_final_vertices = [5 2.5];
error_tolerance = 0.05;
vertical_error = sum(sum((final_vertices - expected_final_vertices).^2,1).^0.5,2);
% note this error is currently about 33% as of the commit this line was
% added
assert(vertical_error <= error_tolerance,['Expected and actual final points are: ',...
    '\n\t(%d,%d)  (expected) \n\t(%d,%d) (actual) \nyielding a vertical error of %d. Error tolerance was %d.\n'],...
    expected_final_vertices(1,1), expected_final_vertices(1,2),...
    final_vertices(1,1),final_vertices(1,2),vertical_error,error_tolerance);


%% Example case 4: tall rectangle
fig_num = 2465;
vertices = [0 0; 0.5 0; 0.5 1; 0 1; 0 0]*10;
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num);

% assert that final vertices are within 5% error of having the same x
% position
final_vertices = new_vertices{end}(1,:);
expected_final_vertices = [2.5 5];
error_tolerance = 0.05;
vertical_error = sum(sum((final_vertices - expected_final_vertices).^2,1).^0.5,2);
% note this error is currently about 33% as of the commit this line was
% added
assert(vertical_error <= error_tolerance,['Expected and actual final points are: ',...
    '\n\t(%d,%d)  (expected) \n\t(%d,%d) (actual) \nyielding a vertical error of %d. Error tolerance was %d.\n'],...
    expected_final_vertices(1,1), expected_final_vertices(1,2),...
    final_vertices(1,1),final_vertices(1,2),vertical_error,error_tolerance);


%% Example case 5: goofy polytope
fig_num = 47464;
vertices = [0 0; 10 0; 5 15; 4 17; 1 13; 0 5; 0 0];
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num);



%% Random polytope calculation
% Set up polytopes
close all;
polytopes = fcn_MapGen_haltonVoronoiTiling([1 100],[1 1]);

edge_cut_step = 0.002;
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box);

% Pick a random polytope
Npolys = length(trim_polytopes);
rand_poly = 1+floor(rand*Npolys);
shrinker = trim_polytopes(rand_poly);

% Do skeleton calculation
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(shrinker.vertices,2727);


fig_num = 11;
figure(fig_num); 
clf;
axis equal;
hold on;

% for edge_cut = edge_cut_step:edge_cut_step:(shrinker.max_radius+edge_cut_step)
for cut = edge_cut_step:edge_cut_step:(shrinker.max_radius+edge_cut_step)
    
    % Find the shape that is less than or equal to the cut
    shape_index = find(cut_distance<=cut,1,'last');
    
    % Grab vertices to start from, cut to start from
    starting_vertices = new_vertices{shape_index};
    starting_cut = cut_distance(shape_index);
    
    % Calculate projection distance
    projection_distance = cut - starting_cut;
    
    % Determine final vertices
    final_vertices = starting_vertices + new_projection_vectors{shape_index}*projection_distance;
    
    % Plot results
    plot(final_vertices(:,1),final_vertices(:,2),'.-','Linewidth',2,'Markersize',20);
end

