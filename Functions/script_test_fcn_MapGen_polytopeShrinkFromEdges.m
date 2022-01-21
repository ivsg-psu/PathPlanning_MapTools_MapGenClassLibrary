% script_test_fcn_MapGen_polytopeShrinkFromEdges
% Tests function: fcn_MapGen_polytopeShrinkFromEdges

% REVISION HISTORY:
% 2021_08_02
% -- first written by S. Brennan using


%% failing case 1: non-normal wall shrinking
fig_num = 675;
% this polytope has a vertical wall
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
% assert that the wall is vertical to start (i.e. x position of 1st and 4th
% vertices are equal)
assert(vertices(1,1) == vertices(4,1));
test_polytope.vertices = vertices;
% fill in other fields from the vertices field
test_polytope = fcn_MapGen_fillPolytopeFieldsFromVertices(test_polytope);
% perform a small edge shrink
edge_cut = 0.1;
shrunk = fcn_MapGen_polytopeShrinkFromEdges(...
    test_polytope,edge_cut,fig_num);
% extract new vertices
new_vertices = shrunk.vertices;
% assert that new vertices are within 5% error of having the same x
% position
error_tolerance = 0.05;
vertical_error = abs(new_vertices(1,1)-new_vertices(4,1))/new_vertices(4,1);
% note this error is currently about 33% as of the commit this line was
% added
assert(vertical_error <= error_tolerance,['Wall should be vertical but x positions of start and end point',...
    ' were %d and %d yielding a vertical error of %d. Error tolerance was %d.\n'],...
    new_vertices(1,1),new_vertices(4,1),vertical_error,error_tolerance);

%% failing case 2: concave polytope produced
fig_num = 657;
% this polytope is convex
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
test_polytope.vertices = vertices;
% fill in other fields from the vertices field
test_polytope = fcn_MapGen_fillPolytopeFieldsFromVertices(test_polytope);
% assert that the polytope is convex to start
[angles, ~, ~] = fcn_MapGen_polytopeFindVertexAngles(test_polytope.vertices);
interior_angles = 180-angles*180/pi
assert(~any(interior_angles>180));
% perform a large edge shrink
edge_cut = 3;
shrunk = fcn_MapGen_polytopeShrinkFromEdges(...
    test_polytope,edge_cut,fig_num);
% extract new vertices
new_vertices = shrunk.vertices;
% assert that the polytope is convex after shrinking
[new_angles, ~, ~] = fcn_MapGen_polytopeFindVertexAngles(new_vertices);
new_interior_angles = 180-new_angles*180/pi
assert(~any(new_interior_angles>180),['All interior angles must be < 180 ',...
    'polytope to be convex']);

%% failing case 3: NaN angles produced
fig_num = 756;
% this polytope is convex
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
test_polytope.vertices = vertices;
% fill in other fields from the vertices field
test_polytope = fcn_MapGen_fillPolytopeFieldsFromVertices(test_polytope);
% assert that the polytope is convex to start
[angles, ~, ~] = fcn_MapGen_polytopeFindVertexAngles(test_polytope.vertices);
interior_angles = 180-angles*180/pi
assert(~any(interior_angles>180));
% perform an even larger edge shrink
edge_cut = 4;
shrunk = fcn_MapGen_polytopeShrinkFromEdges(...
    test_polytope,edge_cut,fig_num);
% extract new vertices
new_vertices = shrunk.vertices;
% assert that the polytope is convex after shrinking
[new_angles, ~, ~] = fcn_MapGen_polytopeFindVertexAngles(new_vertices);
new_interior_angles = 180-new_angles*180/pi
assert(~any(isnan(new_interior_angles)),['Interior angles had NaN values.']);

%% Basic example of vertex calculation - a square
fig_num = 1;
% vertices = [0 0; 1 0; 1 1; 0 1; 0 0]*5;
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
test_polytope.vertices = vertices;

% fill in other fields from the vertices field
test_polytope = fcn_MapGen_fillPolytopeFieldsFromVertices(test_polytope);


edge_cut = 0.1;
fcn_MapGen_polytopeShrinkFromEdges(...
    test_polytope,edge_cut,fig_num);

%% Basic example of vertex calculation - a triangle
close all;

fig_num = 2;
vertices = [0 0; 1 1; 0 1; 0 0]*5;
test_polytope.vertices = vertices;

% fill in other fields from the vertices field
test_polytope = fcn_MapGen_fillPolytopeFieldsFromVertices(test_polytope);

edge_cut = 0.1;
fcn_MapGen_polytopeShrinkFromEdges(...
    test_polytope,edge_cut,fig_num);

%% Basic example of vertex calculation - a triangle with too big a cut
fig_num = 3;
vertices = [0 0; 1 1; 0 1; 0 0];
test_polytope.vertices = vertices;

% fill in other fields from the vertices field
test_polytope = fcn_MapGen_fillPolytopeFieldsFromVertices(test_polytope);

edge_cut = 2;
fcn_MapGen_polytopeShrinkFromEdges(...
    test_polytope,edge_cut,fig_num);

%% Iterations over cuts
close all;
fig_num = 4;
figure(fig_num); clf;

vertices = [0 0; 0.4 0.1; 1 1; 0 1; 0 0]*5;
test_polytope.vertices = vertices;

% fill in other fields from the vertices field
test_polytope = fcn_MapGen_fillPolytopeFieldsFromVertices(test_polytope);

for edge_cut = 0.02:0.02:1
    fcn_MapGen_polytopeShrinkFromEdges(...
        test_polytope,edge_cut,fig_num);
end

%% Random polytope calculation
% Set up polytopes
polytopes = fcn_MapGen_haltonVoronoiTiling([1 100],[1 1]);

edge_cut_step = 0.002;
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box);

% Pick a random polytope
Npolys = length(trim_polytopes);
rand_poly = 1+floor(rand*Npolys);
shrinker = trim_polytopes(rand_poly);

%% Basic example of vertex calculation
fig_num = 11;
figure(fig_num); clf;
% for edge_cut = edge_cut_step:edge_cut_step:(shrinker.max_radius+edge_cut_step)
for edge_cut = edge_cut_step:edge_cut_step:(shrinker.max_radius+edge_cut_step)
    fcn_MapGen_polytopeShrinkFromEdges(...
        shrinker,edge_cut,fig_num);
end
