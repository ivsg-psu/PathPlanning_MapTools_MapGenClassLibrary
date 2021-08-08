% script_test_fcn_MapGen_polytopeShrinkFromEdges
% Tests function: fcn_MapGen_polytopeShrinkFromEdges

% REVISION HISTORY:
% 2021_08_02
% -- first written by S. Brennan using


%% Basic example of vertex calculation - a square
fig_num = 1;
vertices = [0 0; 1 0; 1 1; 0 1; 0 0]*5;
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
