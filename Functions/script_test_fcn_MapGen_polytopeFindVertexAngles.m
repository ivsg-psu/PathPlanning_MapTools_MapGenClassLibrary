% script_test_fcn_MapGen_polytopeFindVertexAngles
% Tests function: fcn_MapGen_polytopeFindVertexAngles

% REVISION HISTORY:
% 2021_08_01
% -- first written by S. Brennan 


%% Basic example of vertex calculation - a square
fig_num = 1;
vertices = [0 0; 1 0; 1 1; 0 1; 0 0];
angles =...
    fcn_MapGen_polytopeFindVertexAngles(...
    vertices,fig_num);
assert(1000*eps>abs(360-sum(angles)*180/pi));

%% Basic example of vertex calculation - a triangle
fig_num = 2;
vertices = [0 0; 1 1; 0 1; 0 0];
angles =...
    fcn_MapGen_polytopeFindVertexAngles(...
    vertices,fig_num);
assert(1000*eps>abs(360-sum(angles)*180/pi));

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
angles =...
    fcn_MapGen_polytopeFindVertexAngles(...
    shrinker.vertices,fig_num);
assert(1000*eps>abs(360-sum(angles)*180/pi));


  