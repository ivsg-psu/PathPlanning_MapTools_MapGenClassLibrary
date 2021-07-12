% script_test_fcn_MapGen_cropPolytopeToRange

%% Test case 1: simple square
fig_num = 1;

verticies = [1 1; 1 -1; -1 -1; -1 1]*0.5;
interior_point = [0 0];
AABB = [0 0 1 1];
[cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies,interior_point,AABB,fig_num);


%% Test case 2: simple triangle
fig_num = 2;

verticies = [1 -0.5; -1 -1; -0.5 1];
interior_point = [0 0];
[cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies,interior_point,AABB,fig_num);


%% Test case 3: vertices enclosing region
fig_num = 3;

verticies = [1 1; 1 -1; -1 -1; -1 1]*2;
interior_point = [0 0];
[cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies,interior_point,AABB,fig_num);


%% Test case 4: vertices all within region
fig_num = 4;

verticies = [ 0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
interior_point = [0.5 0.5];
[cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies,interior_point,AABB,fig_num);

%% Test case 5: infinity in the numbers
fig_num = 6;

verticies = [...
   -0.0457    0.0471
   -0.4075    0.0851
   -2.7148    0.1857
  -11.3670    0.2644
       Inf       Inf
    0.3841  -53.6209
    0.2255   -6.8725
    0.1203   -1.0609
    0.0352   -0.0163
   -0.0081    0.0368];

interior_point = [0 0];
[cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies,interior_point,AABB,fig_num);


%% Points on edges
fig_num = 7;

verticies = [
   0.951225451552411   0.038481963258910
   0.951225451552411                   0
   1.000000000000000   0.055571820180916
   0.972711461970130   0.055571820180916
   0.967138332311524   0.053930059059772];
interior_point = [0.978515625000000   0.035665294924554];
[cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies,interior_point,AABB,fig_num);


%% Points valid except infinity
fig_num = 8;

interior_point = [0.9404    0.0133];

verticies = [
    0.9512    0.0385
       Inf       Inf
    0.9275    0.0315
    0.9318    0.0417];

[cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies,interior_point,AABB,fig_num);

