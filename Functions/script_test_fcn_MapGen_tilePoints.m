% script_test_fcn_MapGen_tilePoints
% Tests function: fcn_MapGen_tilePoints

% REVISION HISTORY:
% 2023_02_24
% -- first written by S. Brennan
 close all
 

%% Basic example call
fig_num = 1;
Npoints = 20;
tile_depth = 1;
AABB = [0 0 2 2];
input_points = 2*rand(Npoints,2);
[tiled_points] = fcn_MapGen_tilePoints(input_points,tile_depth,AABB, fig_num);

% Check sizes
assert(length(tiled_points(1,:))==2); % Is it 2 columns?
assert(length(tiled_points(:,1))==Npoints*9); % Does it have right number of points?



% assert(vertical_error <= error_tolerance,['Wall should be vertical but x positions of start and end point',...
%     ' were %d and %d yielding a vertical error of %d. Error tolerance was %d.\n'],...
%     new_vertices(1,1),new_vertices(4,1),vertical_error,error_tolerance);

%% Call with points outside AABB
fig_num = 1;
Npoints = 200;
tile_depth = 1;
AABB = [0 0 2 2];
input_points = 0.5*randn(Npoints,2)+[7,0] + 1;
[tiled_points] = fcn_MapGen_tilePoints(input_points,tile_depth,AABB, fig_num);


% Check sizes
assert(length(tiled_points(1,:))==2); % Is it 2 columns?
assert(length(tiled_points(:,1))==Npoints*9); % Does it have right number of points?

% assert(vertical_error <= error_tolerance,['Wall should be vertical but x positions of start and end point',...
%     ' were %d and %d yielding a vertical error of %d. Error tolerance was %d.\n'],...
%     new_vertices(1,1),new_vertices(4,1),vertical_error,error_tolerance);

%% Call with AABB bigger than points
fig_num = 1;
Npoints = 200;
tile_depth = 1;
AABB = [0 0 5 5];
input_points = 0.5*randn(Npoints,2)+1;
[tiled_points] = fcn_MapGen_tilePoints(input_points,tile_depth,AABB, fig_num);

% Check sizes
assert(length(tiled_points(1,:))==2); % Is it 2 columns?
assert(length(tiled_points(:,1))==Npoints*9); % Does it have right number of points?

% assert(vertical_error <= error_tolerance,['Wall should be vertical but x positions of start and end point',...
%     ' were %d and %d yielding a vertical error of %d. Error tolerance was %d.\n'],...
%     new_vertices(1,1),new_vertices(4,1),vertical_error,error_tolerance);

%% Call with AABB shifted
fig_num = 1;
Npoints = 200;
tile_depth = 1;
AABB = [0 0 2 2]- 7;
input_points = 0.5*randn(Npoints,2) + 1;
[tiled_points] = fcn_MapGen_tilePoints(input_points,tile_depth,AABB, fig_num);


% Check sizes
assert(length(tiled_points(1,:))==2); % Is it 2 columns?
assert(length(tiled_points(:,1))==Npoints*9); % Does it have right number of points?

% assert(vertical_error <= error_tolerance,['Wall should be vertical but x positions of start and end point',...
%     ' were %d and %d yielding a vertical error of %d. Error tolerance was %d.\n'],...
%     new_vertices(1,1),new_vertices(4,1),vertical_error,error_tolerance);

%% Call with depth of 2
fig_num = 1;
Npoints = 20;
tile_depth = 2;
AABB = [0 0 2 2];
input_points = 2*rand(Npoints,2);
[tiled_points] = fcn_MapGen_tilePoints(input_points,tile_depth,AABB, fig_num);


% Check sizes
assert(length(tiled_points(1,:))==2); % Is it 2 columns?
assert(length(tiled_points(:,1))==Npoints*25); % Does it have right number of points?

% assert(vertical_error <= error_tolerance,['Wall should be vertical but x positions of start and end point',...
%     ' were %d and %d yielding a vertical error of %d. Error tolerance was %d.\n'],...
%     new_vertices(1,1),new_vertices(4,1),vertical_error,error_tolerance);

%% Call with depth of 5
fig_num = 1;
Npoints = 20;
tile_depth = 5;
AABB = [0 0 2 2];
input_points = 2*rand(Npoints,2);
[tiled_points] = fcn_MapGen_tilePoints(input_points,tile_depth,AABB, fig_num);

% Check sizes
assert(length(tiled_points(1,:))==2); % Is it 2 columns?
assert(length(tiled_points(:,1))==Npoints*121); % Does it have right number of points?


% assert(vertical_error <= error_tolerance,['Wall should be vertical but x positions of start and end point',...
%     ' were %d and %d yielding a vertical error of %d. Error tolerance was %d.\n'],...
%     new_vertices(1,1),new_vertices(4,1),vertical_error,error_tolerance);

if 1==0 % Fail cases
    %% Input incorrect, wrong number of arguments
    ;
    input_points = rand(10,1);
    tile_depth = 2;
    fcn_MapGen_tilePoints(input_points,tile_depth);

    %% Input incorrect, wrong number of arguments
    
    input_points = rand(10,1);
    tile_depth = 2;
    AABB = [0 0 1 1];
    fcn_MapGen_tilePoints(input_points,tile_depth,AABB,3,4);

    %% Input incorrect for input_points, only 1 column
    
    input_points = rand(10,1);
    tile_depth = 2;
    AABB = [0 0 1 1];
    fcn_MapGen_tilePoints(input_points,tile_depth,AABB);
    %% Input incorrect for input_points, 3 column
    
    input_points = rand(10,1);
    tile_depth = 2;
    AABB = [0 0 1 1];
    fcn_MapGen_tilePoints(input_points,tile_depth,AABB);
    %% Input incorrect for tile_depth, not strictly positive
    
    input_points = rand(10,2);
    tile_depth = 0;
    AABB = [0 0 1 1];
    fcn_MapGen_tilePoints(input_points,tile_depth,AABB);
    %% Input incorrect for tile_depth, not an integer
    
    input_points = rand(10,2);
    tile_depth = 1.2;
    AABB = [0 0 1 1];
    fcn_MapGen_tilePoints(input_points,tile_depth,AABB);
    %% Input incorrect for AABB, not an 4x1
    
    input_points = rand(10,2);
    tile_depth = 1;
    AABB = [0 0 1];
    fcn_MapGen_tilePoints(input_points,tile_depth,AABB);
    %% Input incorrect for AABB, not an 4x1
    
    input_points = rand(10,2);
    tile_depth = 1;
    AABB = [0 0 1 1 1];
    fcn_MapGen_tilePoints(input_points,tile_depth,AABB);
    %% Input incorrect for AABB, not an 4x1
    
    input_points = rand(10,2);
    tile_depth = 1;
    AABB = [0 0 1 1; 0 0 0 0];
    fcn_MapGen_tilePoints(input_points,tile_depth,AABB);
end
