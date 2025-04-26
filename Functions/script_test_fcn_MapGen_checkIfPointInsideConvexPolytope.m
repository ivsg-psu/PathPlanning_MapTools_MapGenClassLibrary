% script_test_fcn_MapGen_checkIfPointInsideConvexPolytope
% Tests: fcn_MapGen_checkIfPointInsideConvexPolytope

% 
% REVISION HISTORY:
% 
% 2021_07_14 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§


close all;

%% Test case 1: inside, close to bottom
vertices = [0 0; 1 0; 1 1; 0 1; 0 0];
test_point = [0.3 0.2];
fig_num = 1;
in_polytope = fcn_MapGen_checkIfPointInsideConvexPolytope(test_point,vertices,fig_num);
assert(isequal(true,in_polytope))

%% Test case 2: On the vertex
vertices = [0 0; 1 0; 1 1; 0 1; 0 0];
test_point = [0 0];
fig_num = 2;
in_polytope = fcn_MapGen_checkIfPointInsideConvexPolytope(test_point,vertices,fig_num);
assert(isequal(true,in_polytope))

%% Test case 3: On the edge
vertices = [0 0; 1 0; 1 1; 0 1; 0 0];
test_point = [0 0.5];
fig_num = 3;
in_polytope = fcn_MapGen_checkIfPointInsideConvexPolytope(test_point,vertices,fig_num);
assert(isequal(true,in_polytope))

%% Test case 4: outside
vertices = [0 0; 1 0; 1 1; 0 1; 0 0];
test_point = [-1 -1];
fig_num = 4;
in_polytope = fcn_MapGen_checkIfPointInsideConvexPolytope(test_point,vertices,fig_num);
assert(isequal(false,in_polytope))

%% Test case 5: lots of random points
vertices = [0 0; 1 0; 1 1; 0 1; 0 0];
test_points = rand(100,2)*2 - [0.5 0.5];
fig_num = 5;
for ith_point = 1:length(test_points)
    test_point = test_points(ith_point,:);
    in_polytope = fcn_MapGen_checkIfPointInsideConvexPolytope(test_point,vertices,fig_num);
end



