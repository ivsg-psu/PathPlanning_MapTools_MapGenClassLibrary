% script_test_fcn_MapGen_cropVerticesByWallIntersections
% Tests: fcn_MapGen_cropVerticesByWallIntersections

% 
% REVISION HISTORY:
% 
% 2021_07_15 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§




%% Test case 1: Going from out to in
walls = [0 0; 1 0; 1 1; 0 1; 0 0];
test_vertices = [-0.3 0.2; 0.3 0.2];
fig_num = 1;
[cropped_vertices,NwallsHit] = fcn_MapGen_cropVerticesByWallIntersections(test_vertices,walls,fig_num);
assert(all(([0 0.2; 0.3 0.2]-eps*[1 1; 1 1])<cropped_vertices,'all') && all([0 0.2; 0.3 0.2]+eps*[1 1; 1 1]>cropped_vertices,'all'));
assert(isequal(NwallsHit,1))


%% Test case 2: Going from in to out
walls = [0 0; 1 0; 1 1; 0 1; 0 0];
test_vertices = [0.3 0.2; 1.3 0.2];
fig_num = 2;
[cropped_vertices,NwallsHit] = fcn_MapGen_cropVerticesByWallIntersections(test_vertices,walls,fig_num);
assert(all(([0.3 0.2; 1 0.2]-eps*[1 1; 1 1])<cropped_vertices,'all') && all([0.3 0.2; 1 0.2]+eps*[1 1; 1 1]>cropped_vertices,'all'));
assert(isequal(NwallsHit,1))

%% Test case 3: all inside
walls = [0 0; 1 0; 1 1; 0 1; 0 0];
test_vertices = [0.3 0.2; 0.4 0.2];
fig_num = 3;
[cropped_vertices,NwallsHit] = fcn_MapGen_cropVerticesByWallIntersections(test_vertices,walls,fig_num);
assert(all(([0.3 0.2; 0.4 0.2]-eps*[1 1; 1 1])<cropped_vertices,'all') && all([0.3 0.2; 0.4 0.2]+eps*[1 1; 1 1]>cropped_vertices,'all'));
assert(isequal(NwallsHit,0))

%% Test case 4: all outside
walls = [0 0; 1 0; 1 1; 0 1; 0 0];
test_vertices = [-0.3 0.2; -0.4 0.2];
fig_num = 4;
[cropped_vertices,NwallsHit] = fcn_MapGen_cropVerticesByWallIntersections(test_vertices,walls,fig_num);
assert(isequal([],cropped_vertices));
assert(isequal(NwallsHit,0))

%% Test case 5: crossing over
walls = [0 0; 1 0; 1 1; 0 1; 0 0];
test_vertices = [-0.3 0.2; 1.4 0.2];
fig_num = 5;
[cropped_vertices,NwallsHit] = fcn_MapGen_cropVerticesByWallIntersections(test_vertices,walls,fig_num);
true_answer = [1 0.2; 0 0.2];
assert(all((true_answer-eps*[1 1; 1 1])<cropped_vertices,'all') && all(true_answer+eps*[1 1; 1 1]>cropped_vertices,'all'));
assert(isequal(NwallsHit,2))

%% Test case 6: aligned with edge, across
walls = [0 0; 1 0; 1 1; 0 1; 0 0];
test_vertices = [-0.3 0; 1.4 0];
fig_num = 6;
[cropped_vertices,NwallsHit] = fcn_MapGen_cropVerticesByWallIntersections(test_vertices,walls,fig_num);
true_answer = [0 0; 1 0];
assert(all((true_answer-eps*[1 1; 1 1])<cropped_vertices,'all') && all(true_answer+eps*[1 1; 1 1]>cropped_vertices,'all'));
assert(isequal(NwallsHit,3))


%% Test case 7: aligned with edge, inside
walls = [0 0; 1 0; 1 1; 0 1; 0 0];
test_vertices = [0.5 0; 1.4 0];
fig_num = 7;
[cropped_vertices,NwallsHit] = fcn_MapGen_cropVerticesByWallIntersections(test_vertices,walls,fig_num);
true_answer = [0.5 0; 1 0];
assert(all((true_answer-eps*[1 1; 1 1])<cropped_vertices,'all') && all(true_answer+eps*[1 1; 1 1]>cropped_vertices,'all'));
assert(isequal(NwallsHit,2))

%% Test case 8: aligned with edge, corner
walls = [0 0; 1 0; 1 1; 0 1; 0 0];
test_vertices = [1 0; 1.4 0];
fig_num = 8;
[cropped_vertices,NwallsHit] = fcn_MapGen_cropVerticesByWallIntersections(test_vertices,walls,fig_num);
true_answer = [1 0];
assert(all((true_answer-eps*[1 1; 1 1])<cropped_vertices,'all') && all(true_answer+eps*[1 1; 1 1]>cropped_vertices,'all'));
assert(isequal(NwallsHit,2))
