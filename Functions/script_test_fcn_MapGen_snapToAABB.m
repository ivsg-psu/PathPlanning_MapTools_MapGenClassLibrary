% script_test_fcn_MapGen_snapToAABB
% Tests: fcn_MapGen_snapToAABB

% 
% REVISION HISTORY:
% 
% 2021_07_02 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§




%% Test case 1: inside, close to bottom
box = [0 0 1 1];
test_point = [0.3 0.2];
fig_num = 2;
snap_type = 0; 
snap_point = fcn_MapGen_snapToAABB(box,test_point,snap_type, fig_num); 
assert(isequal([0.3 0],snap_point))

%% Test case 2: inside, close to right
box = [0 0 1 1];
test_point = [0.8 0.4];
fig_num = 2;
snap_type = 0; 
snap_point = fcn_MapGen_snapToAABB(box,test_point,snap_type, fig_num); 
assert(isequal([1 0.4],snap_point))

%% Test case 3: inside, close to top
box = [0 0 1 1];
test_point = [0.6 0.9];
fig_num = 2;
snap_type = 0; 
snap_point = fcn_MapGen_snapToAABB(box,test_point,snap_type, fig_num); 
assert(isequal([0.6 1],snap_point))

%% Test case 4: inside, close to left
box = [0 0 1 1];
test_point = [0.2 0.7];
fig_num = 2;
snap_type = 0; 
snap_point = fcn_MapGen_snapToAABB(box,test_point,snap_type, fig_num); 
assert(isequal([0 0.7],snap_point))

%% Test case 1.1: outside, close to bottom
box = [0 0 1 1];
test_point = [0.3 -0.2];
fig_num = 2;
snap_type = 0; 
snap_point = fcn_MapGen_snapToAABB(box,test_point,snap_type, fig_num); 
% assert(isequal([0.3 0],snap_point))
assert(isequal(test_point,snap_point))


%% Test case 2.1: inside, close to right
box = [0 0 1 1];
test_point = [1.2 0.4];
fig_num = 2;
snap_type = 0; 
snap_point = fcn_MapGen_snapToAABB(box,test_point,snap_type, fig_num); 
% assert(isequal([1 0.4],snap_point))
assert(isequal(test_point,snap_point))

%% Test case 3.1: inside, close to top
box = [0 0 1 1];
test_point = [0.6 1.1];
fig_num = 2;
snap_type = 0; 
snap_point = fcn_MapGen_snapToAABB(box,test_point,snap_type, fig_num); 
% assert(isequal([0.6 1],snap_point))
assert(isequal(test_point,snap_point))

%% Test case 4.1: inside, close to left
box = [0 0 1 1];
test_point = [-0.2 0.7];
fig_num = 2;
snap_type = 0; 
snap_point = fcn_MapGen_snapToAABB(box,test_point,snap_type, fig_num); 
% assert(isequal([0 0.7],snap_point))
assert(isequal(test_point,snap_point))


%% Test case 5.1: user-defined vector
box = [0 0 1 1];
test_point = [0.2 0.7; 0.2 0.9];
fig_num = 51;
snap_type = 2;
snap_point = fcn_MapGen_snapToAABB(box,test_point,snap_type, fig_num);
true_point = [0.2 1];
assert(all(abs(true_point-snap_point)<eps))

