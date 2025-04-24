% script_test_fcn_MapGen_snapToAABB
% Tests: fcn_MapGen_snapToAABB

%
% REVISION HISTORY:
%
% 2021_07_14 by Sean Brennan
% -- first write of script
% 2025_04_24 by Sean Brennan
% -- fixed calls to match revised function
%%%%%%%%%%%%%%§

close all;

%% check input arguments?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                     _______                    ___
%  / ____|                   |__   __|                  / _ \
% | (___  _ __   __ _ _ __      | |_   _ _ __   ___    | | | |
%  \___ \| '_ \ / _` | '_ \     | | | | | '_ \ / _ \   | | | |
%  ____) | | | | (_| | |_) |    | | |_| | |_) |  __/   | |_| |
% |_____/|_| |_|\__,_| .__/     |_|\__, | .__/ \___|    \___/
%                    | |            __/ | |
%                    |_|           |___/|_|

% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Snap%20Type%20%20%200
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Test case 0.1: inside, close to bottom
fig_num = 001;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [0.3 0.2];
snap_type = 0;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([0.1667         0],round(snap_point,4)))
assert(isequal(1,wall_number))


%% Test case 0.2: inside, close to right
fig_num = 002;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [0.8 0.4];
snap_type = 0;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([1.0000    0.3333],round(snap_point,4)))
assert(isequal(2,wall_number))

%% Test case 0.3: inside, close to top
fig_num = 003;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [0.6 0.9];
snap_type = 0;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([0.6250    1.0000],round(snap_point,4)))
assert(isequal(3,wall_number))

%% Test case 0.4: inside, close to left
fig_num = 004;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [0.2 0.7];
snap_type = 0;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([0    0.8333],round(snap_point,4)))
assert(isequal(4,wall_number))


%% Test case 0.5: outside, close to bottom
fig_num = 005;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [0.3 -0.2];
snap_type = 0;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([0.3000   -0.2000],round(snap_point,4)))
assert(isequal(1,wall_number))



%% check snap type 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                     _______                   __
%  / ____|                   |__   __|                 /_ |
% | (___  _ __   __ _ _ __      | |_   _ _ __   ___     | |
%  \___ \| '_ \ / _` | '_ \     | | | | | '_ \ / _ \    | |
%  ____) | | | | (_| | |_) |    | | |_| | |_) |  __/    | |
% |_____/|_| |_|\__,_| .__/     |_|\__, | .__/ \___|    |_|
%                    | |            __/ | |
%                    |_|           |___/|_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Snap%20Type%20%20%201
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Test case 1.1: inside, close to bottom
fig_num = 101;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [0.3 0.2];
snap_type = 1;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([0.3         0],round(snap_point,4)))
assert(isequal(1,wall_number))


%% Test case 1.2: inside, close to right
fig_num = 102;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [0.8 0.4];
snap_type = 1;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([1.0000    0.4],round(snap_point,4)))
assert(isequal(2,wall_number))

%% Test case 1.3: inside, close to top
fig_num = 103;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [0.6 0.9];
snap_type = 1;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([0.6    1.0000],round(snap_point,4)))
assert(isequal(3,wall_number))

%% Test case 1.4: inside, close to left
fig_num = 104;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [0.2 0.7];
snap_type = 1;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([0    0.7],round(snap_point,4)))
assert(isequal(4,wall_number))


%% Test case 1.5: outside, close to bottom
fig_num = 105;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [0.3 -0.2];
snap_type = 1;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([0.3000   -0.2000],round(snap_point,4)))
assert(isequal(1,wall_number))

%% Test case 1.6: inside, close to right
fig_num = 106;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [1.2 0.4];
snap_type = 1;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([1.2000    0.4000],round(snap_point,4)))
assert(isequal(2,wall_number))

%% Test case 1.7: inside, close to top
fig_num = 107;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [0.6 1.1];
snap_type = 1;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([0.6000    1.1000],round(snap_point,4)))
assert(isequal(3,wall_number))


%% Test case 1.8: inside, close to left
fig_num = 108;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [-0.2 0.7];
snap_type = 1;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([-0.2000    0.7000],round(snap_point,4)))
assert(isequal(4,wall_number))


%% check snap type 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                     _______                   ___
%  / ____|                   |__   __|                 |__ \
% | (___  _ __   __ _ _ __      | |_   _ _ __   ___       ) |
%  \___ \| '_ \ / _` | '_ \     | | | | | '_ \ / _ \     / /
%  ____) | | | | (_| | |_) |    | | |_| | |_) |  __/    / /_
% |_____/|_| |_|\__,_| .__/     |_|\__, | .__/ \___|   |____|
%                    | |            __/ | |
%                    |_|           |___/|_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Snap%20Type%20%20%202
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Test case 2.9: user-defined vector
fig_num = 209;

axis_aligned_bounding_box = [0 0 1 1];
test_point = [0.2 0.7; 0.2 0.9];
snap_type = 2;

[snap_point, wall_number] = fcn_MapGen_snapToAABB(axis_aligned_bounding_box, test_point, (snap_type), (fig_num) );

assert(isequal([0.2000    1.0000],round(snap_point,4)))
assert(isequal(3,wall_number))


