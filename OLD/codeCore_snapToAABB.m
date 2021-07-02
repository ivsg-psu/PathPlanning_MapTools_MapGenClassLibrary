%% Test case 1: inside, close to bottom
box = [0 0 1 1];
test_point = [0.3 0.2];
fig_num = 2;
snap_point = snapToAABB(box,test_point,fig_num);
assert(isequal([0.3 0],snap_point))

%% Test case 2: inside, close to right
box = [0 0 1 1];
test_point = [0.8 0.4];
fig_num = 2;
snap_point = snapToAABB(box,test_point,fig_num);
assert(isequal([1 0.4],snap_point))

%% Test case 3: inside, close to top
box = [0 0 1 1];
test_point = [0.6 0.9];
fig_num = 2;
snap_point = snapToAABB(box,test_point,fig_num);
assert(isequal([0.6 1],snap_point))

%% Test case 4: inside, close to left
box = [0 0 1 1];
test_point = [0.2 0.7];
fig_num = 2;
snap_point = snapToAABB(box,test_point,fig_num);
assert(isequal([0 0.7],snap_point))

%% Test case 1.1: outside, close to bottom
box = [0 0 1 1];
test_point = [0.3 -0.2];
fig_num = 2;
snap_point = snapToAABB(box,test_point,fig_num);
assert(isequal([0.3 0],snap_point))

%% Test case 2.1: inside, close to right
box = [0 0 1 1];
test_point = [1.2 0.4];
fig_num = 2;
snap_point = snapToAABB(box,test_point,fig_num);
assert(isequal([1 0.4],snap_point))

%% Test case 3.1: inside, close to top
box = [0 0 1 1];
test_point = [0.6 1.1];
fig_num = 2;
snap_point = snapToAABB(box,test_point,fig_num);
assert(isequal([0.6 1],snap_point))

%% Test case 4.1: inside, close to left
box = [0 0 1 1];
test_point = [-0.2 0.7];
fig_num = 2;
snap_point = snapToAABB(box,test_point,fig_num);
assert(isequal([0 0.7],snap_point))

function snap_point = snapToAABB(axis_aligned_bounding_box,test_point,fig_num)

%% main code §
center = [mean([axis_aligned_bounding_box(1,1) axis_aligned_bounding_box(1,3)]),mean([axis_aligned_bounding_box(1,2) axis_aligned_bounding_box(1,4)])];
vector = test_point - center;
angle = atan2(vector(2),vector(1));

snap_point = test_point;
if angle>=-pi/4 && angle<pi/4  % This is the x-max wall
    snap_point(1,1) = axis_aligned_bounding_box(1,3);
elseif angle>=pi/4 && angle<pi*3/4 % This is the y-max wall
    snap_point(1,2) = axis_aligned_bounding_box(1,4);
elseif angle>=-3*pi/4 && angle<(-pi/4) % This is the y-min wall
    snap_point(1,2) = axis_aligned_bounding_box(1,2);
else % This is the x-min wall 
    snap_point(1,1) = axis_aligned_bounding_box(1,1);
end


figure(fig_num);
clf;
hold on;
axis equal
grid on;

% Plot the bounding box
box_outline = [axis_aligned_bounding_box(1,1) axis_aligned_bounding_box(1,2); axis_aligned_bounding_box(1,3) axis_aligned_bounding_box(1,2); axis_aligned_bounding_box(1,3) axis_aligned_bounding_box(1,4); axis_aligned_bounding_box(1,1) axis_aligned_bounding_box(1,4); axis_aligned_bounding_box(1,1) axis_aligned_bounding_box(1,2)];
plot(box_outline(:,1),box_outline(:,2),'-');

% Plot the test point
plot(test_point(:,1),test_point(:,2),'o');

% Plot the snap point
plot(snap_point(:,1),snap_point(:,2),'x');

% Plot the snap point
plot([test_point(:,1) snap_point(1,1)],[test_point(:,2) snap_point(:,2)],'-');
% §
% Debug
%
% Functions §
    
end
