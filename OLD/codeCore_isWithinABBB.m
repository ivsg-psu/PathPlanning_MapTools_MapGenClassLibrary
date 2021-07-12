AABB = [0 0 1 1];
test_points = randn(100,2);
fig_num = 1;
isInside = fcn_MapGen_isWithinABBB(AABB,test_points,fig_num);



function isInside = fcn_MapGen_isWithinABBB(AABB,test_points,fig_num)
% Checks if the points are within the AABB, returning a vector of 1's or
% 0's the same length as the nubmer of rows of points.

%% main code §

% % See: https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection
% % for details on axis-aligned bounding boxes (AABB)

isInside = (test_points(:,1)>AABB(1,1))  & ...
    (test_points(:,2)>AABB(1,2))  & ...
    (test_points(:,1)<AABB(1,3))  & ...
    (test_points(:,2)<AABB(1,4));

% plotting
figure(fig_num);
clf;
hold on;

% Convert axis-aligned bounding box to wall format
walls = [AABB(1,1) AABB(1,2); AABB(1,3) AABB(1,2); AABB(1,3) AABB(1,4); AABB(1,1) AABB(1,4); AABB(1,1) AABB(1,2)];

% Plot the walls
plot(walls(:,1),walls(:,2),'k-');

% Plot the test_points

% plot(...
%     [test_points(:,1); test_points(1,1)],...
%     [test_points(:,2); test_points(1,2)],...
%     '.-');
plot(test_points(:,1), test_points(:,2),'k.');

% Plot the interior points
plot(test_points(isInside,1),test_points(isInside,2),'go');

% §
% Debug
%
% Functions §

end
