% script_test for fcn_MapGen_isCrossingAABB
% Tests: test for fcn_MapGen_isCrossingAABB
%
% REVISION HISTORY:
% 2024_03_20
% -- first written by S. Harnett
% 2025_04_25 - S. Brennan
% -- changed AABB to show that code is not working as expected

close all; 
AABB = [2.5 2.5 3.5 3.5];
figure; hold on; box on;
fill([AABB(1) AABB(1) AABB(3) AABB(3)], [AABB(2) AABB(4) AABB(4) AABB(2)],[0 0 1],'FaceAlpha',0.3);
test_coords = [1 3 5];
T = combinations(test_coords,test_coords);
test_pts = table2array(T);
[isInside] = fcn_MapGen_isCrossingAABB(AABB, test_pts);
for i = 1:size(isInside,1)
    for j = 1:size(isInside,2)
        if isInside(i,j)
            col_str = 'g';
        else
            col_str = 'r';
        end
        plot([test_pts(i,1) test_pts(j,1)],[test_pts(i,2) test_pts(j,2)],strcat(col_str,'--'),'LineWidth',2)
    end
end
