% script_test_fcn_MapGen_isWithinABBB
% Tests: fcn_MapGen_isWithinABBB

%
% REVISION HISTORY:
%
% 2021_07_11 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§



AABB = [0 0 1 1];
test_points = randn(100,2);
fig_num = 1;
isInside = fcn_MapGen_isWithinABBB(AABB,test_points,fig_num);

assert(true);


