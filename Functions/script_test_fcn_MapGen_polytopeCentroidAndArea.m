% script_test_fcn_MapGen_polytopeCentroidAndArea
% Tests: fcn_MapGen_polytopeCentroidAndArea

% 
% REVISION HISTORY:
% 
% 2021_07_02 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§

fig_num = 7;

x = [3; 4; 2; -1; -2; -3; -4; -2; 1; 2; 3];
y = [1; 2; 2; 3; 2; -1; -2; -3; -3; -2; 1];
[Centroid,Area] = fcn_MapGen_polytopeCentroidAndArea([x,y],fig_num);

figure(2);
plot(x,y,'g-','linewidth',2)
hold on
plot(Centroid(:,1),Centroid(:,2),'kx','linewidth',1)

