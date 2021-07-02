%% main code §
x = [3; 4; 2; -1; -2; -3; -4; -2; 1; 2; 3];
y = [1; 2; 2; 3; 2; -1; -2; -3; -3; -2; 1];
[Centroid,Area] = fcn_polytope_calculation_centroid_and_area([x,y]);
plot(x,y,'g-','linewidth',2)
hold on
plot(Centroid(:,1),Centroid(:,2),'kx','linewidth',1)