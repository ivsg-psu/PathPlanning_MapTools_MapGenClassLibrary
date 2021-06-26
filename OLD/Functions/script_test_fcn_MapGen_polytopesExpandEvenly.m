% script_test_fcn_MapGen_polytopesExpandEvenly
% Tests: fcn_MapGen_polytopesExpandEvenly

% 
% REVISION HISTORY:
% 
% 2018_11_17, Adjusted example code on 2021_04_28 by Seth Tau, Rebased on 2021_06_26 by S. Brennan by Seth Tau
% -- first write of script
%%%%%%%%%%%%%%§



xv = [-2 -1 1 2 2 1 -1 -2];
yv = [-1 -2 -2 -1 1 2 2 1];
polytopes.vertices = [[xv xv(1)]' [yv yv(1)]'];
polytopes.xv = xv;
polytopes.yv = yv;

polytopes.distances = sum((polytopes(1).vertices(1:end-1,:)-polytopes(1).vertices(2:end,:)).^2,2).^0.5;
[Cx,Cy,polytopes.area] = fcn_polytope_calculation_centroid_and_area([xv xv(1)],[yv yv(1)]);
polytopes.mean = [Cx, Cy];
polytopes.max_radius = max(sum((polytopes.vertices(1:end-1,:)-ones(length(xv),1)*polytopes.mean).^2,2).^0.5);
delta = 0.01;
exp_dist = 1;

% Call the function
exp_polytopes=fcn_MapGen_polytopesExpandEvenly(polytopes,delta,exp_dist);


% Plot the results
fcn_MapGen_plotPolytopes(polytopes,99,'r-',2);
fcn_MapGen_plotPolytopes(exp_polytopes,99,'b-',2,[-4 4 -4 4],'square');
legend('Original','Expanded')
box on
xlabel('X Position')
ylabel('Y Position')

