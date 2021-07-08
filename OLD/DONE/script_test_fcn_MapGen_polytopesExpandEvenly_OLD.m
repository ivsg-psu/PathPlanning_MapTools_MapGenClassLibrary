% script_test_fcn_MapGen_polytopesExpandEvenly
% Tests: fcn_MapGen_polytopesExpandEvenly

% 
% REVISION HISTORY:
%§ 
% 2018_11_17, Seth Tau 
% -- first write of script
% 2021_04_28, Seth Tau
% -- Adjusted example code , 
% 2021_06_26 S. Brennan
% -- Rebased code

% Prep the workspace
close all;
clear polytopes;
polytopes = fcn_MapGen_generateOneRandomPolytope;

% xv = [-2 -1 1 2 2 1 -1 -2];
% yv = [-1 -2 -2 -1 1 2 2 1];
% polytopes.vertices = [[xv xv(1)]' [yv yv(1)]'];
% polytopes.xv = xv;
% polytopes.yv = yv;
% 
% polytopes.distances = sum((polytopes(1).vertices(1:end-1,:)-polytopes(1).vertices(2:end,:)).^2,2).^0.5;
% [Cx,Cy,polytopes.area] = fcn_MapGen_polytopeCentroidAndArea([xv xv(1)],[yv yv(1)]);
% polytopes.mean = [Cx, Cy];
% polytopes.max_radius = max(sum((polytopes.vertices(1:end-1,:)-ones(length(xv),1)*polytopes.mean).^2,2).^0.5);



polytopes.vertices = [
    1.0000    0.5217
    1.0000    0.5242
    0.9300    0.6329
    0.8472    0.6479
    0.8921    0.5627
    1.0000    0.5217
];            
polytopes.xv = [1 1 0.9300 0.8472 0.8921];
polytopes.yv = [0.5217 0.5242 0.6329 0.6479 0.5627];
polytopes.distances = [
    0.0025
    0.1293
    0.0842
    0.0963
    0.1154];
polytopes.mean = [0.9204 0.5894];
polytopes.area = 0.0079;
polytopes.max_radius = 0.1045;

% Set parameters
delta = 0.01; % Set the delta value (what is this used for?)
exp_dist = 0.04; % Set the expansion distance
fig_num = 222; % Set the figure number

% Call the function
exp_polytopes=fcn_MapGen_polytopesExpandEvenly(polytopes,delta,exp_dist,fig_num);


