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
% delta = 0.01; % Set the delta value (what is this used for?)
exp_dist = 0.04; % Set the expansion distance
fig_num = 221; % Set the figure number

% Call the function
exp_polytopes=fcn_MapGen_polytopesExpandEvenly(polytopes,exp_dist,fig_num);

assert(isequal(round(exp_polytopes.area,4),0.0150));
assert(isequal(round(exp_polytopes.max_radius,4),0.1445));

% Call the OLD function
% exp_polytopes=fcn_MapGen_polytopesExpandEvenly_OLD(polytopes,delta,exp_dist,fig_num);

%% Second test
% Generate a map from a name
map_name = "HST 30 450 SQT 0 1 0 1 SMV 0.02 0.005 1e-6 1234";
plot_flag = 1; disp_name = [1, 0.05 -0.05, 0.5 0.5 0.5, 10];
line_style = '-'; line_width = 2; color = [0 0 1];
axis_limits = [0 1 -0.1 1]; axis_style = 'square';
fill_info = [1 1 0 1 0.5];
fig_num = 7;

[polytopes,fig]=fcn_MapGen_nameToMap(...
    map_name,...
    plot_flag,...
    disp_name,...
    fig_num,...
    line_style,...
    line_width,....
    color,...
    axis_limits,...
    axis_style,...
    fill_info);

% Set expansion parameters
exp_dist = 0.01; % Set the expansion distance
fig_num = 222; % Set the figure number

% Call the function
exp_polytopes=fcn_MapGen_polytopesExpandEvenly(polytopes,exp_dist,fig_num);
