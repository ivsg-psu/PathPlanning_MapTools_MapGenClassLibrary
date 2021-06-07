% script_test_fcn_MapGen_polytopeMapGen
% Tests function: fcn_MapGen_polytopeMapGen

% REVISION HISTORY:
% 2021_06_06 
% -- first written by S. Brennan. 


% generate Voronoi tiling from Halton points
Halton_range = [1 200]; % range of Halton points to use to generate the tiling
% remove the edge polytope that extend past the high and low points
xlow = 0; xhigh = 1; ylow = 0; yhigh = 1;
bounding_box = [xlow ylow; xhigh yhigh];

% shink the polytopes so that they are no longer tiled
des_radius = 0.03; % desired average maximum radius
sigma_radius = 0.002; % desired standard deviation in maximum radii
min_rad = 0.0001; % minimum possible maximum radius for any obstacle
shrink_seed = 1111; % seed used for randomizing the shrinking process
fig_num = 1;

[map_polytopes,all_pts] = ...
    fcn_MapGen_polytopeMapGen(...
    Halton_range,bounding_box,...
    des_radius,sigma_radius,min_rad,shrink_seed,fig_num);


% Initiate the plot
fig = 103; % figure to plot on
line_spec = 'b-'; % edge line plotting
line_width = 2; % linewidth of the edge
axes_limits = [0 1 0 1]; % x and y axes limits
axis_style = 'square'; % plot axes style
fcn_MapGen_plotPolytopes(map_polytopes,fig,line_spec,line_width,axes_limits,axis_style);
