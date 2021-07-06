% script_demo_MapGenLibrary.m
% This is a script that shows the capabilities of the "MapGen" class of
% functions via demonstrations.

% Revision history:
% 2021_06_07:
% -- First write of the function, using the "Vis" library demo script as
% starter
% 2021_06_09
% -- Added other types of point generators

% TO-DO:
% -- finish the growth of polytopes functions, e.g. get these working with the library
% -- add Nick's sensor code 
% -- Add ability to extend halton set to right (e.g. "scrolling" map)
% -- add functions that, given a map, give core statistics (look out limit, linear density, etc - basically make functions to calculate all the pi-values and interpretations we might need)
% -- add prior work on grid-based map generation


close all

%% Set up workspace
clear flag_was_run_before  % Force init to always run?

if ~exist('flag_was_run_before','var')
    
    clc
    close all
    
    % add necessary directories
    addpath([pwd '\Functions'])
    %     addpath([pwd '\GeomClassLibrary\Functions'])
    %     addpath([pwd '\MapGenClassLibrary\Functions'])
    %     addpath([pwd '\Plotting'])
    %     addpath([pwd '\Map_Generation\polytope_generation'])
    %     addpath([pwd '\Map_Generation\polytope_editing'])
    %     addpath([pwd '\Map_Generation\polytope_calculation'])
    
    flag_was_run_before = 1;
end

%% Show how inputs are checked
Twocolumn_of_numbers_test = [4 1; 3 0; 2 5];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers');

%% Show how to generate polytopes from a tiling
fig_num = 10;

% pull halton set
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values

% pick values from halton set
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
seed_points = points_scrambled(low_pt:high_pt,:);
[V,C] = voronoin(seed_points);
% V = V.*stretch;

% fill polytopes from tiling
polytopes = fcn_MapGen_generatePolysFromTiling(seed_points,V,C,stretch,fig_num);

%% Generate a set of polytopes from various pseudo-random sources
close all;

% Generate a set of polytopes from the Sobol set
fig_num = 11;
Sobol_range = [1 1000]; % range of Sobol points to use to generate the tiling
tiled_polytopes = fcn_MapGen_sobolVoronoiTiling(Sobol_range,[1 1],fig_num);
title('Sobel set');


% Generate a set of polytopes from the Halton set
fig_num = 12;
Halton_range = [1 1000]; % range of Halton points to use to generate the tiling
tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);
title('Halton set');

% Generate a set of polytopes from the Latin Hypercube set
fig_num = 13;
Latin_range = [1 1000]; % range of Halton points to use to generate the tiling
tiled_polytopes = fcn_MapGen_latinVoronoiTiling(Latin_range,[1 1],fig_num);
title('Latin Hypercube set');


% Generate a set of polytopes from the Random set
fig_num = 14;
Rand_range = [1 1000]; % range of Halton points to use to generate the tiling
tiled_polytopes = fcn_MapGen_randVoronoiTiling(Rand_range,[1 1],fig_num);
title('Uniform random set');

% Generate a set of polytopes from the Random Normal set
fig_num = 15;
Rand_range = [1 1000]; % range of Halton points to use to generate the tiling
tiled_polytopes = fcn_MapGen_randomNormalVoronoiTiling(Rand_range,[1 1],fig_num);
title('Random normally distributed set');

%% Show how the maps can be trimmed, shrunk, etc
% Generate a set of polytopes from the Halton set
fig_num = 21;
Halton_range = [1 200]; % range of Halton points to use to generate the tiling
tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);

%% Plot the polytopes
fig_num = 22;
line_width = 2;
axis_limits = [0 1 0 1];
fcn_MapGen_plotPolytopes(tiled_polytopes,fig_num,'r',line_width,axis_limits);

%% remove the edge polytopes that extend past the high and low points
fig_num = 23;
xlow = 0.01; xhigh = 0.99; ylow = 0.01; yhigh = 0.99;
bounding_box = [xlow ylow; xhigh yhigh];
trimmed_polytopes = ...
    fcn_MapGen_polytopeCropEdges(tiled_polytopes,bounding_box,fig_num);

%% Shrink to radius
fig_num = 24;
des_rad = 0.03; sigma_radius = 0; min_rad = 0.001;
shrunk_polytopes2=fcn_MapGen_polytopesShrinkToRadius(...
    trimmed_polytopes,des_rad,sigma_radius,min_rad,fig_num);

%% Show how we can shrink one polytope
Npolys = length(trimmed_polytopes);
rand_poly = 1+floor(rand*Npolys);
shrinker = trimmed_polytopes(rand_poly);

fig_num = 31;
orig_radius = shrinker.max_radius;
ratios = (0.99:-0.05:0);

for ith_ratio = 1:length(ratios)
    des_rad = orig_radius*ratios(ith_ratio);
    tolerance = 1e-5; % This is the edge distance below which vertices are merged together in the polytope
    shrunk_polytope =...
        fcn_MapGen_polytopeShrinkToRadius(...
        shrinker,des_rad,tolerance,fig_num);
    pause(0.01);
end

%% Show results of removing tight verticies
fig_num = 32;
figure(fig_num);
clf;
orig_radius = shrinker.max_radius;
ratios = (0.99:-0.05:0);

for ith_ratio = 1:length(ratios)
    des_rad = orig_radius*ratios(ith_ratio);
    tolerance = 0.02;
    shrunk_polytope =...
        fcn_MapGen_polytopeShrinkToRadius(...
        shrinker,des_rad,tolerance);
    cleaned_polytope = fcn_MapGen_polytopeRemoveTightVerticies(...
        shrunk_polytope, tolerance,fig_num);
    pause(0.01);
end

%% Generate plots (all above steps) in just one function call
des_radius = 0.03; % desired average maximum radius
sigma_radius = 0.002; % desired standard deviation in maximum radii
min_rad = 0.0001; % minimum possible maximum radius for any obstacle
shrink_seed = 1111; % seed used for randomizing the shrinking process
fig_num = 5;

[map_polytopes,all_pts] = ...
    fcn_MapGen_polytopeMapGen(...
    Halton_range,bounding_box,...
    des_radius,sigma_radius,min_rad,shrink_seed,fig_num);

%% Generate a map from a name
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


%% Show how to expand one polytope
fig_num = 8;

one_polytope = fcn_MapGen_generateOneRandomPolytope;
exp_polytopes=fcn_MapGen_polytopesExpandEvenly(one_polytope,exp_dist,fig_num);

%% Show how to expand many polytopes
fig_num = 7;
exp_dist = 0.01;
exp_polytopes=fcn_MapGen_polytopesExpandEvenly(polytopes,exp_dist,fig_num);

