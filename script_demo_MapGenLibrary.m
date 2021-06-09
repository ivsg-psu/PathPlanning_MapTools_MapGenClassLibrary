% script_demo_MapGenLibrary.m
% This is a script that shows the capabilities of the "MapGen" class of
% functions via demonstrations.

% Revision history:
% 2021_06_07:
% -- First write of the function, using the "Vis" library demo script as
% starter

% TO-DO:
% -- Add positive input checking to fcn_MapGen_polytopeShrinkToRadius
% -- Add ability to extend halton set to right (e.g. "scrolling" map)
% -- Add other types of point generators



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


%% Generate a set of polytopes from the Sobol set
fig_num = 11;
Sobol_range = [1 1000]; % range of Sobol points to use to generate the tiling
tiled_polytopes = fcn_MapGen_sobolVoronoiTiling(Sobol_range,[1 1],fig_num);


%% Generate a set of polytopes from the Halton set
fig_num = 12;
Halton_range = [1 1000]; % range of Halton points to use to generate the tiling
tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);


%% Generate a set of polytopes from the Latin Hypercube set
fig_num = 13;
Latin_range = [1 1000]; % range of Halton points to use to generate the tiling
tiled_polytopes = fcn_MapGen_latinVoronoiTiling(Latin_range,[1 1],fig_num);


%% Generate a set of polytopes from the Random set
fig_num = 14;
Rand_range = [1 1000]; % range of Halton points to use to generate the tiling
tiled_polytopes = fcn_MapGen_randVoronoiTiling(Rand_range,[1 1],fig_num);


%% Plot the polytopes
fig_num = 2;
line_width = 2;
axis_limits = [0 1 0 1];
fcn_MapGen_plotPolytopes(tiled_polytopes,fig_num,'r',line_width,axis_limits);

%% remove the edge polytopes that extend past the high and low points
fig_num = 3;
xlow = 0; xhigh = 1; ylow = 0; yhigh = 1;
bounding_box = [xlow ylow; xhigh yhigh];
trimmed_polytopes = ...
    fcn_MapGen_polytopeCropEdges(tiled_polytopes,bounding_box,fig_num);

%% Shrink to radius
fig_num = 4;
des_rad = 0.05; sigma_radius = 0.01; min_rad = 0.001;
shrunk_polytopes2=fcn_MapGen_polytopeShrinkToRadius(...
    trimmed_polytopes,des_rad,sigma_radius,min_rad,fig_num);

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

[polytopes,fig]=fcn_MapGen_nameToMap(map_name,plot_flag,disp_name,fig_num,line_style,line_width,color,axis_limits,axis_style,fill_info);

