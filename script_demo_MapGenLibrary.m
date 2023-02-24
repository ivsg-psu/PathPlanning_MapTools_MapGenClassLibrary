% script_demo_MapGenLibrary.m
% This is a script that shows the capabilities of the "MapGen" class of
% functions via demonstrations.

% Revision history:
% 2021_06_07:
% -- First write of the function, using the "Vis" library demo script as
% starter
% 2021_06_09
% -- Added other types of point generators
% 2021_07_06
% -- Updated to include the newer expansion functions
% 2021_07_11
% -- Add ability to extend halton set to right (e.g. "scrolling" map), see
% the function: fcn_MapGen_mixedSetVoronoiTiling
% 2021_07_12
% -- Added ability to determine generic map statistics via the function:
% fcn_MapGen_polytopesStatistics
% 2023_01_15
% -- Added demo of edge-based shrinking
% 2023_02_20
% -- Added code to better support README.md
% 2023_02_21
% -- Added Debug utility library


% TO-DO:
% -- add debug library utility, and switch functions to this
% -- add functions that, given a map, give core statistics (look out limit, linear density, etc - basically make functions to calculate all the pi-values and interpretations we might need)
% -- add prior work on grid-based map generation

% Say what the libraries we need
library_name{1}    = 'DebugTools_v2023_01_29';
library_folders{1} = {'Functions','Data'};
library_url{1}     = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/blob/main/Releases/DebugTools_v2023_01_29.zip?raw=true';

% library_name{2}    = 'CodeX_Functions';
% library_folders{2} = {};
% library_url{2}     = 'https://github.com/ivsg-psu/Errata_Tutorials_CodeExercises/blob/main/Releases/CodeX_2023_02_03.zip?raw=true';

close all

%% Do we need to set up the work space?
if ~exist('flag_was_run_before','var') || isempty(flag_was_run_before) 
    fprintf(1,'Installing folders, getting codes ready ...\n');

    % Dependencies and Setup of the Code
    % This code depends on several other libraries of codes that contain
    % commonly used functions. We check to see if these libraries are installed
    % into our "Utilities" folder, and if not, we install them and then set a
    % flag to not install them again.
    
    % Set up dependent libraries?
    if exist('library_name','var')
        for ith_library = 1:length(library_name)
            dependency_name = library_name{ith_library};
            dependency_subfolders = library_folders{ith_library};
            dependency_url = library_url{ith_library};

            fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url);
            clear dependency_name dependency_subfolders dependency_url
        end
    end

    % Set dependencies for this project? 
    if ~exist('flag_functionsAdded','var') || isempty(flag_functionsAdded)
        fcn_DebugTools_addSubdirectoriesToPath(pwd,{'Functions'});
        flag_functionsAdded = 1;
    end
    
    flag_was_run_before = 1;
    disp('Done setting up folders.');
    
end


% %% Set up workspace
% clear flag_was_run_before  % Force init to always run?
% 
% if ~exist('flag_was_run_before','var')
%     
%     clc
%     close all
%     
%     % add necessary directories
%     addpath([pwd '\Functions'])
%     %     addpath([pwd '\GeomClassLibrary\Functions'])
%     %     addpath([pwd '\MapGenClassLibrary\Functions'])
%     %     addpath([pwd '\Plotting'])
%     %     addpath([pwd '\Map_Generation\polytope_generation'])
%     %     addpath([pwd '\Map_Generation\polytope_editing'])
%     %     addpath([pwd '\Map_Generation\polytope_calculation'])
%     
%     flag_was_run_before = 1;
% end

%% Show how inputs are checked
Twocolumn_of_numbers_test = [4 1; 3 0; 2 5];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%    _____ _             _        _____      _       _                    ____                       _   _                 
%   / ____(_)           | |      |  __ \    | |     | |                  / __ \                     | | (_)                
%  | (___  _ _ __   __ _| | ___  | |__) |__ | |_   _| |_ ___  _ __   ___| |  | |_ __   ___ _ __ __ _| |_ _  ___  _ __  ___ 
%   \___ \| | '_ \ / _` | |/ _ \ |  ___/ _ \| | | | | __/ _ \| '_ \ / _ \ |  | | '_ \ / _ \ '__/ _` | __| |/ _ \| '_ \/ __|
%   ____) | | | | | (_| | |  __/ | |  | (_) | | |_| | || (_) | |_) |  __/ |__| | |_) |  __/ | | (_| | |_| | (_) | | | \__ \
%  |_____/|_|_| |_|\__, |_|\___| |_|   \___/|_|\__, |\__\___/| .__/ \___|\____/| .__/ \___|_|  \__,_|\__|_|\___/|_| |_|___/
%                   __/ |                       __/ |        | |               | |                                         
%                  |___/                       |___/         |_|               |_|                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Show how to check if points are within an Axis-Aligned Bounding Box
AABB = [0 0 1 1]; % Define the axis-aligned bounding box
test_points = randn(100,2);
fig_num = 1;
isInside = fcn_MapGen_isWithinABBB(AABB,test_points,fig_num);

%% Show how to plot a polytope

%% Show how we calculate the polytope centroid and area
% Note: this does NOT have to be a convex polytope, as the example shows.
% Copied from script_test_fcn_MapGen_polytopeCentroidAndArea
% Tests: fcn_MapGen_polytopeCentroidAndArea

fig_num = 2;

x = [3; 4; 2; -1; -2; -3; -4; -2; 1; 2; 3];
y = [1; 2; 2; 3; 2; -1; -2; -3; -3; -2; 1];
[Centroid,Area] = fcn_MapGen_polytopeCentroidAndArea([x,y],fig_num);

assert(isequal(round(Centroid,4),[-0.1462,-0.2222]));
assert(isequal(round(Area,4),28.5));

figure(2);
plot(x,y,'g-','linewidth',2)
hold on
plot(Centroid(:,1),Centroid(:,2),'kx','linewidth',1);


%% Show how, if we have only the vertices of a polytope, we can calculate all other fields 
% Also shows how to plot the polytopes
% Copied from script_test_fcn_MapGen_fillPolytopeFieldsFromVerticies
% Tests fcn_MapGen_fillPolytopeFieldsFromVerticies
% Given a polytoope structure array where the vertices field is filled, 
% calculates the values for all the other fields.


fig_num = 3;
clear polytopes
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes);
line_width = 3;
fcn_MapGen_plotPolytopes(polytopes,fig_num,'r-',line_width);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   _____      _       _                   ______ _      _     _  ____                       _   _                 
%  |  __ \    | |     | |                 |  ____(_)    | |   | |/ __ \                     | | (_)                
%  | |__) |__ | |_   _| |_ ___  _ __   ___| |__   _  ___| | __| | |  | |_ __   ___ _ __ __ _| |_ _  ___  _ __  ___ 
%  |  ___/ _ \| | | | | __/ _ \| '_ \ / _ \  __| | |/ _ \ |/ _` | |  | | '_ \ / _ \ '__/ _` | __| |/ _ \| '_ \/ __|
%  | |  | (_) | | |_| | || (_) | |_) |  __/ |    | |  __/ | (_| | |__| | |_) |  __/ | | (_| | |_| | (_) | | | \__ \
%  |_|   \___/|_|\__, |\__\___/| .__/ \___|_|    |_|\___|_|\__,_|\____/| .__/ \___|_|  \__,_|\__|_|\___/|_| |_|___/
%                 __/ |        | |                                     | |                                         
%                |___/         |_|                                     |_|                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Show a detailed step-by-step process behind construction of obstacle map
fig_num = 1010;
AABB = [0 0 1 1]; % Define the axis-aligned bounding box
scale = max(AABB,[],'all') - min(AABB,[],'all');
new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];



% Fill in the Halton set
% pull halton set
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values


% pick values from halton set
Halton_range = [1 100];
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
seed_points = points_scrambled(low_pt:high_pt,:);
[V,C] = voronoin(seed_points);
stretch = [1 1];

% fill polytopes from tiling
[polytopes,all_vertices] = fcn_MapGen_generatePolysFromTiling(seed_points,V,C,AABB, stretch);


% PLOT THE SEED POINTS
figure(fig_num);
clf;

subplot(2,3,1);
% plot the seed points in red
plot(seed_points(:,1),seed_points(:,2),'r.','Markersize',10);


% number the polytopes at seed points
for ith_poly = 1:length(polytopes)
    text_location = seed_points(ith_poly,:);
    text(text_location(1,1),text_location(1,2),sprintf('%.0d',ith_poly));
end
axis(new_axis);
title('Seed points');

% PLOT THE VORONOI lines with the points
subplot(2,3,2);

% plot the seed points in red
plot(seed_points(:,1),seed_points(:,2),'r.','Markersize',10);
hold on;

% plot all vertices
plot(all_vertices(:,2),all_vertices(:,3),'c','Linewidth',1);

axis(new_axis);
title('Voronoi boundaries');


% PLOT THE VORONOI lines with the points
subplot(2,3,3);

% plot the polytopes on current axis
fcn_MapGen_plotPolytopes(polytopes,gca,'b',2);
hold on;


% plot the seed points in red
plot(seed_points(:,1),seed_points(:,2),'r.','Markersize',10);

% plot all vertices
plot(all_vertices(:,2),all_vertices(:,3),'c','Linewidth',1);


axis(new_axis);
title('AABB imposed')


% plot the means in black
subplot(2,3,4);


% plot the polytopes
fcn_MapGen_plotPolytopes(polytopes,gca,'b',2);
hold on;

temp = zeros(length(polytopes),2);
for ith_poly = 1:length(polytopes)
    temp(ith_poly,:) = polytopes(ith_poly).mean;
end
plot(temp(:,1),temp(:,2),'ko','Markersize',3);


% plot the means in black
temp = zeros(length(polytopes),2);
for ith_poly = 1:length(polytopes)
    temp(ith_poly,:) = polytopes(ith_poly).mean;
end
plot(temp(:,1),temp(:,2),'ko','Markersize',3);

axis(new_axis);
title('Polytope stats (mean)');

% plot the shrink to radius
subplot(2,3,5);

des_rad = 0.05; sigma_radius = 0; min_rad = 0.001;
shrunk_polytopes=fcn_MapGen_polytopesShrinkToRadius(...
    polytopes,des_rad,sigma_radius,min_rad);

% plot the shrunk polytopes
fcn_MapGen_plotPolytopes(shrunk_polytopes,gca,'r',2);
hold on;

axis(new_axis);
title('Shrunk to radius polytopes');

% plot the shrink to edge
subplot(2,3,6);

des_gap_size = 0.05;

shrunk_polytopes2=...
    fcn_MapGen_polytopesShrinkFromEdges(...
    polytopes,des_gap_size);

 % plot the shrunk polytopes
fcn_MapGen_plotPolytopes(shrunk_polytopes2,gca,'r',2);
hold on;

axis(new_axis);
title('Shrunk from edge polytopes');



%% Show typical generation of polytopes from a tiling
% Demos fcn_MapGen_generatePolysFromTiling
% and script_test_fcn_MapGen_generatePolysFromTiling
fig_num = 10;

% pull halton set
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values


% pick values from halton set
Halton_range = [1 100];
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
seed_points = points_scrambled(low_pt:high_pt,:);
[V,C] = voronoin(seed_points);

AABB = [0 0 1 1]; % Define the axis-aligned bounding box
stretch = [1 1];

% fill polytopes from tiling
fcn_MapGen_generatePolysFromTiling(seed_points,V,C,AABB, stretch,fig_num);

%% Show time-dependent generation of polytopes from a tiling
% Moves just one point
fig_num = 11111;
clf;


% pull halton set
rng(1111);
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values
AABB = [0 0 1 1]; % Define the axis-aligned bounding box
stretch = [1 1];


% pick values from halton set
Npoints = 100;
Halton_range = [1 Npoints];
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
original_seed_points = points_scrambled(low_pt:high_pt,:);

% Set the point to move
pointID_to_move = 9;

% Set the velocity vectors
%velocities = 0.01*randn(Npoints,2);
velocities = 0.01;

% Prep the movie
% vidfile = VideoWriter('MovingPolytope.mp4','MPEG-4');
vidfile = VideoWriter('MovingPolytope.avi');
open(vidfile);

moved_seed_points = original_seed_points;
for ith_step = 1:1:100
    moved_seed_points(pointID_to_move) = original_seed_points(pointID_to_move) + velocities*ith_step;

    % Wrap points
    seed_points = mod(moved_seed_points,1);

    [V,C] = voronoin(seed_points);


    % fill polytopes from tiling
    polytopes = fcn_MapGen_generatePolysFromTiling(seed_points,V,C,AABB, stretch);

    des_gap_size = 0.02;

    shrunk_polytopes=...
        fcn_MapGen_polytopesShrinkFromEdges(...
        polytopes,des_gap_size);


    % Plot the results
    figure(fig_num);
    clf;

    line_width = 2;
    color = [0 0 0];
    axis_box = [0 1 0 1];
    fcn_MapGen_plotPolytopes(shrunk_polytopes,gca,'-',line_width,color,axis_box,'square'); %,[1 0 0 0 0.5]);
    hold on;
    highlighted_polytope = shrunk_polytopes(pointID_to_move);

    % Highlight the moving polytope
    % FILL_INFO: a 1-by-5 vector to specify wether or not there is fill, the
    % color of fill, and the opacity of the fill [Y/N, R, G, B, alpha]
    % fill_info = [1 1 0 0 1];
    fcn_MapGen_plotPolytopes(highlighted_polytope,gca,'-',line_width,[1 0 0],axis_box,'square');  % , fill_info);


    frame = getframe(gcf);
    writeVideo(vidfile,frame);


    pause(0.1);
end
close(vidfile)


%% Show time-dependent generation of polytopes from a tiling
% Moves all points
fig_num = 1111;
clf;


% pull halton set
rng(1111);
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values
AABB = [0 0 1 1]; % Define the axis-aligned bounding box
stretch = [1 1];


% pick values from halton set
Npoints = 100;
Halton_range = [1 Npoints];
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
original_seed_points = points_scrambled(low_pt:high_pt,:);

% Set the point to move
% pointID_to_move = 9;

% Set the velocity vectors
velocities = 0.01*randn(Npoints,2);
% velocities = 0.01;

% Prep the movie
% vidfile = VideoWriter('MovingPolytope.mp4','MPEG-4');
vidfile = VideoWriter('MovingPolytope.avi');
open(vidfile);

% Set figure parameters
line_width = 2;
color = [0 0 0];
axis_box = [0 1 0 1];


moved_seed_points = original_seed_points;
for ith_step = 1:1:100
    % moved_seed_points(pointID_to_move) = original_seed_points(pointID_to_move) + velocities*ith_step;
    moved_seed_points = original_seed_points + velocities*ith_step;

    % Wrap points to all exist between 0 and 1
    seed_points = mod(moved_seed_points,1);

    % Calculate the Voronoi diagram
    [V,C] = voronoin(seed_points);

    % fill polytopes from tiling of Voronoi diagram
    polytopes = fcn_MapGen_generatePolysFromTiling(seed_points,V,C,AABB, stretch);

    des_gap_size = 0.02;

    shrunk_polytopes=...
        fcn_MapGen_polytopesShrinkFromEdges(...
        polytopes,des_gap_size);


    % Plot the results
    figure(fig_num);
    clf;

    fcn_MapGen_plotPolytopes(shrunk_polytopes,gca,'-',line_width,color,axis_box,'square'); %,[1 0 0 0 0.5]);
    
    %     % Highlight the moving polytope
    %     % FILL_INFO: a 1-by-5 vector to specify wether or not there is fill, the
    %     % color of fill, and the opacity of the fill [Y/N, R, G, B, alpha]
    %     % fill_info = [1 1 0 0 1];
    %     fcn_MapGen_plotPolytopes(highlighted_polytope,gca,'-',line_width,[1 0 0],axis_box,'square');  % , fill_info);


    frame = getframe(gcf);
    writeVideo(vidfile,frame);


    pause(0.1);
end
close(vidfile)


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

%% Show how to create an overlapping set using different AABBs for each set.
close all;
clear mixedSet

fig_num = 1;
stretch = [1 1];
set_range = [1 100];

rng(1234);

mixedSet(1).name = 'haltonset';
mixedSet(1).settings = set_range;
mixedSet(1).AABB = [0 0 1 1];

mixedSet(2).name = 'rand';
mixedSet(2).settings = set_range;
mixedSet(2).AABB = [0.5 0 0.75 1];

polytopes = fcn_MapGen_mixedSetVoronoiTiling(mixedSet,stretch,fig_num);


%% Generate many test sets of polytopes from the Halton set
for i=1:100:500
    fig_num = 21+i;
    figure(fig_num); clf;
    
    Halton_range = [i i+100]; % range of Halton points to use to generate the tiling
    % Halton_range = [1801 1901];
           
    polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);
    
    % Do statistics, checking that the area is always fully filled and we
    % get 101 polytopes each time
    temp = fcn_MapGen_polytopesStatistics(...
    polytopes);
    title(sprintf('Halton range is: [%.0d %.0d]',i,i+100));
    assert(abs(temp.unoccupancy_ratio)<(1000*eps));
    assert(isequal(101,temp.point_density));
    pause(0.1);
end

%% Show how the maps can be trimmed to a box

% Generate polytopes from the Halton set
Halton_range = [5401 5501];

fig_num = 31;
tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);
    
fig_num = 32;
fcn_MapGen_polytopesStatistics(...
    tiled_polytopes,...
    fig_num);

% Plot the polytopes
fig_num = 33;
line_width = 2;
axis_limits = [0 1 0 1];
fcn_MapGen_plotPolytopes(tiled_polytopes,fig_num,'r',line_width,axis_limits);

% remove the edge polytopes that extend on or past the high and low points
fig_num = 34;
xlow = 0.01; xhigh = 0.99; ylow = 0.01; yhigh = 0.99;
bounding_box = [xlow ylow; xhigh yhigh];
trimmed_polytopes = ...
    fcn_MapGen_polytopeCropEdges(tiled_polytopes,bounding_box,fig_num);


%% Show how the polytopes can be shrunk to a specified radius
% Shrink to radius
fig_num = 24;
des_rad = 0.03; sigma_radius = 0; min_rad = 0.001;
shrunk_polytopes2=fcn_MapGen_polytopesShrinkToRadius(...
    trimmed_polytopes,des_rad,sigma_radius,min_rad,fig_num);

%% Show how different shrinking methods change statistics
fig_num = 555;

% Generate polytopes from the Halton set
Halton_range = [5401 5501];
           
tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);
    
% Grab statistics on original map
fcn_MapGen_polytopesStatistics(...
    tiled_polytopes,...
    fig_num+1);

% Shrink to radius
URHERE
fig_num = 556;
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

%% Show calculation of map statistics
fig_num = 9;
fcn_MapGen_polytopesStatistics(...
    polytopes,...
    fig_num);

%% Remove close points and redo statistics
tolerance = 0.01;
cleaned_polytopes = polytopes;
for ith_poly = 1:length(polytopes)
    cleaned_polytopes(ith_poly) = fcn_MapGen_polytopeRemoveTightVerticies(...
        polytopes(ith_poly), tolerance);    
end

fig_num = 20;
fcn_MapGen_polytopesStatistics(...
    cleaned_polytopes,...
    fig_num);


%% Generating starting map for UGV Error Bubbles and Plotting functions (Nick's work)

% create polytopes
stretch = [200, 200]; % stretch in the x and y directions
Halton_range = [1 1000]; % range of Halton points to use to generate the tiling
polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,stretch);

% Plot the polytopes
fig_num = 22;
line_width = 2;
axis_limits = [0 200 0 200];
axis_stype = 'square';
fcn_MapGen_plotPolytopes(polytopes,fig_num,'b',line_width,axis_limits,axis_stype);


% remove the edge polytopes that extend past the high and low points
fig_num = 23;
bounding_box = [0 0; 200 200];
trimmed_polytopes = ...
    fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);

%shrink polytopes to create space
fig_num = 24;
des_rad = 1; sigma_radius = 0.5; min_rad = 0.25;
shrunk_polytopes2=fcn_MapGen_polytopesShrinkToRadius(...
    trimmed_polytopes,des_rad,sigma_radius,min_rad,fig_num);




% generate error bubbles via fcn_MapGen_ugvSensorErrorBubble
[err] = fcn_MapGen_ugvSensorErrorBubble(shrunk_polytopes2, 0, 5);

% Convert error bounds into polytope structure
error_polytopes = shrunk_polytopes2; % Initialize the structure
for ii=1:length(shrunk_polytopes2)
    error_polytopes(ii).vertices = [err(ii).circ_x(err(ii).bubble)', err(ii).circ_y(err(ii).bubble)'];
end
error_polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(error_polytopes);

%verify
h_fig = figure('name','UGV Positioning Bubbles');
fig_num = h_fig.Number;
ax.ugv1=gca;
% hold on
% for ii=1:length(shrunk_polytopes2)
%     plot(shrunk_polytopes2(ii).vertices(:,1),shrunk_polytopes2(ii).vertices(:,2),'k')
%     plot(ax.ugv1,err(ii).circ_x(err(ii).bubble),err(ii).circ_y(err(ii).bubble),'r')
% end
% hold off

fcn_MapGen_plotPolytopes(shrunk_polytopes2,fig_num,'k',line_width,axis_limits,axis_stype);
fcn_MapGen_plotPolytopes(error_polytopes,fig_num,'r',line_width,axis_limits,axis_stype);
xlabel('X Distance [m]')
ylabel('Y Distance [m]')
title('UGV Positioning Bubbles')
legend('Real Object','Perceived Object')
grid on
axis equal

%% Testing work with contour plots, using function fcn_MapGen_ugvSensorError

%generating R and beta values for a grid of points
x = linspace(0,200,200);
y = linspace(100,-100,200);
[X,Y] = meshgrid(x,y);

R = sqrt(X.^2+Y.^2);  % perceived distance, nearly constant for flat object
beta = rad2deg(atan(Y./X));
kappa = zeros(size(R));

clear err;
[err.x, err.y, err.z] = fcn_MapGen_ugvSensorError({R, beta, kappa}, ...
    {0.08, 0.08, 0.08}, {0.03, 0.03, 0.4}, {0.02, -0.05});

err.R = sqrt(err.x.^2 + err.y.^2);

figure('name','Error in X')
contour(x,y,err.x,'ShowText','on')
h=colorbar;
title('Error in X Dimension [m]')
xlabel('X [m]')
ylabel('Y [m]')

figure('name','Error in Y')
contour(x,y,err.y,'ShowText','on')
h=colorbar;
title('Error in Y Dimension [m]')
xlabel('X [m]')
ylabel('Y [m]')

figure('name','Total Error (Bubble Radius) [m]')
contour(x,y,err.R,'ShowText','on')
h=colorbar;
title('Total Error (Bubble Radius) [m]')
xlabel('X [m]')
ylabel('Y [m]')



%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§


function fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url, varargin)
%% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES - MATLAB package installer from URL
%
% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES installs code packages that are
% specified by a URL pointing to a zip file into a default local subfolder,
% "Utilities", under the root folder. It also adds either the package
% subfoder or any specified sub-subfolders to the MATLAB path.
%
% If the Utilities folder does not exist, it is created.
% 
% If the specified code package folder and all subfolders already exist,
% the package is not installed. Otherwise, the folders are created as
% needed, and the package is installed.
% 
% If one does not wish to put these codes in different directories, the
% function can be easily modified with strings specifying the
% desired install location.
% 
% For path creation, if the "DebugTools" package is being installed, the
% code installs the package, then shifts temporarily into the package to
% complete the path definitions for MATLAB. If the DebugTools is not
% already installed, an error is thrown as these tools are needed for the
% path creation.
% 
% Finally, the code sets a global flag to indicate that the folders are
% initialized so that, in this session, if the code is called again the
% folders will not be installed. This global flag can be overwritten by an
% optional flag input.
%
% FORMAT:
%
%      fcn_DebugTools_installDependencies(...
%           dependency_name, ...
%           dependency_subfolders, ...
%           dependency_url)
%
% INPUTS:
%
%      dependency_name: the name given to the subfolder in the Utilities
%      directory for the package install
%
%      dependency_subfolders: in addition to the package subfoder, a list
%      of any specified sub-subfolders to the MATLAB path. Leave blank to
%      add only the package subfolder to the path. See the example below.
%
%      dependency_url: the URL pointing to the code package.
%
%      (OPTIONAL INPUTS)
%      flag_force_creation: if any value other than zero, forces the
%      install to occur even if the global flag is set.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      This code will automatically get dependent files from the internet,
%      but of course this requires an internet connection. If the
%      DebugTools are being installed, it does not require any other
%      functions. But for other packages, it uses the following from the
%      DebugTools library: fcn_DebugTools_addSubdirectoriesToPath
%
% EXAMPLES:
%
% % Define the name of subfolder to be created in "Utilities" subfolder
% dependency_name = 'DebugTools_v2023_01_18';
%
% % Define sub-subfolders that are in the code package that also need to be
% % added to the MATLAB path after install; the package install subfolder
% % is NOT added to path. OR: Leave empty ({}) to only add 
% % the subfolder path without any sub-subfolder path additions. 
% dependency_subfolders = {'Functions','Data'};
%
% % Define a universal resource locator (URL) pointing to the zip file to
% % install. For example, here is the zip file location to the Debugtools
% % package on GitHub:
% dependency_url = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/blob/main/Releases/DebugTools_v2023_01_18.zip?raw=true';
%
% % Call the function to do the install
% fcn_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url)
%
% This function was written on 2023_01_23 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_01_23:
% -- wrote the code originally

% TO DO
% -- Add input argument checking

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % % Flag to plot the final results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end


%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_check_inputs
    % Are there the right number of inputs?
    narginchk(3,4);
end

%% Set the global variable - need this for input checking
% Create a variable name for our flag. Stylistically, global variables are
% usually all caps.
flag_varname = upper(cat(2,'flag_',dependency_name,'_Folders_Initialized'));

% Make the variable global
eval(sprintf('global %s',flag_varname));

if nargin==4
    if varargin{1}
        eval(sprintf('clear global %s',flag_varname));
    end
end

%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~exist(flag_varname,'var') || isempty(eval(flag_varname))
    fprintf(1,'Attempting install of the sub-utility: %s\n',dependency_name);
    % Save the root directory, so we can get back to it after some of the
    % operations below. We use the Print Working Directory command (pwd) to
    % do this. Note: this command is from Unix/Linux world, but is so
    % useful that MATLAB made their own!
    root_directory_name = pwd;

    % Does the directory "Utilities" exist?
    utilities_folder_name = fullfile(root_directory_name,'Utilities');
    if ~exist(utilities_folder_name,'dir')
        fprintf(1,'Main Utility subfolder does not exist - attempting creation of this folder.\n');
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(root_directory_name,'Utilities');

        % Did it work?
        if ~success_flag
            error('Unable to make the Utilities directory. Reason: %s with message ID: %s\n',error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The Utilities directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',error_message, message_ID);
        else
            fprintf(1,'Main Utility subfolder was successfully created.\n');
        end
    else
        fprintf(1,'Main Utility subfolder already exists - skipping creation of this folder.\n');
    end

    % Does the directory for the dependency folder exist?
    dependency_folder_name = fullfile(root_directory_name,'Utilities',dependency_name);
    fprintf(1,'Attempting creation of the sub-utility folder for: %s\n',dependency_name);
    if ~exist(dependency_folder_name,'dir')

        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(utilities_folder_name,dependency_name);

        % Did it work?
        if ~success_flag
            error('Unable to make the dependency directory: %s. Reason: %s with message ID: %s\n',dependency_name, error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The %s directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',dependency_name, error_message, message_ID);
        else
            fprintf(1,'Sub-utility folder was successfully created.\n');
        end
    else
        fprintf(1,'Sub-utility subfolder already exists - skipping creation of this folder.\n');
    end

    % Do the subfolders exist?
    flag_allFoldersThere = 1;
    if isempty(dependency_subfolders)
        flag_allFoldersThere = 0;
        fprintf(1,'It appears the sub-utility has no folder sub-dependencies. None are checked: %s\n',dependency_name);
    else
        fprintf(1,'It appears the sub-utility has sub-folders. Checking if subfolders are already there for: %s\n',dependency_name);

        for ith_folder = 1:length(dependency_subfolders)
            subfolder_name = dependency_subfolders{ith_folder};
            fprintf(1,'\t Checking: %s\n',subfolder_name);
            
            % Create the entire path
            subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
            
            % Check if the folder and file exists that is typically created when
            % unzipping.
            if ~exist(subfunction_folder,'dir')
                flag_allFoldersThere = 0;
            end
        end
    end

    % Do we need to unzip the files?
    fprintf(1,'Checking if download process needed for for: %s\n',dependency_name);
    if flag_allFoldersThere==0
        fprintf(1,'Some dowload is needed for: %s\n',dependency_name);
        % Files do not exist yet - try unzipping them.
        save_file_name = tempname(root_directory_name);
        zip_file_name = websave(save_file_name,dependency_url);
  
        % Is the file there?
        if ~exist(zip_file_name,'file')
            error('The zip file: %s for dependency: %s did not download correctly. This is usually because permissions are restricted on the current directory. Check the code install (see README.md) and try again.\n',zip_file_name, dependency_name);
        else
            fprintf(1,'Download files found and copied locally for: %s\n',dependency_name);
        end

        % Try unzipping
        unzip(zip_file_name, dependency_folder_name);

        % Did this work?
        fprintf(1,'Checking that unzip worked for: %s\n',dependency_name);
        flag_allFoldersThere = 1;
        if ~isempty(dependency_subfolders)
            for ith_folder = 1:length(dependency_subfolders)
                subfolder_name = dependency_subfolders{ith_folder};
                
                % Create the entire path
                subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
                
                % Check if the folder and file exists that is typically created when
                % unzipping.
                if ~exist(subfunction_folder,'dir')
                    flag_allFoldersThere = 0;
                end
            end
        end
        
        if flag_allFoldersThere==0
            error('The necessary dependency: %s has an error in install, or error performing an unzip operation. Check the code install (see README.md) and try again.\n',dependency_name);
        else
            % Clean up the zip file
            delete(zip_file_name);
            fprintf(1,'Installation complete: %s\n',dependency_name);
        end
    else
        fprintf(1,'No downloads seem to be needed for: %s\n',dependency_name);
    end


    % For path creation, if the "DebugTools" package is being installed, the
    % code installs the package, then shifts temporarily into the package to
    % complete the path definitions for MATLAB. If the DebugTools is not
    % already installed, an error is thrown as these tools are needed for the
    % path creation.
    %
    % In other words: DebugTools is a special case because folders not
    % added yet, and we use DebugTools for adding the other directories
    if strcmp(dependency_name(1:10),'DebugTools')
        fprintf(1,'Adding path for Debug tools via file commands. Version installed is: %s\n',dependency_name);

        debugTools_function_folder = fullfile(root_directory_name, 'Utilities', dependency_name,'Functions');

        % Move into the folder, run the function, and move back
        cd(debugTools_function_folder);
        fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        cd(root_directory_name);
    else
        try
            fprintf(1,'Adding path for package to MATLAB environment: %s\n',dependency_folder_name);
            fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        catch
            error('Package installer requires DebugTools package to be installed first. Please install that before installing this package');
        end
        fprintf(1,'Path successfully added for: %s\n',dependency_folder_name);
    end


    % Finally, the code sets a global flag to indicate that the folders are
    % initialized.  Check this using a command "exist", which takes a
    % character string (the name inside the '' marks, and a type string -
    % in this case 'var') and checks if a variable ('var') exists in matlab
    % that has the same name as the string. The ~ in front of exist says to
    % do the opposite. So the following command basically means: if the
    % variable named 'flag_CodeX_Folders_Initialized' does NOT exist in the
    % workspace, run the code in the if statement. If we look at the bottom
    % of the if statement, we fill in that variable. That way, the next
    % time the code is run - assuming the if statement ran to the end -
    % this section of code will NOT be run twice.
    fprintf(1,'Setting global flag so that install does not repeat for: %s\n\n',dependency_name);
    eval(sprintf('%s = 1;',flag_varname));
end

%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plots

    % Nothing to do!



end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function fcn_DebugTools_installDependencies




