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
% 2023_03_13 
% -- Merged changes allowing for repeated tiling polytope field
% 2023_04_27
% -- Updated the installer to latest version of Debug Tools (fixes bug with
% Git archives of zips).
% 2023_05_04 - sbrennan@psu.edu
% -- Cleared the path variable, in case variable of same name shadows
% command. This was causing errors in some codes.
% -- functionalized the clear process
% 2025_04_25 - Sean Brennan
% -- Updated DebugTools_v2024_12_18 dependency
% -- Added global flags for setting test conditions and plotting
% -- Deprecated fcn_MapGen_checkInputsToFunctions, converted to fcn_DebugTools_checkInputsToFunctions
% 2025_06_25 - Sean Brennan
% -- upgraded debug tools to version DebugTools_v2025_06_24
% 2025_07_03 - Sean Brennan
% -- added minor comments in demo script
% 2025_07_07 - Sean Brennan
% -- started updating headers and test scripts. Only have 
% fcn_MapGen_generatePolysFromTiling done so far.
% 2025_07_11 - Sean Brennan
% -- updated DebugTools library
% -- added PathClass library to use this function (better) rather than
%    % fcn_MapGen_findIntersectionOfSegments
% -- for all the tiling variants, deprecated the following
%    % fcn_MapGen_sobolVoronoiTiling
%    % fcn_MapGen_latinVoronoiTiling
%    % (etc)
%    % Merged these to use fcn_MapGen_mixedSetVoronoiTiling and renamed
%    % this to fcn_MapGen_voronoiTiling. Created output from voronoiTiling
%    % that preserves the polytopes and seedPoints for each generator
%    % function
% -- Fixed a bug where corners of AABBs are not being tiled in voronoiTiling
% -- Fixed bug in script_test...voronoiTiling. Error was thrown due to how
%    % corners of AAB were handled, where seedPoint for polytope was assumed
%    % to always be inside the polytope. Added a catch case to fix - see
%    % above
% -- Deprecated fcn_MapGen_checkIfPointInsideConvexPolytope, using
%    % inpolytope instead


% TO-DO:
% -- add debug library utility, and switch functions to this. Done?
% -- add functions that, given a map, give core statistics (look out limit, linear density, etc - basically make functions to calculate all the pi-values and interpretations we might need)
% -- add prior work on grid-based map generation
% 2025_07_03 - Sean Brennan
% -- need codes to generate non-convex obstacles randomly
%    Possible approach: generate convex polytopes, and then carve
%    subpolytopes out of these
% -- need codes to generate 3D obstacles randomly via Halton set
% 2025_07_11 - Sean Brennan
% -- in fcn_MapGen_generatePolysFromTiling, seems all arguments are
%    optional. Need to fix this
% -- make sure all function calls internal to functions have -1 speed set
%    % for figure number
% -- rewrite plotPolytopes using variable input arguments (see plotRoad
%    % library?). 
%    % Then, fix call in fcn_MapGen_voronoiTiling to plot both all Voronoi
%    cells and then all each individual Voronoi cell for each generator,
%    with colors matched. Be sure to pass out plot handle (not figure
%    handle) from plotting function. Also be sure to label DisplayName of
%    plot to allow legends
% -- rename fcn_MapGen_fillPolytopeFieldsFromVertices to be:
%    % fcn_MapGen_polytopesFillFromVertices
% -- rename fcn_MapGen_generateOneRandomPolytope to be:
%    % fcn_MapGen_polytopeGenerateOneRandomPoly
% -- copy example script out of
% fcn_MapGen_generatePolysFromVoronoiAABBWithTiling into here!
% 2025_07_09 - S. Brennan and K. Hayes
% --  need a tool to check if polytope is convex in 
%     % fcn_MapGen_fillPolytopeFieldsFromVertices. This is causing some of
%     % the codes in Bounded_AStar to break

clear library_name library_folders library_url

ith_library = 1;
library_name{ith_library}    = 'DebugTools_v2025_07_10';
library_folders{ith_library} = {'Functions','Data'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/archive/refs/tags/DebugTools_v2025_07_10.zip';

ith_library = ith_library+1;
library_name{ith_library}    = 'PathClass_v2025_07_06';
library_folders{ith_library} = {'Functions'};                                
library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary/archive/refs/tags/PathClass_v2025_07_06.zip?raw=true';

% ith_library = ith_library+1;
% library_name{ith_library}    = 'GPSClass_v2023_04_21';
% library_folders{ith_library} = {''};
% library_url{ith_library}     = 'https://github.com/ivsg-psu/FieldDataCollection_GPSRelatedCodes_GPSClass/archive/refs/tags/GPSClass_v2023_04_21.zip';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'GetUserInputPath_v2023_02_01';
% library_folders{ith_library} = {''};
% library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_PathTools_GetUserInputPath/archive/refs/tags/GetUserInputPath_v2023_02_01.zip?raw=true';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'AlignCoordinates_2023_03_29';
% library_folders{ith_library} = {'Functions'};
% library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_GeomTools_AlignCoordinates/archive/refs/tags/AlignCoordinates_2023_03_29.zip?raw=true';


%% Clear paths and folders, if needed
if 1==1
    clear flag_MapGen_Folders_Initialized
    fcn_INTERNAL_clearUtilitiesFromPathAndFolders;

end

%% Do we need to set up the work space?
if ~exist('flag_MapGen_Folders_Initialized','var')
    this_project_folders = {'Functions','testFixtures','Data'}; % {'Functions','Data'};
    fcn_INTERNAL_initializeUtilities(library_name,library_folders,library_url,this_project_folders);  
    flag_MapGen_Folders_Initialized = 1;
end

%% Set environment flags for input checking in HSOV library
% These are values to set if we want to check inputs or do debugging
setenv('MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS','1');
setenv('MATLABFLAG_MAPGEN_FLAG_DO_DEBUG','0');

%% Set environment flags for input checking in Geometry library
% setenv('MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS','0');
% setenv('MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG','0');

% %% Set environment flags that define the ENU origin
% % This sets the "center" of the ENU coordinate system for all plotting
% % functions
% 
% % Location for Test Track base station
% setenv('MATLABFLAG_PLOTROAD_REFERENCE_LATITUDE','40.86368573');
% setenv('MATLABFLAG_PLOTROAD_REFERENCE_LONGITUDE','-77.83592832');
% setenv('MATLABFLAG_PLOTROAD_REFERENCE_ALTITUDE','344.189');
% 
% 
% %% Set environment flags for plotting
% % These are values to set if we are forcing image alignment via Lat and Lon
% % shifting, when doing geoplot. This is added because the geoplot images
% % are very, very slightly off at the test track, which is confusing when
% % plotting data
% setenv('MATLABFLAG_PLOTROAD_ALIGNMATLABLLAPLOTTINGIMAGES_LAT','-0.0000008');
% setenv('MATLABFLAG_PLOTROAD_ALIGNMATLABLLAPLOTTINGIMAGES_LON','0.0000054');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%    _____      _   _   _                _____ _             _           _ 
%   / ____|    | | | | (_)              / ____| |           | |         | |
%  | |  __  ___| |_| |_ _ _ __   __ _  | (___ | |_ __ _ _ __| |_ ___  __| |
%  | | |_ |/ _ \ __| __| | '_ \ / _` |  \___ \| __/ _` | '__| __/ _ \/ _` |
%  | |__| |  __/ |_| |_| | | | | (_| |  ____) | || (_| | |  | ||  __/ (_| |
%   \_____|\___|\__|\__|_|_| |_|\__, | |_____/ \__\__,_|_|   \__\___|\__,_|
%                                __/ |                                     
%                               |___/                                      
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Getting%20Started
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Welcome to the MapGen library!')

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
clear polytopes
line_width = 3;
fig_num = 111;
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];
fcn_MapGen_plotPolytopes(polytopes,fig_num,'r-',line_width);

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
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes,fig_num);
assert(isequal(round(polytopes(1).max_radius,4),2.8284));



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
% using fcn_MapGen_generatePolysFromVoronoiAABB
fig_num = 1010;
AABB = [0 0 1 1]; % Define the axis-aligned bounding box to be within the unit square
scale = max(AABB,[],'all') - min(AABB,[],'all');
new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];



% Fill in the Halton set
% pull halton set
halton_points = haltonset(2); % Construct the halton set in 2 dimensions
points_scrambled = scramble(halton_points,'RR2'); % scramble values


% pick values from halton set
Halton_range = [1 100]; % The range of points to keep
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
% Pull out only a small number of these points
seed_points = points_scrambled(low_pt:high_pt,:);
[V,C] = voronoin(seed_points);
stretch = [1 1]; % How much to stretch the tiling after creation

% fill polytopes from tiling
[polytopes,all_vertices] = fcn_MapGen_generatePolysFromVoronoiAABB(seed_points,V,C,AABB, stretch);


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
% Demos fcn_MapGen_generatePolysFromVoronoiAABB
% and script_test_fcn_MapGen_generatePolysFromVoronoiAABB
fig_num = 10;

% pull halton set
halton_points = haltonset(2); % Generate the 2-D Halton set
points_scrambled = scramble(halton_points,'RR2'); % scramble values


% pick values from halton set
Halton_range = [1 100]; % Define number of points we want (e.g 1 to 100)
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
seed_points = points_scrambled(low_pt:high_pt,:);
[V,C] = voronoin(seed_points); % Calculate the Voronoi diagram

AABB = [0 0 1 1]; % Define the axis-aligned bounding box
stretch = [1 1]; % Define the stretch factor on each axis

% fill polytopes from tiling, and give a figure number to show results
fcn_MapGen_generatePolysFromVoronoiAABB(seed_points,V,C,AABB, stretch,fig_num);

%% (DEPRECATED) Show time-dependent generation of polytopes from a tiling
% % Moves just one point
% fig_num = 11111;
% clf;
% 
% 
% % pull halton set
% rng(1111);
% halton_points = haltonset(2);
% points_scrambled = scramble(halton_points,'RR2'); % scramble values
% AABB = [0 0 1 1]; % Define the axis-aligned bounding box
% stretch = [1 1];
% 
% 
% % pick values from halton set
% Npoints = 100;
% Halton_range = [1 Npoints];
% low_pt = Halton_range(1,1);
% high_pt = Halton_range(1,2);
% original_seed_points = points_scrambled(low_pt:high_pt,:);
% 
% % Set the point to move
% pointID_to_move = 9;
% 
% % Set the velocity vectors
% %velocities = 0.01*randn(Npoints,2);
% velocities = 0.01;
% 
% % Prep the movie
% % vidfile = VideoWriter('MovingPolytope.mp4','MPEG-4');
% vidfile = VideoWriter('MovingPolytope.avi');
% open(vidfile);
% 
% moved_seed_points = original_seed_points;
% for ith_step = 1:1:100
%     moved_seed_points(pointID_to_move) = original_seed_points(pointID_to_move) + velocities*ith_step;
% 
%     % Wrap points
%     seed_points = mod(moved_seed_points,1);
% 
%     [V,C] = voronoin(seed_points);
% 
% 
%     % fill polytopes from tiling
%     polytopes = fcn_MapGen_generatePolysFromVoronoiAABB(seed_points,V,C,AABB, stretch);
% 
%     des_gap_size = 0.02;
% 
%     shrunk_polytopes=...
%         fcn_MapGen_polytopesShrinkFromEdges(...
%         polytopes,des_gap_size);
% 
% 
%     % Plot the results
%     figure(fig_num);
%     clf;
% 
%     line_width = 2;
%     color = [0 0 0];
%     axis_box = [0 1 0 1];
%     fcn_MapGen_plotPolytopes(shrunk_polytopes,gca,'-',line_width,color,axis_box,'square'); %,[1 0 0 0 0.5]);
%     hold on;
%     highlighted_polytope = shrunk_polytopes(pointID_to_move);
% 
%     % Highlight the moving polytope
%     % FILL_INFO: a 1-by-5 vector to specify wether or not there is fill, the
%     % color of fill, and the opacity of the fill [Y/N, R, G, B, alpha]
%     % fill_info = [1 1 0 0 1];
%     fcn_MapGen_plotPolytopes(highlighted_polytope,gca,'-',line_width,[1 0 0],axis_box,'square');  % , fill_info);
% 
% 
%     frame = getframe(gcf);
%     writeVideo(vidfile,frame);
% 
% 
%     pause(0.1);
% end
% close(vidfile)


% %% Show time-dependent generation of polytopes from a tiling
% % Moves all points
% fig_num = 1111;
% clf;
% 
% 
% % pull halton set
% rng(1111);
% halton_points = haltonset(2);
% points_scrambled = scramble(halton_points,'RR2'); % scramble values
% AABB = [0 0 1 1]; % Define the axis-aligned bounding box
% stretch = [1 1];
% 
% 
% % pick values from halton set
% Npoints = 100;
% Halton_range = [1 Npoints];
% low_pt = Halton_range(1,1);
% high_pt = Halton_range(1,2);
% original_seed_points = points_scrambled(low_pt:high_pt,:);
% 
% % Set the point to move
% % pointID_to_move = 9;
% 
% % Set the velocity vectors
% velocities = 0.01*randn(Npoints,2);
% % velocities = 0.01;
% 
% % Prep the movie
% % vidfile = VideoWriter('MovingPolytope.mp4','MPEG-4');
% vidfile = VideoWriter('MovingPolytope.avi');
% open(vidfile);
% 
% % Set figure parameters
% line_width = 2;
% color = [0 0 0];
% axis_box = [0 1 0 1];
% 
% 
% moved_seed_points = original_seed_points;
% for ith_step = 1:1:100
%     % moved_seed_points(pointID_to_move) = original_seed_points(pointID_to_move) + velocities*ith_step;
%     moved_seed_points = original_seed_points + velocities*ith_step;
% 
%     % Wrap points to all exist between 0 and 1
%     seed_points = mod(moved_seed_points,1);
% 
%     % Calculate the Voronoi diagram
%     [V,C] = voronoin(seed_points);
% 
%     % fill polytopes from tiling of Voronoi diagram
%     polytopes = fcn_MapGen_generatePolysFromVoronoiAABB(seed_points,V,C,AABB, stretch);
% 
%     des_gap_size = 0.02;
% 
%     shrunk_polytopes=...
%         fcn_MapGen_polytopesShrinkFromEdges(...
%         polytopes,des_gap_size);
% 
% 
%     % Plot the results
%     figure(fig_num);
%     clf;
% 
%     fcn_MapGen_plotPolytopes(shrunk_polytopes,gca,'-',line_width,color,axis_box,'square'); %,[1 0 0 0 0.5]);
%     
%     %     % Highlight the moving polytope
%     %     % FILL_INFO: a 1-by-5 vector to specify wether or not there is fill, the
%     %     % color of fill, and the opacity of the fill [Y/N, R, G, B, alpha]
%     %     % fill_info = [1 1 0 0 1];
%     %     fcn_MapGen_plotPolytopes(highlighted_polytope,gca,'-',line_width,[1 0 0],axis_box,'square');  % , fill_info);
% 
% 
%     frame = getframe(gcf);
%     writeVideo(vidfile,frame);
% 
% 
%     pause(0.1);
% end
% close(vidfile)

%% Generate a Voronoi tiling that is a true tile
% This is done by actually tiling around the central seed points copies of
% these same points, then calculating the Voronoi diagram for this

fig_num = 112;

% Pull halton set
rng(1111);
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values
unit_AABB = [0 0 1 1]; % Define the axis-aligned bounding box
stretch = [1 1];

% Pick values from halton set
Npoints = 100;
Halton_range = [1 Npoints];
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
original_seed_points = points_scrambled(low_pt:high_pt,:);

% Call the function
polytopes = fcn_MapGen_generatePolysFromVoronoiAABBWithTiling(original_seed_points,unit_AABB, stretch,fig_num);
assert(isstruct(polytopes));



%% Show time-dependent generation of polytopes from a tiling
% Moves just one point
fig_num = 11111;
clf;


% pull halton set
rng(1111);
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values
% AABB = [0 0 1 1]; % Define the axis-aligned bounding box
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

flag_make_video = 1;
% Prep the movie?

des_gap_size = 0.02;

moved_seed_points = original_seed_points;
Nsteps = 100;
im{Nsteps} = [];
for ith_step = 1:Nsteps
    moved_seed_points(pointID_to_move) = original_seed_points(pointID_to_move) + velocities*ith_step;
    
    % Wrap points
    seed_points = mod(moved_seed_points,1);
       
    % fill polytopes from tiling
    polytopes = fcn_MapGen_generatePolysFromVoronoiAABBWithTiling(seed_points,unit_AABB, stretch);
        
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
    
    if flag_make_video
        frame = getframe(gcf);
        im{ith_step} = frame2im(frame);
    end
    
    

    pause(0.1);
end

if flag_make_video
    delay_time = 1/25;
    filename = "testAnimated.gif"; % Specify the output file name
    for idx = 1:Nsteps
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",delay_time);
        else
            imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",delay_time);
        end
    end
end


%% Show time-dependent generation of polytopes from a tiling
% Moves all points
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
velocities = 0.001*randn(Npoints,2);
% velocities = 0.01;


flag_make_video = 1;
% Prep the movie?

des_gap_size = 0.02;

Nsteps = 1000;
im{Nsteps} = [];
for ith_step = 1:Nsteps
    fprintf(1,'Iteration: %.0d\n',ith_step);
    % moved_seed_points(pointID_to_move) = original_seed_points(pointID_to_move) + velocities*ith_step;
    moved_seed_points = original_seed_points + velocities*ith_step;
    
    % Wrap points
    seed_points = mod(moved_seed_points,1);
       
    % fill polytopes from tiling
    polytopes = fcn_MapGen_generatePolysFromVoronoiAABBWithTiling(seed_points,unit_AABB, stretch);
        
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
    
    if flag_make_video
        frame = getframe(gcf);
        im{ith_step} = frame2im(frame);
    end
    
    

    pause(0.1);
end

if flag_make_video
    delay_time = 1/25;
    filename = "testAnimated.gif"; % Specify the output file name
    for idx = 1:Nsteps
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",delay_time);
        else
            imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",delay_time);
        end
    end
end



%% Generate a set of polytopes from various pseudo-random sources
close all;

% Generate a set of polytopes from the Sobol set
fig_num = 333;
figure(fig_num);
clf;
hold on;
Numpoints = 100;

subplot(2,3,1);
% Sobol_range = [1 Numpoints]; % range of Sobol points to use to generate the tiling
% tiled_polytopes = fcn_MapGen_sobolVoronoiTiling(Sobol_range,[1 1],fig_num);

seedGeneratorNames = 'sobolset';
seedGeneratorRanges = [1 Numpoints];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[~] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num)); 

title('Sobol set');
legend('off');

% Generate a set of polytopes from the Halton set
subplot(2,3,2);
seedGeneratorNames = 'haltonset';
seedGeneratorRanges = [1 Numpoints];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[~] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num)); 

%fcn_MapGen_plotPolytopes(tiled_polytopes,gca,'-',line_width,[0 0 1],axis_box,'square');  % , fill_info);
title('Halton set');
legend('off');

% Generate a set of polytopes from the Latin Hypercube set
subplot(2,3,3);
seedGeneratorNames = 'latin';
seedGeneratorRanges = [1 Numpoints];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[~] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num)); 

title('Latin Hypercube set');
legend('off');

% Generate a set of polytopes from the Random set
subplot(2,3,4);
seedGeneratorNames = 'randn';
seedGeneratorRanges = [1 Numpoints];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[~] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num)); 

title('Uniform random set');
legend('off');

% Generate a set of polytopes from the Random Normal set
subplot(2,3,5);
seedGeneratorNames = 'rand';
seedGeneratorRanges = [1 Numpoints];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[~] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num)); 

title('Random normally distributed set');
legend('off');

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

polytopes = fcn_MapGen_mixedSetVoronoiTiling(mixedSet,stretch,fig_num); %#ok<NASGU>


%% Generate many test sets of polytopes from the Halton set
for i=1:100:500
    fig_num = 21+i;
    figure(fig_num); clf;
    
    seedGeneratorNames = 'haltonset';
    seedGeneratorRanges = [i i+100];
    AABBs = [0 0 1 1];
    mapStretchs = [1 1];
    [~] = fcn_MapGen_voronoiTiling(...
        seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
        seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
        (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
        (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
        (fig_num)); 

    % Do statistics, checking that the area is always fully filled and we
    % get 101 polytopes each time

    % temp = fcn_MapGen_polytopesStatistics(...
    % polytopes);
    % title(sprintf('Halton range is: [%.0d %.0d]',i,i+100));
    % assert(abs(temp.unoccupancy_ratio)<(1000*eps));
    % assert(isequal(101,temp.point_density));
    % pause(0.1);
end

%% Show how the maps can be trimmed to a box

% Generate polytopes from the Halton set

fig_num = 31;
seedGeneratorNames = 'haltonset';
seedGeneratorRanges = [5401 5501];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[tiled_polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num));


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
shrunk_polytopes2 = fcn_MapGen_polytopesShrinkToRadius(...
    trimmed_polytopes,des_rad,sigma_radius,min_rad,fig_num); %#ok<NASGU>

%% Show how different shrinking methods change statistics
fig_num = 555;

% Generate polytopes from the Halton set        
seedGeneratorNames = 'haltonset';
seedGeneratorRanges = [5401 5501];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[tiled_polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num));

    
% Grab statistics on original map
fcn_MapGen_polytopesStatistics(...
    tiled_polytopes,...
    fig_num+1);

% Shrink to radius
fig_num = 556;
des_rad = 0.03; sigma_radius = 0; min_rad = 0.001;
shrunk_polytopes2=fcn_MapGen_polytopesShrinkToRadius(...
    trimmed_polytopes,des_rad,sigma_radius,min_rad,fig_num); %#ok<NASGU>


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
exp_polytopes=fcn_MapGen_polytopesExpandEvenly(one_polytope,exp_dist,fig_num); %#ok<NASGU>

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
seedGeneratorNames = 'haltonset';
seedGeneratorRanges = [1 1000]; % range of Halton points to use to generate the tiling
AABBs = [0 0 1 1];
mapStretchs = [200 200]; % stretch in the x and y directions
[polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num));

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
% h=colorbar;
title('Error in X Dimension [m]')
xlabel('X [m]')
ylabel('Y [m]')

figure('name','Error in Y')
contour(x,y,err.y,'ShowText','on')
% h=colorbar;
title('Error in Y Dimension [m]')
xlabel('X [m]')
ylabel('Y [m]')

figure('name','Total Error (Bubble Radius) [m]')
contour(x,y,err.R,'ShowText','on')
% h=colorbar;
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

%% function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
% Clear out the variables
clear global flag* FLAG*
clear flag*
clear path

% Clear out any path directories under Utilities
path_dirs = regexp(path,'[;]','split');
utilities_dir = fullfile(pwd,filesep,'Utilities');
for ith_dir = 1:length(path_dirs)
    utility_flag = strfind(path_dirs{ith_dir},utilities_dir);
    if ~isempty(utility_flag)
        rmpath(path_dirs{ith_dir});
    end
end

% Delete the Utilities folder, to be extra clean!
if  exist(utilities_dir,'dir')
    [status,message,message_ID] = rmdir(utilities_dir,'s');
    if 0==status
        error('Unable remove directory: %s \nReason message: %s \nand message_ID: %s\n',utilities_dir, message,message_ID);
    end
end

end % Ends fcn_INTERNAL_clearUtilitiesFromPathAndFolders

%% fcn_INTERNAL_initializeUtilities
function  fcn_INTERNAL_initializeUtilities(library_name,library_folders,library_url,this_project_folders)
% Reset all flags for installs to empty
clear global FLAG*

fprintf(1,'Installing utilities necessary for code ...\n');

% Dependencies and Setup of the Code
% This code depends on several other libraries of codes that contain
% commonly used functions. We check to see if these libraries are installed
% into our "Utilities" folder, and if not, we install them and then set a
% flag to not install them again.

% Set up libraries
for ith_library = 1:length(library_name)
    dependency_name = library_name{ith_library};
    dependency_subfolders = library_folders{ith_library};
    dependency_url = library_url{ith_library};

    fprintf(1,'\tAdding library: %s ...',dependency_name);
    fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url);
    clear dependency_name dependency_subfolders dependency_url
    fprintf(1,'Done.\n');
end

% Set dependencies for this project specifically
fcn_DebugTools_addSubdirectoriesToPath(pwd,this_project_folders);

disp('Done setting up libraries, adding each to MATLAB path, and adding current repo folders to path.');
end % Ends fcn_INTERNAL_initializeUtilities


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
% 2023_04_20:
% -- improved error handling
% -- fixes nested installs automatically

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
    % Save the root directory, so we can get back to it after some of the
    % operations below. We use the Print Working Directory command (pwd) to
    % do this. Note: this command is from Unix/Linux world, but is so
    % useful that MATLAB made their own!
    root_directory_name = pwd;

    % Does the directory "Utilities" exist?
    utilities_folder_name = fullfile(root_directory_name,'Utilities');
    if ~exist(utilities_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(root_directory_name,'Utilities');

        % Did it work?
        if ~success_flag
            error('Unable to make the Utilities directory. Reason: %s with message ID: %s\n',error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The Utilities directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',error_message, message_ID);
        end

    end

    % Does the directory for the dependency folder exist?
    dependency_folder_name = fullfile(root_directory_name,'Utilities',dependency_name);
    if ~exist(dependency_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(utilities_folder_name,dependency_name);

        % Did it work?
        if ~success_flag
            error('Unable to make the dependency directory: %s. Reason: %s with message ID: %s\n',dependency_name, error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The %s directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',dependency_name, error_message, message_ID);
        end

    end

    % Do the subfolders exist?
    flag_allFoldersThere = 1;
    if isempty(dependency_subfolders{1})
        flag_allFoldersThere = 0;
    else
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

    % Do we need to unzip the files?
    if flag_allFoldersThere==0
        % Files do not exist yet - try unzipping them.
        save_file_name = tempname(root_directory_name);
        zip_file_name = websave(save_file_name,dependency_url);
        % CANT GET THIS TO WORK --> unzip(zip_file_url, debugTools_folder_name);

        % Is the file there?
        if ~exist(zip_file_name,'file')
            error(['The zip file: %s for dependency: %s did not download correctly.\n' ...
                'This is usually because permissions are restricted on ' ...
                'the current directory. Check the code install ' ...
                '(see README.md) and try again.\n'],zip_file_name, dependency_name);
        end

        % Try unzipping
        unzip(zip_file_name, dependency_folder_name);

        % Did this work? If so, directory should not be empty
        directory_contents = dir(dependency_folder_name);
        if isempty(directory_contents)
            error(['The necessary dependency: %s has an error in install ' ...
                'where the zip file downloaded correctly, ' ...
                'but the unzip operation did not put any content ' ...
                'into the correct folder. ' ...
                'This suggests a bad zip file or permissions error ' ...
                'on the local computer.\n'],dependency_name);
        end

        % Check if is a nested install (for example, installing a folder
        % "Toolsets" under a folder called "Toolsets"). This can be found
        % if there's a folder whose name contains the dependency_name
        flag_is_nested_install = 0;
        for ith_entry = 1:length(directory_contents)
            if contains(directory_contents(ith_entry).name,dependency_name)
                if directory_contents(ith_entry).isdir
                    flag_is_nested_install = 1;
                    install_directory_from = fullfile(directory_contents(ith_entry).folder,directory_contents(ith_entry).name);
                    install_files_from = fullfile(directory_contents(ith_entry).folder,directory_contents(ith_entry).name,'*.*');
                    install_location_to = fullfile(directory_contents(ith_entry).folder);
                end
            end
        end

        if flag_is_nested_install
            [status,message,message_ID] = movefile(install_files_from,install_location_to);
            if 0==status
                error(['Unable to move files from directory: %s\n ' ...
                    'To: %s \n' ...
                    'Reason message: %s\n' ...
                    'And message_ID: %s\n'],install_files_from,install_location_to, message,message_ID);
            end
            [status,message,message_ID] = rmdir(install_directory_from);
            if 0==status
                error(['Unable remove directory: %s \n' ...
                    'Reason message: %s \n' ...
                    'And message_ID: %s\n'],install_directory_from,message,message_ID);
            end
        end

        % Make sure the subfolders were created
        flag_allFoldersThere = 1;
        if ~isempty(dependency_subfolders{1})
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
         % If any are not there, then throw an error
        if flag_allFoldersThere==0
            error(['The necessary dependency: %s has an error in install, ' ...
                'or error performing an unzip operation. The subfolders ' ...
                'requested by the code were not found after the unzip ' ...
                'operation. This suggests a bad zip file, or a permissions ' ...
                'error on the local computer, or that folders are ' ...
                'specified that are not present on the remote code ' ...
                'repository.\n'],dependency_name);
        else
            % Clean up the zip file
            delete(zip_file_name);
        end

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
        debugTools_function_folder = fullfile(root_directory_name, 'Utilities', dependency_name,'Functions');

        % Move into the folder, run the function, and move back
        cd(debugTools_function_folder);
        fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        cd(root_directory_name);
    else
        try
            fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        catch
            error(['Package installer requires DebugTools package to be ' ...
                'installed first. Please install that before ' ...
                'installing this package']);
        end
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

