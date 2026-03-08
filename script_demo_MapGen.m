% script_demo_MapGen.m
% This is a script that shows the capabilities of the "MapGen" class of
% functions via demonstrations.

% REVISION HISTORY:
% 
% 2021_06_07:
% - First write of the function, using the "Vis" library demo script as
%   % starter
% 
% 2021_06_09
% - Added other types of point generators
% 
% 2021_07_06
% - Updated to include the newer expansion functions
% 
% 2021_07_11
% - Add ability to extend halton set to right (e.g. "scrolling" map), see
%   % the function: fcn_MapGen_mixedSetVoronoiTiling
% 
% 2021_07_12
% - Added ability to determine generic map statistics via the function:
%   % fcn_MapGen_statsPolytopes
% 
% 2023_01_15
% - Added demo of edge-based shrinking
% 
% 2023_02_20
% - Added code to better support README.md
% 
% 2023_02_21
% - Added Debug utility library
% 
% 2023_03_13 
% - Merged changes allowing for repeated tiling polytope field
% 
% 2023_04_27
% - Updated the installer to latest version of Debug Tools (fixes bug with
%   % Git archives of zips).
% 
% 2023_05_04 - sbrennan@psu.edu
% - Cleared the path variable, in case variable of same name shadows
%   % command. This was causing errors in some codes.
% - functionalized the clear process
% 
% 2025_04_25 - Sean Brennan
% - Updated DebugTools_v2024_12_18 dependency
% - Added global flags for setting test conditions and plotting
% - Deprecated fcn_MapGen_checkInputsToFunctions, converted to fcn_DebugTools_checkInputsToFunctions
% 
% 2025_06_25 - Sean Brennan
% - upgraded debug tools to version DebugTools_v2025_06_24
% 
% 2025_07_03 - Sean Brennan
% - added minor comments in demo script
% 
% 2025_07_07 - Sean Brennan
% - started updating headers and test scripts. Only have 
%   % fcn_MapGen_generatePolysFromTiling done so far.
% 
% 2025_07_11 - Sean Brennan
% - updated DebugTools library to DebugTools_v2025_07_15
% - added PathClass library to use fcn_Path_findSensorHitOnWall function
%   % rather than
%   % fcn_MapGen_findIntersectionOfSegments
%   % the Path library version is FAR better in speed, debugged, support,
%   % etc. Deprecated fcn_MapGen_findIntersectionOfSegments
% - for all the tiling variants, deprecated the following
%   % fcn_MapGen_sobolVoronoiTiling
%   % fcn_MapGen_latinVoronoiTiling
%   % (etc)
%   % Merged these to use fcn_MapGen_mixedSetVoronoiTiling and renamed
%   % this to fcn_MapGen_generatePolysFromSeedGeneratorNames. Created output from voronoiTiling
%   % that preserves the polytopes and seedPoints for each generator
%   % function
% - Fixed a bug where corners of AABBs are not being tiled in voronoiTiling
% - Fixed bug in script_test...voronoiTiling. Error was thrown due to how
%   % corners of AAB were handled, where seedPoint for polytope was assumed
%   % to always be inside the polytope. Added a catch case to fix - see
%   % above
% - Deprecated fcn_MapGen_checkIfPointInsideConvexPolytope, using
%   % inpolytope instead
% - rewrote plotPolytopes using variable input arguments (see plotRoad
%   % library) to allow arbitrary formatting. Deprecated old version.
% - added option in fcn_MapGen_generatePolysFromSeedGeneratorNames to plot both all Voronoi
%   % cells and then all each individual Voronoi cell for each generator,
%   % with colors matched. 
% - fixed bugs in fcn_MapGen_verticesCropToAABB where infinite vertices
%   % not treated correctly
% - fully used DebugTools library utility for input checking
% - make sure all function calls internal to functions have -1 speed set
%   % for figure number
% - renamed fcn_MapGen_generateOneRandomPolytope to be:
%   % fcn_MapGen_polytopeGenerateOneRandomPoly
% - renamed polytopes function:
%   % fcn_MapGen_fillPolytopeFieldsFromVertices to fcn_MapGen_polytopesFillFieldsFromVertices
%   % fcn_MapGen_increasePolytopeVertexCount to fcn_MapGen_polytopesIncreaseVertexCount
% - renamed AABB functions:
%   % fcn_MapGen_convertAABBtoWalls to fcn_MapGen_AABBConvertToWalls
%   % fcn_MapGen_snapToAABB to fcn_MapGen_AABBsnapTo
%   % fcn_MapGen_isWithinABBB to fcn_MapGen_AABBisWithin
%   % fcn_MapGen_projectVectorToAABB to fcn_MapGen_AABBprojectVectorTo
% - renamed vertices functions:
%   % fcn_MapGen_cropPolytopeToRange to fcn_MapGen_verticesCropToAABB
%   % fcn_MapGen_cropVerticesByWallIntersections to fcn_MapGen_verticesCropToWallIntersections
%   % fcn_MapGen_removeInfiniteVertices to fcn_MapGen_verticesRemoveInfinite
% 
% 2025_07_17 by Sean Brennan
% - renamed more vertices functions:
%   % fcn_MapGen_tilePoints to fcn_MapGen_verticesTiling
% - renamed stats functions:
%   % fcn_MapGen_calculateConvexHullOverlapRatio to fcn_MapGen_statsConvexHullOverlapRatio
%   % fcn_MapGen_polytopesStatistics to fcn_MapGen_statsPolytopes
% - renamed polytopes function:
%   % fcn_MapGen_flattenPolytopeMap to fcn_MapGen_polytopesFlattenMap
% - renamed generator functions:
%   % fcn_MapGen_voronoiTiling to fcn_MapGen_generatePolysFromSeedGeneratorNames
%   % fcn_MapGen_nameToMap to fcn_MapGen_generatePolysFromName
%
% 2025_07_18 - S. Brennan
% - Added GridMapGen tools including:
%   % * fcn_GridMapGen_dilationOccupancyStats: calcs occupancy stats
%   % * fcn_GridMapGen_dilateByN: dilates a matrix by N cells
%   % * fcn_GridMapGen_dilateOccupancyByN: dilates occupancy by N cells 
%   % * fcn_GridMapGen_generateRandomOccupancyMap: generates a random occupancy map
% - Added GridMapGen demos including:
%   % * script_demo_dilateCompareDilateByNSpeeds - shows the dilation
%   %     toolsets
%   % * script_demo_generateRandomOccupancyAnimated - animates random
%   %     blobs
% 
% 2025_07_26 - S. Brennan
% - Added path generation examples from course material
% 
% 2025_07_26 by Sean Brennan
% - fixed fcn_MapGen_polytopeCropEdges
%   % * renamed to fcn_MapGen_polytopesDeleteByAABB, for consistency
%   % * changed input to AABB style, for consistency
% - fixed fcn_MapGen_polytopeMapGen 
%   % * renamed inputs to remove underscores
%   % * changed boundingBox input to AABB format, for consistency
% 
% 2025_07_28 by Sean Brennan
% - added and tested script_test_all_GridMapGen_functions to GridMapGen
% 
% 2025_07_28 by Sean Brennan
% - fcn_MapGen_generatePolysFromSeedGeneratorNames, fixed missing defaults
%
% 2025_07_29 by Sean Brennan
% - fcn_MapGen_polytopeShrinkEvenly created 
%   % * copied from fcn_MapGen_polytopesShrinkEvenly, which was misnamed
%   % * added test script, script_test_fcn_MapGen_polytopeShrinkEvenly
%   % DEPRECATED: fcn_MapGen_polytopeShrinkFromEdges, which is same fcn
% - fcn_MapGen_polytopesShrinkEvenly created
%   % * function fcn_MapGen_polytopesShrinkToGapSize created from prior ver
%   % * prior version was mis-named!
%   % * script_demo_fcn_MapGen_polytopesShrinkToGapSize created from prior
%   % * added test script: script_test_fcn_MapGen_polytopesShrinkEvenly
%
% (new release)
%
% 2025_10_02 by Sean Brennan
% - Updated DebugTools_v2025_09_26b
% - Updated PathClass_v2025_08_03
% - In fcn_MapGen_polytopesShrinkEvenly
%   % Fixed header docstrings to correctly show dependence on 
%   % fcn_MapGen_polytopeShrinkEvenly
% - In fcn_MapGen_polytopeShrinkEvenly  
%   % Fixed header docstrings to correctly show dependence on
%
% 2025_11_06 by Sean Brennan
% - In fcn_MapGen_plotPolytopes
%   % * removed duplicate figure() call within Inputs area (not needed)
% - Moved GUI script to different directory - doesn't belong in Functions
% - Slight updates to README.md, this needs more work
% - Updated script_test_all_functions to latest version
%   % * copied from DebugTools
% - Updated DebugTools_v2025_11_04c
% - Renamed this main function to match repo shortname more clearly
% - Added fcn_MapGen_loadTestMap
%   % * Pulled this function out of BoundedAStar
% - Removed deprecation warnings on fcn_MapGen_polytopeFindVertexSkeleton
% - Updated script_test_all_functions
% (new release)
%
% 2025_11_18 by Sean Brennan
% - Updated revision history to be in Markdown format
% - In fcn_MapGen_plotPolytopes
%   % * updated docstrings to better explain polytopes structure
% (new release)
%
% 2025_11_20 - S. Brennan
% - Set up new installer
% - In all functions
%   % * Updated rev history to be in Markdown format
%   % * Replaced figNum with figNum
%   % * Standardized all  % REVISION HISTORY: and % TO-DO: fields
% - In fcn_MapGen_plotPolytopes
%   % * fixed bug where polytope plotting not plotting closed form in case
%   %   % where input polytopes are not closed off
%   % * Updated rev history to be in Markdown format
%   % * Updated rescaling on plotting to use userdata rather than children
%   % * Updated auto axes to use padding methods from MATLAB
%   % * Replaced fig_+num with figNum
% - Updated script_test_all_functions - ran script to confirm repo
% - Fixed main script issues
%   % * deprecated plotting being used
%   % * deprecated shrinking calls
%   % * deprecated multiple tiling types
% - Fixed deprecated warnings to point to correct files:
%   % * fcn_MapGen_mixedSetVoronoiTiling
%   % * fcn_MapGen_latinVoronoiTiling
%   % * fcn_MapGen_randVoronoiTiling
%   % * fcn_MapGen_randomNormalVoronoiTiling
%   % * fcn_MapGen_sobolVoronoiTiling
%   % * fcn_MapGen_voronoiTiling
% - Fixed bug in fcn_MapGen_polytopeShrinkEvently
%   % * preallocation of polytope structure was producing
%   %   % dissimilar structure fields
% (new release)


% TO-DO:
% - Modify statPoly code to give core statistics including:
%   % look out limit, 
%   % linear density, etc
%   % basically make function to calculate all the pi-values and 
%   % interpretations we might need for publications
% - add prior work on grid-based map generation
% 
% 2025_07_03 - Sean Brennan
% - need codes to generate non-convex obstacles randomly
%   % Possible approach: generate convex polytopes, and then carve
%   % subpolytopes out of these
% - need codes to generate 3D obstacles randomly via Halton set
% 
% 2025_07_09 - S. Brennan and K. Hayes
% - need a tool to check if polytope is convex in 
%   % fcn_MapGen_polytopesFillFieldsFromVertices. This is causing some of
%   % the codes in Bounded_AStar to break.
%   % UPDATE: 2025_11_06 - S. Brennan - there's a new function in Geom
%   % now to do this: See fcn_geometry_findPolytopeOrientations
% 
% 2025_07_11 - Sean Brennan
% - in fcn_MapGen_generatePolysFromTiling, seems all arguments are
%   % optional. Need to fix this
% - compare example script out of
%   % fcn_MapGen_generatePolysFromVoronoiAABBWithTiling to code in this main
%   % demo. remove the code if not being used anymore.
% - rewrite fcn_MapGen_polytopesIncreaseVertexCount to use linspace between
%   % X and Y vertices rather than interpolation. MUCH simpler.
% - rework fcn_MapGen_polytopesPredictLengthCostRatio - code is very messy
% - need to finish script_test_fcn_MapGen_AABBprojectVectorTo
% - figure out difference between fcn_MapGen_generatePolysFromTiling and
%   % fcn_MapGen_generatePolysFromVoronoiAABB - they look the same
% 
% 2025_07_29 - S. Brennan
% - Need to add a test script for script_test_fcn_MapGen_polytopesShrinkToGapSize
% 
% 2025_07_29 by Sean Brennan
% - fcn_MapGen_polytopeFindVertexSkeleton needs to be replaced with VSkel
%   % library function
% 
% 2025_11_06 - S. 
% - Need to finish updating README. Note: the load function (see 2025_11_06
%   % revision notes, needs to be updated.
% - Need to finish VSkel library and then deprecate
%   % fcn_MapGen_polytopeFindVertexSkeleton (it's deprecated now, but warning
%   % is shut off), then search/replace for fcn_MapGen_polytopeFindVertexSkeleton

%% Check which files contain key strings?
if 1==0
    clc
    functionsDirectoryQuery = fullfile(pwd,'Functions','*.*');
    % Use the following instead, if wish to do subdirectories
    % directoryQuery = fullfile(pwd,'Functions','**','*.*');

    fileListFunctionsFolder = dir(functionsDirectoryQuery); %cat(2,'.',filesep,filesep,'script_test_fcn_*.m'));

    % Filter out directories from the list
    fileListFunctionsFolderNoDirectories = fileListFunctionsFolder(~[fileListFunctionsFolder.isdir]);

    % Make sure there's not figNum
    queryString = '% --';
    flagsStringWasFoundInFilesRaw = fcn_DebugTools_directoryStringQuery(fileListFunctionsFolderNoDirectories, queryString, (-1));
    % flagsStringWasFoundInFiles = fcn_INTERNAL_removeFromList(flagsStringWasFoundInFilesRaw, fileListFunctionsFolderNoDirectories,'script_test_all_functions');
    if sum(flagsStringWasFoundInFilesRaw)>0
        fcn_DebugTools_directoryStringQuery(fileListFunctionsFolderNoDirectories, queryString, 1);
        error('A "%s" was found in one of the functions - see listing above.', queryString);
    end
end


%% Make sure we are running out of root directory
st = dbstack; 
thisFile = which(st(1).file);
[filepath,name,ext] = fileparts(thisFile);
cd(filepath);

%% Clear paths and folders, if needed
if 1==1
    clear flag_MapGen_Folders_Initialized
end

if 1==0
    fcn_INTERNAL_clearUtilitiesFromPathAndFolders;
end

if 1==0
    % Resets all paths to factory default
    restoredefaultpath;
end

%% Install dependencies
% Define a universal resource locator (URL) pointing to the repos of
% dependencies to install. Note that DebugTools is always installed
% automatically, first, even if not listed:
clear dependencyURLs dependencySubfolders
ith_repo = 0;

% ith_repo = ith_repo+1;
% dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/PathPlanning_MapTools_MapGenClassLibrary';
% dependencySubfolders{ith_repo} = {'Functions','testFixtures','GridMapGen'};

ith_repo = ith_repo+1;
dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary';
dependencySubfolders{ith_repo} = {'Functions','Data'};

% ith_repo = ith_repo+1;
% dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/FieldDataCollection_VisualizingFieldData_PlotRoad';
% dependencySubfolders{ith_repo} = {'Functions','Data'};
% 
% ith_repo = ith_repo+1;
% dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/PathPlanning_GeomTools_GeomClassLibrary';
% dependencySubfolders{ith_repo} = {'Functions','Data'};

%% Do we need to set up the work space?
if ~exist('flag_MapGen_Folders_Initialized','var')

    % Clear prior global variable flags
    clear global FLAG_*

    % Navigate to the Installer directory
    currentFolder = pwd;
    cd('Installer');
    % Create a function handle
    func_handle = @fcn_DebugTools_autoInstallRepos;

    % Return to the original directory
    cd(currentFolder);

    % Call the function to do the install
    func_handle(dependencyURLs, dependencySubfolders, (0), (-1));

    % Add this function's folders to the path
    this_project_folders = {'Functions','testFixtures','Data', 'GridMapGen'};
    fcn_DebugTools_addSubdirectoriesToPath(pwd,this_project_folders)

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
% Single polytope operations are functions that usually start with
% "polytope"

%% Show how to check if points are within an Axis-Aligned Bounding Box
AABB = [0 0 1 1]; % Define the axis-aligned bounding box
test_points = randn(100,2);
figNum = 1;
isInside = fcn_MapGen_AABBisWithin(AABB,test_points,figNum);


%% fcn_MapGen_plotPolytopes
figNum = 10001;
titleString = sprintf('DEMO case: basic plotting demo');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

clear polytopes
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];

plotFormat.LineWidth = 3;
plotFormat.MarkerSize = 10;
plotFormat.LineStyle = '-';

fillFormat = [];

% Call the function
h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat),(fillFormat),(figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(ishandle(h_plot));

% Check variable sizes
assert(size(h_plot,1)==1);
assert(size(h_plot,2)==1);

% Check variable values
assert(isequal(h_plot.Parent.Parent.Number, figNum));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% Show how we calculate the polytope centroid and area
% Note: this does NOT have to be a convex polytope, as the example shows.
% Copied from script_test_fcn_MapGen_polytopeCentroidAndArea
% Tests: fcn_MapGen_polytopeCentroidAndArea

figNum = 2;

x = [3; 4; 2; -1; -2; -3; -4; -2; 1; 2; 3];
y = [1; 2; 2; 3; 2; -1; -2; -3; -3; -2; 1];
[Centroid,Area] = fcn_MapGen_polytopeCentroidAndArea([x,y],figNum);

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


figNum = 3;
clear polytopes
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];
polytopes = fcn_MapGen_polytopesFillFieldsFromVertices(polytopes,figNum);
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
figNum = 1010;
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
figure(figNum);
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
plotFormat.LineWidth = 2;
plotFormat.MarkerSize = 10;
plotFormat.LineStyle = '-';
plotFormat.Color = [0 0 1];

fillFormat = [1 0 0 0 0.5];

% Call the function
fcn_MapGen_plotPolytopes(polytopes, (plotFormat),(fillFormat),(figNum));

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
% plot the polytopes on current axis
plotFormat.LineWidth = 2;
plotFormat.MarkerSize = 10;
plotFormat.LineStyle = '-';
plotFormat.Color = [0 0 1];

fillFormat = [1 0 0 0 0.5];

% Call the function
fcn_MapGen_plotPolytopes(polytopes, (plotFormat),(fillFormat),(figNum));

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
plotFormat.LineWidth = 2;
plotFormat.MarkerSize = 10;
plotFormat.LineStyle = '-';
plotFormat.Color = [1 0 0];

fillFormat = [1 0 0 0 0.5];

% Call the function
fcn_MapGen_plotPolytopes(polytopes, (plotFormat),(fillFormat),(figNum));
hold on;

axis(new_axis);
title('Shrunk to radius polytopes');

% plot the shrink to edge
subplot(2,3,6);

des_gap_size = 0.05;


shrunk_polytopes2 = ...
    fcn_MapGen_polytopesShrinkEvenly(...
    polytopes,...
    des_gap_size/2,...
    (figNum));


% plot the shrunk polytopes
plotFormat.LineWidth = 2;
plotFormat.MarkerSize = 10;
plotFormat.LineStyle = '-';
plotFormat.Color = [1 0 0];

fillFormat = [1 0 0 0 0.5];

% Call the function
h_plot = fcn_MapGen_plotPolytopes(shrunk_polytopes2, (plotFormat),(fillFormat),(figNum));
hold on;

axis(new_axis);
title('Shrunk from edge polytopes');



%% Show typical generation of polytopes from a tiling
% Demos fcn_MapGen_generatePolysFromVoronoiAABB
% and script_test_fcn_MapGen_generatePolysFromVoronoiAABB
figNum = 10;

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
fcn_MapGen_generatePolysFromVoronoiAABB(seed_points,V,C,AABB, stretch,figNum);

%% (DEPRECATED) Show time-dependent generation of polytopes from a tiling
% % Moves just one point
% figNum = 11111;
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
%     figure(figNum);
%     clf;
% 
% 
% % plot the polytopes on current axis
% plotFormat.LineWidth = 2;
% plotFormat.MarkerSize = 10;
% plotFormat.LineStyle = '-';
% plotFormat.Color = [0 0 0];
% 
% fillFormat = [1 0 0 0 0.5];
% 
% % Call the function
% h_plot = fcn_MapGen_plotPolytopes(shrunk_polytopes, (plotFormat),(fillFormat),(figNum));

%     hold on;
%     highlighted_polytope = shrunk_polytopes(pointID_to_move);
% 
%     % Highlight the moving polytope
% % plot the polytopes on current axis
% plotFormat.LineWidth = 2;
% plotFormat.MarkerSize = 10;
% plotFormat.LineStyle = '-';
% plotFormat.Color = [1 0 0];
% 
% fillFormat = [1 1 0 0 1];
% 
% % Call the function
% h_plot = fcn_MapGen_plotPolytopes(highlighted_polytope, (plotFormat),(fillFormat),(figNum));
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
% figNum = 1111;
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
%     figure(figNum);
%     clf;
% 
% 
% % plot the shrunk polytopes
% plotFormat.LineWidth = 2;
% plotFormat.MarkerSize = 10;
% plotFormat.LineStyle = '-';
% plotFormat.Color = [1 0 0];
% 
% fillFormat = [1 0 0 0 0.5];
% 
% % Call the function
% h_plot = fcn_MapGen_plotPolytopes(shrunk_polytopes2, (plotFormat),(fillFormat),(figNum));

%     
%     %     % Highlight the moving polytope

% % plot the shrunk polytopes
% plotFormat.LineWidth = 2;
% plotFormat.MarkerSize = 10;
% plotFormat.LineStyle = '-';
% plotFormat.Color = [1 0 0];
% 
% fillFormat = [1 0 0 0 0.5];
% 
% % Call the function
% h_plot = fcn_MapGen_plotPolytopes(highlighted_polytope, (plotFormat),(fillFormat),(figNum));

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

figNum = 112;

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
polytopes = fcn_MapGen_generatePolysFromVoronoiAABBWithTiling(original_seed_points,unit_AABB, stretch,figNum);
assert(isstruct(polytopes));



%% Show time-dependent generation of polytopes from a tiling
% Moves just one point
figNum = 11111;

% WARNING: VERY slow
if 1==0
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



        % shrunk_polytopes=...
        %     fcn_MapGen_polytopesShrinkFromEdges(...
        %     polytopes,des_gap_size);

        shrunk_polytopes = ...
            fcn_MapGen_polytopesShrinkEvenly(...
            polytopes,...
            des_gap_size/2,...
            (-1));

        % Create even cost
        for ith_poly = 1:length(polytopes)
            shrunk_polytopes(ith_poly).cost = 0.7;
        end

        % Plot the results
        figure(figNum);
        clf;

        % plot the polytopes
        plotFormat.LineWidth = 2;
        plotFormat.MarkerSize = 10;
        plotFormat.LineStyle = '-';
        plotFormat.Color = [0 0 0];

        fillFormat = [1 0 0 1 0.5];


        % Call the function
        fcn_MapGen_plotPolytopes(shrunk_polytopes, (plotFormat),(fillFormat),(figNum));

        hold on;
        highlighted_polytope = shrunk_polytopes(pointID_to_move);

        % Highlight the moving polytope
        % FILL_INFO: a 1-by-5 vector to specify wether or not there is fill, the
        % color of fill, and the opacity of the fill [Y/N, R, G, B, alpha]
        % fill_info = [1 1 0 0 1];

        % plot the polytopes
        plotFormat.LineWidth = 2;
        plotFormat.MarkerSize = 10;
        plotFormat.LineStyle = '-';
        plotFormat.Color = [1 0 0];

        fillFormat = [1 0 0 0 0.5];

        % Call the function
        fcn_MapGen_plotPolytopes(highlighted_polytope, (plotFormat),(fillFormat),(figNum));

        axis([0 1 0 1]);

        if flag_make_video
            frame = getframe(gcf);
            im{ith_step} = frame2im(frame);
        end



        pause(0.1);
    end

    if flag_make_video
        delay_time = 1/25;
        filename = "testAnimated_onePoly_25Hz.gif"; % Specify the output file name
        fullFileNameWithPath = fullfile(pwd,'Images',filename);
        for idx = 1:Nsteps
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,fullFileNameWithPath,"gif","LoopCount",Inf,"DelayTime",delay_time);
            else
                imwrite(A,map,fullFileNameWithPath,"gif","WriteMode","append","DelayTime",delay_time);
            end
        end
    end
end

%% Show time-dependent generation of polytopes from a tiling
% Moves all points
% WARNING: very slow

if 1==0
    figNum = 11111;
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
    velocities = 0.002*randn(Npoints,2);
    % velocities = 0.01;


    flag_make_video = 1;
    % Prep the movie?

    des_gap_size = 0.02;

    Nsteps = 500;
    im{Nsteps} = [];
    for ith_step = 1:Nsteps
        fprintf(1,'Iteration: %.0d\n',ith_step);
        % moved_seed_points(pointID_to_move) = original_seed_points(pointID_to_move) + velocities*ith_step;
        moved_seed_points = original_seed_points + velocities*ith_step;

        % Wrap points
        seed_points = mod(moved_seed_points,1);

        % fill polytopes from tiling
        polytopes = fcn_MapGen_generatePolysFromVoronoiAABBWithTiling(seed_points,unit_AABB, stretch);

        % shrunk_polytopes=...
        %     fcn_MapGen_polytopesShrinkFromEdges(...
        %     polytopes,des_gap_size);

        shrunk_polytopes = ...
            fcn_MapGen_polytopesShrinkEvenly(...
            polytopes,...
            des_gap_size/2,...
            (figNum));


        % Create even cost
        for ith_poly = 1:length(polytopes)
            shrunk_polytopes(ith_poly).cost = 0.7;
        end

        % Plot the results
        figure(figNum);
        clf;


        % plot the polytopes
        plotFormat.LineWidth = 2;
        plotFormat.MarkerSize = 10;
        plotFormat.LineStyle = '-';
        plotFormat.Color = [0 0 0];

        fillFormat = [1 0 0 1 0.5];

        % Call the function
        fcn_MapGen_plotPolytopes(shrunk_polytopes, (plotFormat),(fillFormat),(figNum));

        hold on;
        highlighted_polytope = shrunk_polytopes(pointID_to_move);

        % Highlight the moving polytope
        % FILL_INFO: a 1-by-5 vector to specify wether or not there is fill, the
        % color of fill, and the opacity of the fill [Y/N, R, G, B, alpha]
        % fill_info = [1 1 0 0 1];

        % plot the polytopes
        plotFormat.LineWidth = 2;
        plotFormat.MarkerSize = 10;
        plotFormat.LineStyle = '-';
        plotFormat.Color = [1 0 0];

        fillFormat = [1 0 0 0 0.5];

        % Call the function
        fcn_MapGen_plotPolytopes(highlighted_polytope, (plotFormat),(fillFormat),(figNum));

        axis([0 1 0 1]);

        if flag_make_video
            frame = getframe(gcf);
            im{ith_step} = frame2im(frame);
        end



        pause(0.1);
    end

    if flag_make_video
        delay_time = 1/25;
        filename = "testAnimated_allPoly_25Hz.gif"; % Specify the output file name
        fullFileNameWithPath = fullfile(pwd,'Images',filename);
        for idx = 1:Nsteps
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,fullFileNameWithPath,"gif","LoopCount",Inf,"DelayTime",delay_time);
            else
                imwrite(A,map,fullFileNameWithPath,"gif","WriteMode","append","DelayTime",delay_time);
            end
        end
    end
end


%% Generate a set of polytopes from various pseudo-random sources
close all;

rng(1);

% Generate a set of polytopes from the Sobol set
figNum = 333;
figure(figNum);
clf;
hold on;
Numpoints = 100;

subplot(2,3,1);
% Sobol_range = [1 Numpoints]; % range of Sobol points to use to generate the tiling
% tiled_polytopes = fcn_MapGen_sobolVoronoiTiling(Sobol_range,[1 1],figNum);

seedGeneratorNames = 'sobolset';
seedGeneratorRanges = [1 Numpoints];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[~] = fcn_MapGen_generatePolysFromSeedGeneratorNames(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (figNum)); 

title('Sobol set');
legend('off');

% Generate a set of polytopes from the Halton set
subplot(2,3,2);
seedGeneratorNames = 'haltonset';
seedGeneratorRanges = [1 Numpoints];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[~] = fcn_MapGen_generatePolysFromSeedGeneratorNames(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (figNum)); 


% % plot the polytopes
% plotFormat.LineWidth = 2;
% plotFormat.MarkerSize = 10;
% plotFormat.LineStyle = '-';
% plotFormat.Color = [0 0 1];
%
% fillFormat = [1 0 0 0 0.5];
% 
% % Call the function
% h_plot = fcn_MapGen_plotPolytopes(tiled_polytopes, (plotFormat),(fillFormat),(figNum));
% 

title('Halton set');
legend('off');

% Generate a set of polytopes from the Latin Hypercube set
subplot(2,3,3);
seedGeneratorNames = 'latin';
seedGeneratorRanges = [1 Numpoints];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[~] = fcn_MapGen_generatePolysFromSeedGeneratorNames(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (figNum)); 

title('Latin Hypercube set');
legend('off');

% Generate a set of polytopes from the Random set
subplot(2,3,4);
seedGeneratorNames = 'randn';
seedGeneratorRanges = [1 Numpoints];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[~] = fcn_MapGen_generatePolysFromSeedGeneratorNames(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (figNum)); 

title('Uniform random set');
legend('off');

% Generate a set of polytopes from the Random Normal set
subplot(2,3,5);
seedGeneratorNames = 'rand';
seedGeneratorRanges = [1 Numpoints];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[~] = fcn_MapGen_generatePolysFromSeedGeneratorNames(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (figNum)); 

title('Random normally distributed set');
legend('off');

%% fcn_MapGen_generatePolysFromSeedGeneratorNames - example of multiple tilings on same map
% This is case 10001 in the test script
figNum = 10001;
titleString = sprintf('fcn_MapGen_generatePolysFromSeedGeneratorNames - example of multiple tilings on same map');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;


mapStretch = [1 1];
set_range = [1 100];

rng(1234);

Nsets = 3;
ith_set = 0;
seedGeneratorNames  = cell(Nsets,1);
seedGeneratorRanges = cell(Nsets,1);
AABBs               = cell(Nsets,1);
mapStretchs        = cell(Nsets,1);

ith_set = ith_set+1;
seedGeneratorNames{ith_set,1} = 'haltonset';
seedGeneratorRanges{ith_set,1} = set_range;
AABBs{ith_set,1} = [0 0 1 1];
mapStretchs{ith_set,1} = mapStretch;

ith_set = ith_set+1;
seedGeneratorNames{ith_set} = 'randn';
seedGeneratorRanges{ith_set,1} = set_range;
AABBs{ith_set,1} = [1 0 2 1];
mapStretchs{ith_set,1} = mapStretch;

ith_set = ith_set+1;
seedGeneratorNames{ith_set} = 'randn';
seedGeneratorRanges{ith_set,1} = set_range;
AABBs{ith_set,1} = [1 1 2 2];
mapStretchs{ith_set,1} = mapStretch;

[polytopes] = fcn_MapGen_generatePolysFromSeedGeneratorNames(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(polytopes));
assert(isfield(polytopes,'vertices'));
assert(isfield(polytopes,'xv'));
assert(isfield(polytopes,'yv'));
assert(isfield(polytopes,'distances'));
assert(isfield(polytopes,'mean'));
assert(isfield(polytopes,'area'));
assert(isfield(polytopes,'max_radius'));
assert(isfield(polytopes,'min_radius'));
assert(isfield(polytopes,'mean_radius'));
assert(isfield(polytopes,'radii'));
assert(isfield(polytopes,'cost'));
assert(isfield(polytopes,'parent_poly_id'));

% Check variable sizes
NinSet = set_range(2)-set_range(1)+1;
assert(length(polytopes)== NinSet * Nsets); 

% Check variable values
assert(size(polytopes(1).vertices,2) == 2);
assert(size(polytopes(1).xv,2) >= 2);
assert(size(polytopes(1).yv,2) >= 2);
assert(size(polytopes(1).distances,2) == 1);
assert(isequal(size(polytopes(1).mean), [1 2]));
assert(isequal(size(polytopes(1).area), [1 1]));
assert(isequal(size(polytopes(1).max_radius), [1 1]));
assert(isequal(size(polytopes(1).min_radius), [1 1]));
assert(isequal(size(polytopes(1).mean_radius), [1 1]));
assert(size(polytopes(1).radii,2) == 1);
assert(isequal(size(polytopes(1).cost), [1 1]));
assert(isempty(polytopes(1).parent_poly_id));

% Check variable values
% (these change randomly)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%% Generate many test sets of polytopes from the Halton set
for i=1:100:500
    figNum = 21+i;
    figure(figNum); clf;
    
    seedGeneratorNames = 'haltonset';
    seedGeneratorRanges = [i i+100];
    AABBs = [0 0 1 1];
    mapStretchs = [1 1];
    [~] = fcn_MapGen_generatePolysFromSeedGeneratorNames(...
        seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
        seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
        (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
        (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
        (figNum)); 

    % Do statistics, checking that the area is always fully filled and we
    % get 101 polytopes each time

    % temp = fcn_MapGen_statsPolytopes(...
    % polytopes);
    % title(sprintf('Halton range is: [%.0d %.0d]',i,i+100));
    % assert(abs(temp.unoccupancy_ratio)<(1000*eps));
    % assert(isequal(101,temp.point_density));
    % pause(0.1);
end


%% fcn_MapGen_polytopesDeleteByAABB - allows maps to be trimmed to a box
% This is case 10001
figNum = 10001;
titleString = sprintf('fcn_MapGen_polytopesDeleteByAABB - allows maps to be trimmed to a box');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

seedGeneratorNames = 'haltonset';
seedGeneratorRanges = [1 1000];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[polytopes] = fcn_MapGen_generatePolysFromSeedGeneratorNames(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (-1));

boundingBox = [0.25,0.25, 0.75,0.75];

% Call the function
trimmedPolytopes = fcn_MapGen_polytopesDeleteByAABB(polytopes, boundingBox,(figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(trimmedPolytopes));
assert(isfield(trimmedPolytopes,'vertices'));
assert(isfield(trimmedPolytopes,'xv'));
assert(isfield(trimmedPolytopes,'yv'));
assert(isfield(trimmedPolytopes,'distances'));
assert(isfield(trimmedPolytopes,'mean'));
assert(isfield(trimmedPolytopes,'area'));
assert(isfield(trimmedPolytopes,'max_radius'));
assert(isfield(trimmedPolytopes,'min_radius'));
assert(isfield(trimmedPolytopes,'mean_radius'));
assert(isfield(trimmedPolytopes,'radii'));
assert(isfield(trimmedPolytopes,'cost'));
assert(isfield(trimmedPolytopes,'parent_poly_id'));

% Check variable sizes
assert(length(trimmedPolytopes)<=length(polytopes)); 

% Check variable values
% (can't - randomly generated!)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% Show how the polytopes can be shrunk to a specified radius
% Shrink to radius
figNum = 24;
des_rad = 0.03; sigma_radius = 0; min_rad = 0.001;
shrunk_polytopes2 = fcn_MapGen_polytopesShrinkToRadius(...
    trimmed_polytopes,des_rad,sigma_radius,min_rad,figNum); %#ok<NASGU>

%% Show how different shrinking methods change statistics
figNum = 555;
figure(figNum); clf;

% Generate polytopes from the Halton set        
seedGeneratorNames = 'haltonset';
seedGeneratorRanges = [5401 5501];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[tiled_polytopes] = fcn_MapGen_generatePolysFromSeedGeneratorNames(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (figNum));

    
% Grab statistics on original map
figure(figNum+10);
clf;
fcn_MapGen_statsPolytopes(...
    tiled_polytopes,...
    figNum+10);

% Shrink to radius+10
figNum = 556;
figure(figNum); clf;

des_rad = 0.03; sigma_radius = 0; min_rad = 0.001;
shrunk_polytopes2=fcn_MapGen_polytopesShrinkToRadius(...
    trimmed_polytopes,des_rad,sigma_radius,min_rad,figNum);


%% fcn_MapGen_polytopeShrinkToRadius - shrinks one polytope
% see script_test_fcn_MapGen_polytopeShrinkToRadius
% This is case 10002 from script
titleString = sprintf('fcn_MapGen_polytopeShrinkToRadius - shrinks one polytope');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

seedGeneratorNames = 'haltonset';
seedGeneratorRanges = [1 100];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[polytopes] = fcn_MapGen_generatePolysFromSeedGeneratorNames(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (-1));


bounding_box = [0,0, 1,1];
trim_polytopes = fcn_MapGen_polytopesDeleteByAABB(polytopes,bounding_box,-1);

% Pick a random polytope
Npolys = length(trim_polytopes);
rand_poly = 1+floor(rand*Npolys);
shrinker = trim_polytopes(rand_poly);

orig_radius = shrinker.max_radius;
ratios = (0.99:-0.05:0);

for ith_ratio = 1:length(ratios)
    newRadius = orig_radius*ratios(ith_ratio);

    % Call the function
    shrunkPolytope = fcn_MapGen_polytopeShrinkToRadius(shrinker, newRadius, (figNum));

    sgtitle(titleString, 'Interpreter','none');

    % Check variable types
    assert(isstruct(shrunkPolytope));
    assert(isfield(shrunkPolytope,'vertices'));
    assert(isfield(shrunkPolytope,'xv'));
    assert(isfield(shrunkPolytope,'yv'));
    assert(isfield(shrunkPolytope,'distances'));
    assert(isfield(shrunkPolytope,'mean'));
    assert(isfield(shrunkPolytope,'area'));
    assert(isfield(shrunkPolytope,'max_radius'));
    assert(isfield(shrunkPolytope,'min_radius'));
    assert(isfield(shrunkPolytope,'mean_radius'));
    assert(isfield(shrunkPolytope,'radii'));
    assert(isfield(shrunkPolytope,'cost'));
    assert(isfield(shrunkPolytope,'parent_poly_id'));

    % Check variable sizes
    assert(isequal(length(shrinker),length(shrunkPolytope)));

    % Check variable values
    assert(isequal(round(shrunkPolytope.max_radius,4),round(newRadius,4)));

    pause(0.01);
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% fcn_MapGen_polytopeRemoveTightVerticies - removes tight verticies
figNum = 32;
titleString = sprintf('fcn_MapGen_polytopeRemoveTightVerticies - removes tight verticies');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

orig_radius = shrinker.max_radius;
ratios = (0.99:-0.05:0);

for ith_ratio = 1:length(ratios)
    des_rad = orig_radius*ratios(ith_ratio);
    tolerance = 0.02;
    shrunk_polytope =...
        fcn_MapGen_polytopeShrinkToRadius(...
        shrinker,des_rad, (figNum));
    cleaned_polytope = fcn_MapGen_polytopeRemoveTightVerticies(...
        shrunk_polytope, tolerance,figNum);
    pause(0.01);
end

%% fcn_MapGen_polytopeMapGen - Generate polytope maps in just one function call
% See script_test_fcn_MapGen_polytopeMapGen - this is case 10001 from that
% script

figNum = 10004;
titleString = sprintf('fcn_MapGen_polytopeMapGen - Generate polytope maps in just one function call');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% generate Voronoi tiling from Halton points
haltonRange = [1 200]; % range of Halton points to use to generate the tiling
% remove the edge polytope that extend past the high and low points
xlow = 0; xhigh = 1; ylow = 0; yhigh = 1;
boundingBox = [xlow ylow, xhigh yhigh];

% shink the polytopes so that they are no longer tiled
des_radius = 0.03; % desired average maximum radius
sigma_radius = 0.002; % desired standard deviation in maximum radii
min_rad = 0.0001; % minimum possible maximum radius for any obstacle
shrink_seed = 1111; % seed used for randomizing the shrinking process

% Call the function
[map_polytopes,all_pts,mu_rad_final,sigma_rad_final] = ...
    fcn_MapGen_polytopeMapGen(...
    haltonRange,boundingBox,...
    des_radius,sigma_radius,min_rad,shrink_seed,(figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(map_polytopes));
assert(isfield(map_polytopes,'vertices'));
assert(isfield(map_polytopes,'xv'));
assert(isfield(map_polytopes,'yv'));
assert(isfield(map_polytopes,'distances'));
assert(isfield(map_polytopes,'mean'));
assert(isfield(map_polytopes,'area'));
assert(isfield(map_polytopes,'max_radius'));
assert(isfield(map_polytopes,'min_radius'));
assert(isfield(map_polytopes,'mean_radius'));
assert(isfield(map_polytopes,'radii'));
assert(isfield(map_polytopes,'cost'));
assert(isfield(map_polytopes,'parent_poly_id'));
assert(isnumeric(all_pts));
assert(isnumeric(mu_rad_final));
assert(isnumeric(sigma_rad_final));

% Check variable sizes
assert(isequal(200,length(map_polytopes))); 
assert(size(all_pts,1)>200);
assert(size(all_pts,2)==5);
assert(isequal(size(mu_rad_final),[1 1]));
assert(isequal(size(sigma_rad_final),[1 1]));

% Check variable values
assert(isequal(round(des_radius,3),round(mu_rad_final,3)));
assert(isequal(round(sigma_radius,3),round(sigma_rad_final,3)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% fcn_MapGen_generatePolysFromName - Generate polytope maps from a name
% See script_test_fcn_MapGen_generatePolysFromName - this is case 10001 from that
% script

figNum = 10004;
titleString = sprintf('fcn_MapGen_generatePolysFromName - Generate polytope maps from a name');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;
map_name = "HST 1 100 SQT 0 1 0 1 SMV 0.01 0.001 1e-6 1111";
plot_flag = 1; 
disp_name = 0; 

line_style = 'r-';
line_width = 2;

% Call the function
[polytopes, h_fig] = fcn_MapGen_generatePolysFromName(map_name, plot_flag, disp_name, (figNum), (line_style), (line_width));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(polytopes));
assert(isfield(polytopes,'vertices'));
assert(isfield(polytopes,'xv'));
assert(isfield(polytopes,'yv'));
assert(isfield(polytopes,'distances'));
assert(isfield(polytopes,'mean'));
assert(isfield(polytopes,'area'));
assert(isfield(polytopes,'max_radius'));
assert(isfield(polytopes,'min_radius'));
assert(isfield(polytopes,'mean_radius'));
assert(isfield(polytopes,'radii'));
assert(isfield(polytopes,'cost'));
assert(isfield(polytopes,'parent_poly_id'));
assert(ishandle(h_fig));

% Check variable sizes
Npolys = 100;
assert(isequal(Npolys,length(polytopes))); 
assert(isequal(size(h_fig),[1 1]));

% Check variable values
assert(isequal(h_fig.Number,figNum));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%% fcn_MapGen_polytopesExpandEvenly - expand one polytope
% See script_test_fcn_MapGen_polytopesExpandEvenly - this is case 10001 from that
% script

figNum = 10001;
titleString = sprintf('fcn_MapGen_polytopesExpandEvenly - expand one polytope');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

polytopes = fcn_MapGen_polytopeFillEmptyPoly(-1);

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
expansionDistance = 0.04; % Set the expansion distance

% Call the function
expandedPolytopes = fcn_MapGen_polytopesExpandEvenly(polytopes, expansionDistance, figNum);

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(expandedPolytopes));
assert(isfield(expandedPolytopes,'vertices'));
assert(isfield(expandedPolytopes,'xv'));
assert(isfield(expandedPolytopes,'yv'));
assert(isfield(expandedPolytopes,'distances'));
assert(isfield(expandedPolytopes,'mean'));
assert(isfield(expandedPolytopes,'area'));
assert(isfield(expandedPolytopes,'max_radius'));
assert(isfield(expandedPolytopes,'min_radius'));
assert(isfield(expandedPolytopes,'mean_radius'));
assert(isfield(expandedPolytopes,'radii'));
assert(isfield(expandedPolytopes,'cost'));
assert(isfield(expandedPolytopes,'parent_poly_id'));

% Check variable sizes
assert(isequal(length(polytopes),length(expandedPolytopes))); 

% Check variable values
assert(isequal(round(expandedPolytopes.area,4),0.0150));
assert(isequal(round(expandedPolytopes.max_radius,4),0.1445));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%% fcn_MapGen_polytopesExpandEvenly - expand many polytope
% See script_test_fcn_MapGen_polytopesExpandEvenly - this is case 10002 from that
% script

figNum = 10002;
titleString = sprintf('DEMO case: expansion of a polytope field');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

map_name = "HST 30 450 SQT 0 1 0 1 SMV 0.02 0.005 1e-6 1234";
plot_flag = 1; disp_name = [1, 0.05 -0.05, 0.5 0.5 0.5, 10];
line_style = '-'; line_width = 2; color = [0 0 1];
axis_limits = [0 1 -0.1 1]; axis_style = 'square';
fill_info = [1 1 0 1 0.5];
figNum = 7;

[polytopes,~] =fcn_MapGen_generatePolysFromName(...
    map_name,...
    plot_flag,...
    disp_name,...
    -1,...
    line_style,...
    line_width,....
    color,...
    axis_limits,...
    axis_style,...
    fill_info);

% Set expansion parameters
expansionDistance = 0.01; % Set the expansion distance

% Call the function
expandedPolytopes = fcn_MapGen_polytopesExpandEvenly(polytopes, expansionDistance, figNum);

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(expandedPolytopes));
assert(isfield(expandedPolytopes,'vertices'));
assert(isfield(expandedPolytopes,'xv'));
assert(isfield(expandedPolytopes,'yv'));
assert(isfield(expandedPolytopes,'distances'));
assert(isfield(expandedPolytopes,'mean'));
assert(isfield(expandedPolytopes,'area'));
assert(isfield(expandedPolytopes,'max_radius'));
assert(isfield(expandedPolytopes,'min_radius'));
assert(isfield(expandedPolytopes,'mean_radius'));
assert(isfield(expandedPolytopes,'radii'));
assert(isfield(expandedPolytopes,'cost'));
assert(isfield(expandedPolytopes,'parent_poly_id'));

% Check variable sizes
assert(isequal(length(polytopes),length(expandedPolytopes))); 

% Check variable values
% Too many values

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% fcn_MapGen_statsPolytopes - find polytope field statistics
% See script_test_fcn_MapGen_statsPolytopes - this is case 10001 from that
% script
figNum = 10001;
titleString = sprintf('fcn_MapGen_statsPolytopes - find polytope field statistics');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

range = [200 220];
seedGeneratorNames = 'haltonset';
seedGeneratorRanges = range;
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[polytopes] = fcn_MapGen_generatePolysFromSeedGeneratorNames(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (-1));

% Call the function
polyMapStats = fcn_MapGen_statsPolytopes(polytopes, (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(polyMapStats));

% Check variable sizes
% Too many

% Check variable values
% Too many

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%%%%% Remove close points and redo statistics
tolerance = 0.01;
cleaned_polytopes = polytopes;
for ith_poly = 1:length(polytopes)
    cleaned_polytopes(ith_poly) = fcn_MapGen_polytopeRemoveTightVerticies(...
        polytopes(ith_poly), tolerance);    
end

figNum = 20;
figure(figNum);clf;

fcn_MapGen_statsPolytopes(...
    cleaned_polytopes,...
    figNum);




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
