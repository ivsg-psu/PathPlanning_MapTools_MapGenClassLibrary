% script_test_fcn_MapGen_polytopesShrinkToRadius
% Tests function: fcn_MapGen_polytopesShrinkToRadius

% REVISION HISTORY:
% 2021_06_09
% -- first written by S. Brennan using
% script_test_fcn_MapGen_polytopeCropEdges as a template

close all;

%% Basic example of uniform shrinking
fig_num = 2;
figure(fig_num);
clf;

% Set up variables
polytopes = fcn_MapGen_haltonVoronoiTiling([1 100]);


% Fill in test data
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,-1);


des_rad = 0.05; sigma_radius = 0; min_rad = 0.001;
shrunk_polytopes=...
    fcn_MapGen_polytopesShrinkToRadius(...
    trim_polytopes,des_rad,sigma_radius,min_rad,fig_num);
field_stats = fcn_MapGen_polytopesStatistics(shrunk_polytopes);

assert(isstruct(shrunk_polytopes));
assert(isequal(round(field_stats.average_max_radius,4),round(des_rad,4)));

%% Basic example of non-uniform shrinking
fig_num = 3;
figure(fig_num);
clf;

% Set up variables
polytopes = fcn_MapGen_haltonVoronoiTiling([1 100]);


% Fill in test data
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,-1);

des_rad = 0.05; sigma_radius = 0.01; min_rad = 0.001;
shrunk_polytopes = fcn_MapGen_polytopesShrinkToRadius(...
    trim_polytopes,des_rad,sigma_radius,min_rad,fig_num);
assert(isstruct(shrunk_polytopes));

%% Warning thrown because of truncation
% This happens, for example, where there is a large standard deviation with small radius
fig_num = 4;
figure(fig_num);
clf;

% Set up variables
polytopes = fcn_MapGen_haltonVoronoiTiling([1 100]);


% Fill in test data
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,-1);

des_rad = 0.001; sigma_radius = 0.01; min_rad = 0.0001;
shrunk_polytopes=...
    fcn_MapGen_polytopesShrinkToRadius(...
    trim_polytopes,des_rad,sigma_radius,min_rad,fig_num);
assert(isstruct(shrunk_polytopes));

%% Error cases follow

if 1==0
    %% Error thrown because minimum radius less than zero
    fig_num = 4;
    des_rad = 0.05; sigma_radius = 0.01; min_rad = -0.001;
    shrunk_polytopes=...
        fcn_MapGen_polytopesShrinkToRadius(...
        trim_polytopes,des_rad,sigma_radius,min_rad,fig_num); %#ok<*NASGU>
end
