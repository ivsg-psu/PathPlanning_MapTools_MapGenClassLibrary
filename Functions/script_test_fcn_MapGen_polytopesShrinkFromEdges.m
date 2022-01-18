% script_fcn_MapGen_polytopeShrinkFromEdges
% Tests function: fcn_MapGen_polytopeShrinkFromEdges

% REVISION HISTORY:
% 2022_01_17
% -- first written by S. Harnett using
% script_test_fcn_MapGen_polytopesShrinkToRadius as a template

%% Set up variables
close all;
polytopes = fcn_MapGen_haltonVoronoiTiling([1 100]);

fig_num = 1;
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);

%% Basic example of uniform shrinking
fig_num = 2;
des_rad = 0.05; sigma_radius = 0; min_rad = 0.001;
shrunk_polytopes1=...
    fcn_MapGen_polytopesShrinkFromEdges(...
    trim_polytopes,des_rad,sigma_radius,min_rad,fig_num);

assert(isequal(round(shrunk_polytopes1(1).max_radius,4),round(des_rad,4)));

%% Basic example of non-uniform shrinking
fig_num = 3;
des_rad = 0.05; sigma_radius = 0.01; min_rad = 0.001;
shrunk_polytopes2=fcn_MapGen_polytopesShrinkFromEdges(...
    trim_polytopes,des_rad,sigma_radius,min_rad,fig_num);


%% Warning thrown because of truncation
% This happens, for example, where there is a large standard deviation with small radius
fig_num = 4;
des_rad = 0.001; sigma_radius = 0.01; min_rad = 0.0001;
shrunk_polytopes1=...
    fcn_MapGen_polytopesShrinkFromEdges(...
    trim_polytopes,des_rad,sigma_radius,min_rad,fig_num);


%% Error cases follow

if 1==0
    %% Error thrown because minimum radius less than zero
    fig_num = 4;
    des_rad = 0.05; sigma_radius = 0.01; min_rad = -0.001;
    shrunk_polytopes1=...
        fcn_MapGen_polytopesShrinkFromEdges(...
        trim_polytopes,des_rad,sigma_radius,min_rad,fig_num); %#ok<*NASGU>
end
