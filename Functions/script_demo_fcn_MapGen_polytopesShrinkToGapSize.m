% script_demo_fcn_MapGen_polytopesShrinkToGapSize
% Repeatedly runs function: fcn_MapGen_polytopesShrinkToGapSize
% to generate data and compare commanded gap size and observed gap size

% REVISION HISTORY:
% 2022_01_17 by Sean Brennan, sbrennan@psu.edu
% - first written by S. Harnett using
%   % script_test_fcn_MapGen_polytopesShrinkToRadius as a template
% 
% 2025_04_28 by Sean Brennan, sbrennan@psu.edu
% - renamed from script_test_* by S. Harnett
% - see script_test_* for single test case
% 
% 2025_11_20 by Sean Brennan, sbrennan@psu.edu
% - Updated rev history to be in Markdown format
% - Replaced fig_+num with figNum

% TO-DO:
% 
% 2025_11_20 by Sean Brennan, sbrennan@psu.edu
% - fill in to-do items here.


close all;


%% Set up variables
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


figNum = 1;
bounding_box = [0,0, 1,1];
trim_polytopes = fcn_MapGen_polytopesDeleteByAABB(polytopes,bounding_box,figNum);
pre_shrink_stats = fcn_MapGen_statsPolytopes(trim_polytopes, -1);
R_bar_initial = pre_shrink_stats.average_max_radius;
des_gap_size_all = [];
R_bar_final_all = [];
rho_all = [];

% loop through gap sizes from rather small to quite large
% The following code is VERY slow - uncomment to see statistics in error
% using the shink-from-edges approach.
for ith_size = linspace(0.001,0.08,5)
    figNum = 2;
    des_gap_size = ith_size;
    shrunk_polytopes1=...
        fcn_MapGen_polytopesShrinkToGapSize(...
        trim_polytopes,des_gap_size, -1);
    field_stats = fcn_MapGen_statsPolytopes(shrunk_polytopes1, -1);
    rho = field_stats.linear_density_mean;
    R_bar_final = field_stats.average_max_radius;
    R_bar_final_all = [R_bar_final_all; R_bar_final];
    des_gap_size_all = [des_gap_size_all; des_gap_size];
    rho_all = [rho_all; rho];
end

figure
hold on
box on
grid on
plot(des_gap_size_all,des_gap_size_all-2*(R_bar_final_all - R_bar_initial))
xlabel('desired or commanded gap size')
ylabel('difference between commanded gap size and gap size estimate from average max radius change')

figure
hold on
box on
grid on
plot(des_gap_size_all./R_bar_initial*100,2*(R_bar_final_all - R_bar_initial))
xlabel('desired or commanded gap size as percent of initial average max radius')
ylabel('gap size estimate from average max radius change')
