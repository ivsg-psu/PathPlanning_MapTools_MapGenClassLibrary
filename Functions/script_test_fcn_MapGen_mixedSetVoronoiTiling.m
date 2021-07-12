% script_test_fcn_MapGen_mixedSetVoronoiTiling
% Tests function: fcn_MapGen_mixedSetVoronoiTiling

% REVISION HISTORY:
% 2021_06_06 
% -- first written by S. Brennan. 

%'haltonset','sobolset','lhsdesign','rand','randn'

close all;
clear mixedSet

fig_num = 1;
stretch = [1 1];
set_range = [1 100];

rng(1234);

mixedSet(1).name = 'haltonset';
mixedSet(1).settings = set_range;
mixedSet(1).AABB = [0 0 1 1];

mixedSet(2).name = 'randn';
mixedSet(2).settings = set_range;
mixedSet(2).AABB = [1 0 2 1];

mixedSet(3).name = 'randn';
mixedSet(3).settings = set_range;
mixedSet(3).AABB = [1 1 2 2];

polytopes = fcn_MapGen_mixedSetVoronoiTiling(mixedSet,stretch,fig_num);



%% Show that the stretch works
fig_num = 200;
stretch = [100 200];
polytopes = fcn_MapGen_mixedSetVoronoiTiling(mixedSet,stretch,fig_num);

%% Animate a set moving sideways
close all;
fig_num = 333;
stretch = [1 1];
set_range = [1 100];

rng(1234);

clear mixedSet
for ith_set = 1:10
    mixedSet(ith_set).name = 'haltonset'; %#ok<*SAGROW>
    mixedSet(ith_set).settings = set_range+[100 100]*(ith_set-1);
    mixedSet(ith_set).AABB = [ith_set-1 0 ith_set 1];
end

polytopes = fcn_MapGen_mixedSetVoronoiTiling(mixedSet,stretch);

% Plot the polytopes
line_width = 2;
axis_limits = [0 1 0 1];
axis_stype = 'equal';
fcn_MapGen_plotPolytopes(polytopes,fig_num,'b',line_width,axis_limits,axis_stype);

step = 0.01;
for ith_step = 0:step:9
    axis([ith_step ith_step+1 0 1]);
    pause(0.02);
end

%% Create an overlapping set
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



