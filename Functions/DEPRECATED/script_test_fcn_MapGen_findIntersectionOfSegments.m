% script_test_fcn_MapGen_findIntersectionOfSegments
% This is a script to exercise the function: fcn_MapGen_findIntersectionOfSegments.m
% This function was written on 2021_06_05 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Modification history:
%      2021_06_05
%      -- wrote function, adapted from script_test_fcn_MapGen_findIntersectionOfSegments.m
% 2025_04_26 - S. Brennan
% -- Started adding assertions, better print statements. Made it about 1/3
% of way. Need to finish

close all

%% Simple test 1 - a simple intersection
fig_num = 1;
figure(fig_num);
clf;

wall_start = [0 10];
wall_end   = [10 10];
sensor_vector_start = [2 1];
sensor_vector_end   = [5 15];
flag_search_type = 0;

[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
title('Simple intersection result');

assert(isequal(size(distance),[1 1]));
assert(isequal(size(location),[1 2]));
assert(isequal(round(distance,4),9.2043));
assert(isequal(round(location,4),[3.9286,10]));

%% Simple test 2 - no intersections
fig_num = 2;
figure(fig_num);
clf;

wall_start = [-4 10];
wall_end   = [2 10];
sensor_vector_start = [0 0];
sensor_vector_end   = [5 12];
flag_search_type = 0;

[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
title('No intersection result');

assert(isequal(size(distance),[1 1]));
assert(isequal(size(location),[1 2]));
%assert(isequal(round(distance,4),9.2043));
assert(all(isnan(distance)));
%assert(isequal(round(location,4),[3.9286,10]));
assert(all(isnan(location)));

%% Simple test 3 - multiple intersections
fig_num = 3;
figure(fig_num);
clf;

wall_start = [0 10; 10 10; 0 6; 10 6];
wall_end = [10 10; 0 6; 10 6; 0 2];

sensor_vector_start = [0 0];
sensor_vector_end   = [5 12];
flag_search_type = 0;

[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
title('Multiple intersections result');


assert(isequal(size(distance),[1 1]));
assert(isequal(size(location),[1 2]));
assert(isequal(round(distance,4),2.6000));
% assert(all(isnan(distance)));
assert(isequal(round(location,4),[1.0000    2.4000]));
% assert(all(isnan(location)));

%% Simple test 4 - intersection through a vertex
fig_num = 4;
figure(fig_num);
clf;

wall_start = [0 5; 4 5];
wall_end = [4 5; 8 2];
sensor_vector_start = [4 0];
sensor_vector_end   = [4 8];
flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
title('Intersection through a vertex result');


assert(isequal(size(distance),[1 1]));
assert(isequal(size(location),[1 2]));
assert(isequal(round(distance,4),5.0000));
% assert(all(isnan(distance)));
assert(isequal(round(location,4),[4 5]));
% assert(all(isnan(location)));

%% Simple test 5 - intersection at start of sensor
fig_num = 5;
figure(fig_num);
clf;

path = [0 5; 4 5; 8 2];
wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [4 5];
sensor_vector_end   = [4 8];
flag_search_type = 0;

[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
title('Intersection at start of sensor result');


assert(isequal(size(distance),[1 1]));
assert(isequal(size(location),[1 2]));
assert(isequal(round(distance,4),0.0000));
% assert(all(isnan(distance)));
assert(isequal(round(location,4),[4 5]));
% assert(all(isnan(location)));

%% Simple test 6 - intersection at end of sensor
fig_num = 6;
figure(fig_num);
clf;

path = [0 5; 4 5; 8 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [4 0];
sensor_vector_end   = [4 5];
flag_search_type = 0;

[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
title('Intersection at end of sensor result');


assert(isequal(size(distance),[1 1]));
assert(isequal(size(location),[1 2]));
assert(isequal(round(distance,4),5.0000));
% assert(all(isnan(distance)));
assert(isequal(round(location,4),[4 5]));
% assert(all(isnan(location)));



%% Simple test - identically overlapping colinear
fig_num = 6;
figure(fig_num);
clf;
title_string = 'Simple test - identically overlapping colinear';
fprintf(1,'%s: \n',title_string);

path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [0 10];
sensor_vector_end   = [10 10];

flag_search_type = 0;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);
print_more_results(distance,location,path_segments);

title(sprintf(title_string));


assert(isequal(size(distance),[1 1]));
assert(isequal(size(location),[1 2]));
assert(isequal(round(distance,4),0.0000));
% assert(all(isnan(distance)));
assert(isequal(round(location,4),[0 10]));
% assert(all(isnan(location)));


%% Simple test - partially overlapping colinear 2
fig_num = 8;
figure(fig_num);
clf;
title_string = 'Simple test - partially overlapping colinear 2';
fprintf(1,'%s: \n',title_string);

path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-2 10];
sensor_vector_end   = [10 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Simple test - partially overlapping colinear 3
fig_num = 9;
figure(fig_num);
clf;
title_string = 'Simple test - partially overlapping colinear 3';
fprintf(1,'%s: \n',title_string);


path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-2 10];
sensor_vector_end   = [5 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Simple test - partially overlapping colinear 4
fig_num = 10;
figure(fig_num);
clf;
title_string = 'Simple test - partially overlapping colinear 4';
fprintf(1,'%s: \n',title_string);

path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [3 10];
sensor_vector_end   = [5 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Simple test - partially overlapping colinear 5
fig_num = 11;
figure(fig_num);
clf;
title_string = 'Simple test - partially overlapping colinear 5';
fprintf(1,'%s: \n',title_string);

path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [3 10];
sensor_vector_end   = [15 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Simple test - super overlapping colinear 1
fig_num = 12;
figure(fig_num);
clf;
title_string = 'Simple test - super overlapping colinear 1';
fprintf(1,'%s: \n',title_string);

path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [15 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Simple test - end overlapping colinear 1
fig_num = 13;
figure(fig_num);
clf;
title_string = 'Simple test - end overlapping colinear 1';
fprintf(1,'%s: \n',title_string);

path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [0 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Simple test 14 - end overlapping colinear 2
fig_num = 14;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [10 10];
sensor_vector_end   = [15 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);


%% Simple test 27 - identically overlapping colinear BACKWARDS
fig_num = 27;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [0 10];
sensor_vector_start   = [10 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Simple test - partially overlapping colinear 1 BACKWARDS
fig_num = 28;
figure(fig_num);
clf;
title_string = 'Simple test - partially overlapping colinear 1 BACKWARDS';
fprintf(1,'%s: \n',title_string);

path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-2 10];
sensor_vector_start   = [10 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Simple test 29 - partially overlapping colinear 1 BACKWARDS
fig_num = 29;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);


fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-2 10];
sensor_vector_start   = [5 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Simple test 30 - partially overlapping colinear 1 BACKWARDS
fig_num = 30;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);



fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [3 10];
sensor_vector_start   = [5 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Simple test 31 - partially overlapping colinear 1 BACKWARDS
fig_num = 31;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);



fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [3 10];
sensor_vector_start   = [15 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Simple test 32 - super overlapping colinear 1 BACKWARDS
fig_num = 32;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);



fprintf(1,'Super overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-3 10];
sensor_vector_start   = [15 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);


%% Simple test 33 - end overlapping colinear 1 BACKWARDS
fig_num = 33;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);



fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-3 10];
sensor_vector_start   = [0 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Simple test 34 - end overlapping colinear 2 BACKWARDS
fig_num = 34;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);



fprintf(1,'End overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [10 10];
sensor_vector_start   = [15 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);



%% Simple test 35 - non overlapping colinear 1
fig_num = 35;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);



fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [-1 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Simple test 36 - non overlapping colinear 2
fig_num = 36;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);



fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [13 10];
sensor_vector_end   = [15 10];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);


%% Advanced test 1 - intersection beyond a sensor's range with flag
fig_num = 901;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);



fprintf(1,'Intersection beyond sensor range result: \n');
path = [0 5; 4 5; 8 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [4 0];
sensor_vector_end   = [4 2];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%%
fig_num = 902;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

path = [0 5; 4 5; 8 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [4 0];
sensor_vector_end   = [4 2];
flag_search_type = 1;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%%
fig_num = 903;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);


% Test the negative condition
path = [0 5; 4 5; 8 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [4 6];
sensor_vector_end   = [4 8];

flag_search_type = 0;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%%
fig_num = 904;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

path = [0 5; 4 5; 8 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [4 6];
sensor_vector_end   = [4 8];
fig_debugging = 2346;
flag_search_type = 1;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%   __  __       _ _   _ _    _ _ _
%  |  \/  |     | | | (_) |  | (_) |
%  | \  / |_   _| | |_ _| |__| |_| |_
%  | |\/| | | | | | __| |  __  | | __|
%  | |  | | |_| | | |_| | |  | | | |_
%  |_|  |_|\__,_|_|\__|_|_|  |_|_|\__|
%
%


%% Advanced test 2 - multiple intersections
fig_num = 2001;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Multiple intersections reporting all results: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [0 0];
sensor_vector_end   = [5 12];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced test 3 - multiple intersections possible, but no hits
fig_num = 2002;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Multiple intersections possible but no hits, reporting all results: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [0 0];
sensor_vector_end   = [0.5 1.2];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced test 4 - multiple intersections possible, but few hits
fig_num = 2003;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Multiple intersections possible but few hits, reporting all results: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [0 0];
sensor_vector_end   = [2.5 6];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);



%   __  __       _ _   _ _    _ _ _    ____                 _                   _
%  |  \/  |     | | | (_) |  | (_) |  / __ \               | |                 (_)
%  | \  / |_   _| | |_ _| |__| |_| |_| |  | |_   _____ _ __| | __ _ _ __  _ __  _ _ __   __ _
%  | |\/| | | | | | __| |  __  | | __| |  | \ \ / / _ \ '__| |/ _` | '_ \| '_ \| | '_ \ / _` |
%  | |  | | |_| | | |_| | |  | | | |_| |__| |\ V /  __/ |  | | (_| | |_) | |_) | | | | | (_| |
%  |_|  |_|\__,_|_|\__|_|_|  |_|_|\__|\____/  \_/ \___|_|  |_|\__,_| .__/| .__/|_|_| |_|\__, |
%                                                                  | |   | |             __/ |
%                                                                  |_|   |_|            |___/


%% Advanced Multihit Overlapping test - identically overlapping colinear
fig_num = 2004;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'identically overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [0 10];
sensor_vector_end   = [10 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);
print_more_results(distance,location,path_segments);


%% Advanced Multihit Overlapping  test 8 - partially overlapping colinear 1
fig_num = 2005;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-2 10];
sensor_vector_end   = [10 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 9 - partially overlapping colinear 1
fig_num = 2006;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-2 10];
sensor_vector_end   = [5 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 10 - partially overlapping colinear 1
fig_num = 2007;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [3 10];
sensor_vector_end   = [5 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 11 - partially overlapping colinear 1
fig_num = 2008;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [3 10];
sensor_vector_end   = [15 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 12 - super overlapping colinear 1
fig_num = 2009;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Super overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [15 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 13 - end overlapping colinear 1
fig_num = 2010;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [0 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 14 - end overlapping colinear 2
fig_num = 2014;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [10 10];
sensor_vector_end   = [15 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);






%% Advanced Multihit Overlapping  test 27 - identically overlapping colinear BACKWARDS
fig_num = 2027;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'identically overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [0 10];
sensor_vector_start   = [10 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 28 - partially overlapping colinear 1 BACKWARDS
fig_num = 2028;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [-2 10];
sensor_vector_start   = [10 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 29 - partially overlapping colinear 1 BACKWARDS
fig_num = 2029;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [-2 10];
sensor_vector_start   = [5 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 30 - partially overlapping colinear 1 BACKWARDS
fig_num = 2030;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [3 10];
sensor_vector_start   = [5 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 31 - partially overlapping colinear 1 BACKWARDS
fig_num = 2031;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [3 10];
sensor_vector_start   = [15 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 32 - super overlapping colinear 1 BACKWARDS
fig_num = 2032;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Super overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [-3 10];
sensor_vector_start   = [15 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);


%% Advanced Multihit Overlapping  test 33 - end overlapping colinear 1 BACKWARDS
fig_num = 2033;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [-3 10];
sensor_vector_start   = [0 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 34 - end overlapping colinear 2 BACKWARDS
fig_num = 2034;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'End overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [10 10];
sensor_vector_start   = [15 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);



%% Advanced Multihit Overlapping  test 15 - non overlapping colinear 1
fig_num = 2035;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [-1 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 15 - non overlapping colinear 2
fig_num = 2036;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_start = [13 10];
sensor_vector_end   = [15 10];

flag_search_type = 2;
[distance,location] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_results(distance,location);



%   __  __       _ _   _ _    _ _ _   __  __       _ _ _   _____      _   _
%  |  \/  |     | | | (_) |  | (_) | |  \/  |     | (_) | |  __ \    | | | |
%  | \  / |_   _| | |_ _| |__| |_| |_| \  / |_   _| |_| |_| |__) |_ _| |_| |__
%  | |\/| | | | | | __| |  __  | | __| |\/| | | | | | | __|  ___/ _` | __| '_ \
%  | |  | | |_| | | |_| | |  | | | |_| |  | | |_| | | | |_| |  | (_| | |_| | | |
%  |_|  |_|\__,_|_|\__|_|_|  |_|_|\__|_|  |_|\__,_|_|_|\__|_|   \__,_|\__|_| |_|
%
%


%% Advanced Multihit Overlapping test - identically overlapping colinear
fig_num = 3001;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'identically overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [0 10];
sensor_vector_end   = [10 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 8 - partially overlapping colinear 1
fig_num = 3002;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-2 10];
sensor_vector_end   = [10 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 9 - partially overlapping colinear 1
fig_num = 3003;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-2 10];
sensor_vector_end   = [5 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 10 - partially overlapping colinear 1
fig_num = 3004;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [3 10];
sensor_vector_end   = [5 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 11 - partially overlapping colinear 1
fig_num = 3005;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [3 10];
sensor_vector_end   = [15 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 12 - super overlapping colinear 1
fig_num = 3006;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Super overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [15 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 13 - end overlapping colinear 1
fig_num = 3007;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [0 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 14 - end overlapping colinear 2
fig_num = 3008;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [10 10];
sensor_vector_end   = [15 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);






%% Advanced Multihit Overlapping  test 3009 - identically overlapping colinear BACKWARDS
fig_num = 3009;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'identically overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [0 10];
sensor_vector_start   = [10 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 3010 - partially overlapping colinear 1 BACKWARDS
fig_num = 3010;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-2 10];
sensor_vector_start   = [10 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 3011 - partially overlapping colinear 1 BACKWARDS
fig_num = 3011;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-2 10];
sensor_vector_start   = [5 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 3012 - partially overlapping colinear 1 BACKWARDS
fig_num = 3012;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [3 10];
sensor_vector_start   = [5 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 3013 - partially overlapping colinear 1 BACKWARDS
fig_num = 3013;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [3 10];
sensor_vector_start   = [15 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 3014 - super overlapping colinear 1 BACKWARDS
fig_num = 3014;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Super overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-3 10];
sensor_vector_start   = [15 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);


%% Advanced Multihit Overlapping  test 3015 - end overlapping colinear 1 BACKWARDS
fig_num = 3015;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-3 10];
sensor_vector_start   = [0 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 3016 - end overlapping colinear 2 BACKWARDS
fig_num = 3016;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'End overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [10 10];
sensor_vector_start   = [15 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);



%% Advanced Multihit Overlapping  test 3017 - non overlapping colinear 1
fig_num = 3017;
figure(fig_num);
clf;
title_string = 'Partially overlapping colinear result';
fprintf(1,'%s: \n',title_string);

fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [-1 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test - non overlapping colinear 2
fig_num = 3018;
figure(fig_num);
clf;
title_string = 'Advanced Multihit Overlapping  test - non overlapping colinear 2';
fprintf(1,'%s: \n',title_string);

path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [13 10];
sensor_vector_end   = [15 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%% Advanced Random Multihit - 3019
fig_num = 3019;
figure(fig_num);
clf;
title_string = 'Advanced Random Multihit';
fprintf(1,'%s: \n',title_string);


Num_walls = 10;
wall_start = 10*rand(Num_walls,2);
wall_end   = 10*rand(Num_walls,2);
sensor_vector_start = [0 0];
sensor_vector_end   = [10 10];

flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);


%% Advanced Random Single Hit - 3020

fig_num = 3020;
figure(fig_num);
clf;
title_string = 'Advanced Random Single Hit';
fprintf(1,'%s: \n',title_string);

Num_walls = 10;
wall_start = 10*rand(Num_walls,2);
wall_end   = 10*rand(Num_walls,2);
sensor_vector_start = [0 0];
sensor_vector_end   = [10 10];

flag_search_type = 0;
[distance,location,path_segments] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);
print_more_results(distance,location,path_segments);

%%
function print_results(distance,location)
fprintf(1,'Distance \t Location X \t Location Y \n');
if ~isempty(distance)
    for i_result = 1:length(distance(:,1))
        fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f\n',distance(i_result),location(i_result,1),location(i_result,2));
    end
end
end

%%
function print_more_results(distance,location,path_segments)
fprintf(1,'Distance \t Location X \t Location Y \t PathSegment \n');
if ~isempty(distance)
    for i_result = 1:length(distance(:,1))
        fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f \t\t %.0d\n',distance(i_result),location(i_result,1),location(i_result,2),path_segments(i_result));
    end
end
end
