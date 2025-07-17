% script_test_fcn_MapGen_voronoiTiling
% Tests function: fcn_MapGen_voronoiTiling

% REVISION HISTORY:
% 2021_06_06
% -- first written by S. Brennan for fcn_MapGen_mixedSetVoronoiTiling
% 2025_07_11 - S. Brennan, sbrennan@psu.edu
% -- updated script testing to standard form

%% Set up the workspace
close all

%% Code demos start here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                              ____   __    _____          _
%  |  __ \                            / __ \ / _|  / ____|        | |
%  | |  | | ___ _ __ ___   ___  ___  | |  | | |_  | |     ___   __| | ___
%  | |  | |/ _ \ '_ ` _ \ / _ \/ __| | |  | |  _| | |    / _ \ / _` |/ _ \
%  | |__| |  __/ | | | | | (_) \__ \ | |__| | |   | |___| (_) | (_| |  __/
%  |_____/ \___|_| |_| |_|\___/|___/  \____/|_|    \_____\___/ \__,_|\___|
%
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Demos%20Of%20Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 1

close all;
fprintf(1,'Figure: 1XXXXXX: DEMO cases\n');

%% DEMO case: basic example of multiple tilings on same map
fig_num = 10001;
titleString = sprintf('DEMO case: basic example of multiple tilings on same map');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;


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

[polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num));

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
assert(isequal(get(gcf,'Number'),fig_num));

%% DEMO case: Create overlapping sets
fig_num = 10002;
titleString = sprintf('DEMO case: Create overlapping sets');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;


mapStretch = [1 1];
set_range = [1 100];

rng(1234);

Nsets = 2;
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
seedGeneratorNames{ith_set} = 'rand';
seedGeneratorRanges{ith_set,1} = set_range;
AABBs{ith_set,1} = [0.5 0 0.75 1];
mapStretchs{ith_set,1} = mapStretch;

[polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num));

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
assert(isequal(get(gcf,'Number'),fig_num));

%% DEMO case: Animate a set moving sideways
close all;
fig_num = 10003;
titleString = sprintf('DEMO case: Animate a set moving sideways');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;


mapStretch = [1 1];
set_range = [1 100];

rng(1234);

Nsets = 3;
ith_set = 0;
seedGeneratorNames  = cell(Nsets,1);
seedGeneratorRanges = cell(Nsets,1);
AABBs               = cell(Nsets,1);
mapStretchs        = cell(Nsets,1);

% mixedSet(ith_set).name = 'haltonset'; %#ok<*SAGROW>
% mixedSet(ith_set).settings = set_range+[100 100]*(ith_set-1);
% mixedSet(ith_set).AABB = [ith_set-1 0 ith_set 1];
for ith_generator = 1:Nsets
    ith_set = ith_set+1;
    seedGeneratorNames{ith_set,1} = 'haltonset';
    seedGeneratorRanges{ith_set,1} = set_range+[100 100]*(ith_set-1);
    AABBs{ith_set,1} = [ith_set-1 0 ith_set 1];
    mapStretchs{ith_set,1} = mapStretch;
end

[polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num));

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
assert(isequal(get(gcf,'Number'),fig_num));

%%
if 1==1
    % Animate the polytopes
    axis([0 1 0 1]);
    axis equal;
    step = 0.01;
    for ith_step = 0:step:2
        axis([ith_step ith_step+1 0 1]);
        pause(0.02);
    end
end

%% Test cases start here. These are very simple, usually trivial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  _______ ______  _____ _______ _____
% |__   __|  ____|/ ____|__   __/ ____|
%    | |  | |__  | (___    | | | (___
%    | |  |  __|  \___ \   | |  \___ \
%    | |  | |____ ____) |  | |  ____) |
%    |_|  |______|_____/   |_| |_____/
%
%
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=TESTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 2

close all;
fprintf(1,'Figure: 2XXXXXX: TEST mode cases\n');

%% TEST case: Show that the mapStretch works with multiple generators
fig_num = 20001;
titleString = sprintf('TEST case: Show that the mapStretch works with multiple generators');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

mapStretch = [100 200];
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

[polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num));

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
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: Single set test
fig_num = 20001;
titleString = sprintf('TEST case: Single set test');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

mapStretch = [1 1];
set_range = [1 100];

rng(1234);

Nsets = 1;
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

[polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num));

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
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: Single set test with sobolset
fig_num = 20002;
titleString = sprintf('TEST case: Single set test with sobolset');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

mapStretch = [1 1];
set_range = [1 100];

rng(123);

Nsets = 1;
ith_set = 0;
seedGeneratorNames  = cell(Nsets,1);
seedGeneratorRanges = cell(Nsets,1);
AABBs               = cell(Nsets,1);
mapStretchs        = cell(Nsets,1);

ith_set = ith_set+1;
seedGeneratorNames{ith_set,1} = 'sobolset';
seedGeneratorRanges{ith_set,1} = set_range;
AABBs{ith_set,1} = [0 0 1 1];
mapStretchs{ith_set,1} = mapStretch;

[polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num));

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
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: Single set test with non cell inputs
fig_num = 20003;
titleString = sprintf('TEST case: Single set test with non cell inputs');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;



rng(1234);

% seedGeneratorNames  = cell(Nsets,1);
% seedGeneratorRanges = cell(Nsets,1);
% AABBs               = cell(Nsets,1);
% mapStretchs        = cell(Nsets,1);
Nsets = 1;
mapStretch = [1 1];
set_range = [1 100];
seedGeneratorNames = 'haltonset';
seedGeneratorRanges = set_range;
AABBs = [0 0 1 1];
mapStretchs = mapStretch;

[polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (fig_num));

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
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: Test of all generator types
fig_num = 20004;
titleString = sprintf('TEST case: Test of all generator types');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

allGeneratorTypes = {'haltonset', 'sobolset','lhsdesign','rand','randn'};

mapStretch = [1 1];
set_range = [1 21];

% BUG: Produces error! set_range = [1 20];

rng(1234);

Nsets = 1; 

seedGeneratorNames  = cell(Nsets,1);
seedGeneratorRanges = cell(Nsets,1);
AABBs               = cell(Nsets,1);
mapStretchs        = cell(Nsets,1);

for ith_set = 1:length(allGeneratorTypes)

    seedGeneratorNames{1,1} = allGeneratorTypes{ith_set};
    seedGeneratorRanges{1,1} = set_range;
    AABBs{1,1} = [0 0 1 1];
    mapStretchs{1,1} = mapStretch;

    figure(fig_num);
    nexttile;

    [polytopes] = fcn_MapGen_voronoiTiling(...
        seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
        seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
        (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
        (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
        (fig_num));
    title(allGeneratorTypes{ith_set})

    sgtitle(titleString, 'Interpreter','none');
    legend('off');
    axis equal

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
    assert(isequal(get(gcf,'Number'),fig_num));
end


%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 8

close all;
fprintf(1,'Figure: 8XXXXXX: FAST mode cases\n');

%% Basic example - NO FIGURE
fig_num = 80001;
fprintf(1,'Figure: %.0f: FAST mode, empty fig_num\n',fig_num);
figure(fig_num); close(fig_num);

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

[polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    ([]));

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

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

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

[polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (-1));

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


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

mapStretch = [1 1];
set_range = [1 100];

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

Niterations = 5;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations

    rng(1234);

    % Call the function
    [polytopes] = fcn_MapGen_voronoiTiling(...
        seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
        seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
        (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
        (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
        ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    rng(1234);

    % Call the function
    [polytopes] = fcn_MapGen_voronoiTiling(...
        seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
        seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
        (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
        (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
        (-1));
end
fast_method = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

% Plot results as bar chart
figure(373737);
clf;
hold on;

X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% BUG cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____  _    _  _____
% |  _ \| |  | |/ ____|
% | |_) | |  | | |  __    ___ __ _ ___  ___  ___
% |  _ <| |  | | | |_ |  / __/ _` / __|/ _ \/ __|
% | |_) | |__| | |__| | | (_| (_| \__ \  __/\__ \
% |____/ \____/ \_____|  \___\__,_|___/\___||___/
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=BUG%20cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All bug case figures start with the number 9

% close all;

%% BUG

%% Fail conditions
if 1==0
    %

end


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


% %% fcn_INTERNAL_loadExampleData
% function [seed_points, V, C] = fcn_INTERNAL_loadExampleData
%
%
% % pull halton set
% halton_points = haltonset(2);
% points_scrambled = scramble(halton_points,'RR2'); % scramble values
%
% % pick values from halton set
% Halton_range = [1801 1901];
% low_pt = Halton_range(1,1);
% high_pt = Halton_range(1,2);
% seed_points = points_scrambled(low_pt:high_pt,:);
% [V,C] = voronoin(seed_points);
% % V = V.*stretch;
% end % Ends fcn_INTERNAL_loadExampleData
