% script_test_fcn_MapGen_generatePolysFromVoronoiAABB
% Tests: fcn_MapGen_generatePolysFromVoronoiAABB

%
% REVISION HISTORY:
%
% 2021_07_02 by Sean Brennan
% -- first write of script
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

%% DEMO case: basic demo filling polytopes from tiling
fig_num = 10001;
titleString = sprintf('DEMO case: basic demo filling polytopes from tiling');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

[seedPoints, V, C] = fcn_INTERNAL_loadExampleData;

AABB = [0 0 1 1];
stretch = [1 1];

% Call the function
polytopes = fcn_MapGen_generatePolysFromVoronoiAABB(seedPoints,V,C,AABB, stretch, (fig_num));

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
Npolytopes = length(seedPoints(:,1));
assert(isequal(Npolytopes,length(polytopes))); 

% Check variable values
% (cannot check - random)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

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
% 
% %% TEST case: simple crossing at origin
% fig_num = 20001;
% titleString = sprintf('TEST case: simple crossing at origin');
% fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
% figure(fig_num); clf;


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

[seedPoints, V, C] = fcn_INTERNAL_loadExampleData;

AABB = [0 0 1 1];
stretch = [1 1];

% Call the function
polytopes = fcn_MapGen_generatePolysFromVoronoiAABB(seedPoints,V,C,AABB, stretch, ([]));

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
Npolytopes = length(seedPoints(:,1));
assert(isequal(Npolytopes,length(polytopes))); 

% Check variable values
% (cannot check - random)

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

[seedPoints, V, C] = fcn_INTERNAL_loadExampleData;

AABB = [0 0 1 1];
stretch = [1 1];

% Call the function
polytopes = fcn_MapGen_generatePolysFromVoronoiAABB(seedPoints,V,C,AABB, stretch, (-1));

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
Npolytopes = length(seedPoints(:,1));
assert(isequal(Npolytopes,length(polytopes))); 

% Check variable values
% (cannot check - random)

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

[seedPoints, V, C] = fcn_INTERNAL_loadExampleData;

AABB = [0 0 1 1];
stretch = [1 1];

Niterations = 20;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    polytopes = fcn_MapGen_generatePolysFromVoronoiAABB(seedPoints,V,C,AABB, stretch, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    polytopes = fcn_MapGen_generatePolysFromVoronoiAABB(seedPoints,V,C,AABB, stretch, (-1));
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

%% fcn_INTERNAL_loadExampleData
function [seedPoints, V, C] = fcn_INTERNAL_loadExampleData
% pull halton set
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values

% pick values from halton set
Halton_range = [1801 1901];
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
seedPoints = points_scrambled(low_pt:high_pt,:);
[V,C] = voronoin(seedPoints);

end % Ends fcn_INTERNAL_loadExampleData