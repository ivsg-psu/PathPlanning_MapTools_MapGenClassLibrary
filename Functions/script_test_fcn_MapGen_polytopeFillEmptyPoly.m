% script_test_fcn_MapGen_polytopeFillEmptyPoly
% Tests: fcn_MapGen_polytopeFillEmptyPoly

% Revision history
% 2021_07_11 by Sean Brennan
% -- first write of script

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

%% DEMO case: basic call to function
fig_num = 10001;
titleString = sprintf('DEMO case: basic call to function');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Call the function
emptyPolytope = fcn_MapGen_polytopeFillEmptyPoly((fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(emptyPolytope));
assert(isfield(emptyPolytope,'vertices'));
assert(isfield(emptyPolytope,'xv'));
assert(isfield(emptyPolytope,'yv'));
assert(isfield(emptyPolytope,'distances'));
assert(isfield(emptyPolytope,'mean'));
assert(isfield(emptyPolytope,'area'));
assert(isfield(emptyPolytope,'max_radius'));
assert(isfield(emptyPolytope,'min_radius'));
assert(isfield(emptyPolytope,'mean_radius'));
assert(isfield(emptyPolytope,'radii'));
assert(isfield(emptyPolytope,'cost'));
assert(isfield(emptyPolytope,'parent_poly_id'));

% Check variable sizes
assert(isequal(1,length(emptyPolytope))); 

% Check variable values
assert(isempty(emptyPolytope.vertices));
assert(isempty(emptyPolytope.xv));
assert(isempty(emptyPolytope.yv));
assert(isempty(emptyPolytope.distances));
assert(isempty(emptyPolytope.mean));
assert(isempty(emptyPolytope.area));
assert(isempty(emptyPolytope.max_radius));
assert(isempty(emptyPolytope.min_radius));
assert(isempty(emptyPolytope.mean_radius));
assert(isempty(emptyPolytope.radii));
assert(isempty(emptyPolytope.cost));
assert(isempty(emptyPolytope.parent_poly_id));

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

% %% TEST case: This one returns nothing since there is no portion of the path in criteria
% fig_num = 20001;
% titleString = sprintf('TEST case: This one returns nothing since there is no portion of the path in criteria');
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

% Call the function
emptyPolytope = fcn_MapGen_polytopeFillEmptyPoly(([]));

% Check variable types
assert(isstruct(emptyPolytope));
assert(isfield(emptyPolytope,'vertices'));
assert(isfield(emptyPolytope,'xv'));
assert(isfield(emptyPolytope,'yv'));
assert(isfield(emptyPolytope,'distances'));
assert(isfield(emptyPolytope,'mean'));
assert(isfield(emptyPolytope,'area'));
assert(isfield(emptyPolytope,'max_radius'));
assert(isfield(emptyPolytope,'min_radius'));
assert(isfield(emptyPolytope,'mean_radius'));
assert(isfield(emptyPolytope,'radii'));
assert(isfield(emptyPolytope,'cost'));
assert(isfield(emptyPolytope,'parent_poly_id'));

% Check variable sizes
assert(isequal(1,length(emptyPolytope))); 

% Check variable values
assert(isempty(emptyPolytope.vertices));
assert(isempty(emptyPolytope.xv));
assert(isempty(emptyPolytope.yv));
assert(isempty(emptyPolytope.distances));
assert(isempty(emptyPolytope.mean));
assert(isempty(emptyPolytope.area));
assert(isempty(emptyPolytope.max_radius));
assert(isempty(emptyPolytope.min_radius));
assert(isempty(emptyPolytope.mean_radius));
assert(isempty(emptyPolytope.radii));
assert(isempty(emptyPolytope.cost));
assert(isempty(emptyPolytope.parent_poly_id));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Call the function
emptyPolytope = fcn_MapGen_polytopeFillEmptyPoly((-1));

% Check variable types
assert(isstruct(emptyPolytope));
assert(isfield(emptyPolytope,'vertices'));
assert(isfield(emptyPolytope,'xv'));
assert(isfield(emptyPolytope,'yv'));
assert(isfield(emptyPolytope,'distances'));
assert(isfield(emptyPolytope,'mean'));
assert(isfield(emptyPolytope,'area'));
assert(isfield(emptyPolytope,'max_radius'));
assert(isfield(emptyPolytope,'min_radius'));
assert(isfield(emptyPolytope,'mean_radius'));
assert(isfield(emptyPolytope,'radii'));
assert(isfield(emptyPolytope,'cost'));
assert(isfield(emptyPolytope,'parent_poly_id'));

% Check variable sizes
assert(isequal(1,length(emptyPolytope))); 

% Check variable values
assert(isempty(emptyPolytope.vertices));
assert(isempty(emptyPolytope.xv));
assert(isempty(emptyPolytope.yv));
assert(isempty(emptyPolytope.distances));
assert(isempty(emptyPolytope.mean));
assert(isempty(emptyPolytope.area));
assert(isempty(emptyPolytope.max_radius));
assert(isempty(emptyPolytope.min_radius));
assert(isempty(emptyPolytope.mean_radius));
assert(isempty(emptyPolytope.radii));
assert(isempty(emptyPolytope.cost));
assert(isempty(emptyPolytope.parent_poly_id));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

% Fill in seed points, V, and C
[seed_points, V, C] = fcn_INTERNAL_loadExampleData;

% fill polytopes from tiling
AABB = [0 0 1 1];
stretch = [1 1];
flag_removeEdgePolytopes = 1; % do NOT fill in polytopes to edge

Niterations = 10;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    emptyPolytope = fcn_MapGen_polytopeFillEmptyPoly(([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    emptyPolytope = fcn_MapGen_polytopeFillEmptyPoly((-1));
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
function [seed_points, V, C] = fcn_INTERNAL_loadExampleData


% pull halton set
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values

% pick values from halton set
Halton_range = [1801 1901];
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
seed_points = points_scrambled(low_pt:high_pt,:);
[V,C] = voronoin(seed_points);
% V = V.*stretch;
end % Ends fcn_INTERNAL_loadExampleData