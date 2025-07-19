% script_test_fcn_MapGen_polytopesExpandEvenly
% Tests: fcn_MapGen_polytopesExpandEvenly

%
% REVISION HISTORY:
%§
% 2018_11_17, Seth Tau
% -- first write of script
% 2021_04_28, Seth Tau
% -- Adjusted example code ,
% 2021_06_26 S. Brennan
% -- Rebased code
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

%% DEMO case: expansion of a single polytope
fig_num = 10001;
titleString = sprintf('DEMO case: expansion of a single polytope');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

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
expandedPolytopes = fcn_MapGen_polytopesExpandEvenly(polytopes, expansionDistance, fig_num);

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
assert(isequal(get(gcf,'Number'),fig_num));

%% DEMO case: expansion of a polytope field
fig_num = 10002;
titleString = sprintf('DEMO case: expansion of a polytope field');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

map_name = "HST 30 450 SQT 0 1 0 1 SMV 0.02 0.005 1e-6 1234";
plot_flag = 1; disp_name = [1, 0.05 -0.05, 0.5 0.5 0.5, 10];
line_style = '-'; line_width = 2; color = [0 0 1];
axis_limits = [0 1 -0.1 1]; axis_style = 'square';
fill_info = [1 1 0 1 0.5];
fig_num = 7;

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
expandedPolytopes = fcn_MapGen_polytopesExpandEvenly(polytopes, expansionDistance, fig_num);

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
expandedPolytopes = fcn_MapGen_polytopesExpandEvenly(polytopes, expansionDistance, []);

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

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

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
expandedPolytopes = fcn_MapGen_polytopesExpandEvenly(polytopes, expansionDistance, -1);

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

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

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

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    expandedPolytopes = fcn_MapGen_polytopesExpandEvenly(polytopes, expansionDistance, []);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    expandedPolytopes = fcn_MapGen_polytopesExpandEvenly(polytopes, expansionDistance, []);
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