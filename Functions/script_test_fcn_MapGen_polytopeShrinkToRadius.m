% script_test_fcn_MapGen_polytopesShrinkToRadius
% Tests function: fcn_MapGen_polytopesShrinkToRadius

% REVISION HISTORY:
% 
% 2021_06_09 by Sean Brennan, sbrennan@psu.edu
% - first written by S. Brennan using
%    % script_test_fcn_MapGen_polytopesDeleteByAABB as a template
% 
% 2025_07_11 by Sean Brennan, sbrennan@psu.edu
% - updated script testing to standard form
% 
% 2025_11_20 by Sean Brennan, sbrennan@psu.edu
% - Updated rev history to be in Markdown format
% - Replaced fig_+num with figNum

% TO-DO:
% 
% 2025_11_20 by Sean Brennan, sbrennan@psu.edu
% - fill in to-do items here.


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

%% DEMO case: Basic example of radius-specified shrinking
figNum = 10001;
titleString = sprintf('DEMO case: Basic example of radius-specified shrinking');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

shrinker = fcn_INTERNAL_loadExampleData;

orig_radius = shrinker.max_radius;
ratio = 0.5;
newRadius = orig_radius*ratio;

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

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% DEMO case: Iterative example of radius-specified shrinking
figNum = 10002;
titleString = sprintf('DEMO case: Iterative example of radius-specified shrinking');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

shrinker = fcn_INTERNAL_loadExampleData;

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
% figNum = 20001;
% titleString = sprintf('TEST case: simple crossing at origin');
% fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
% figure(figNum); clf;


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
figNum = 80001;
fprintf(1,'Figure: %.0f: FAST mode, empty figNum\n',figNum);
figure(figNum); close(figNum);

shrinker = fcn_INTERNAL_loadExampleData;

orig_radius = shrinker.max_radius;
ratio = 0.5;
newRadius = orig_radius*ratio;

% Call the function
shrunkPolytope = fcn_MapGen_polytopeShrinkToRadius(shrinker, newRadius, ([]));

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

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: FAST mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

shrinker = fcn_INTERNAL_loadExampleData;

orig_radius = shrinker.max_radius;
ratio = 0.5;
newRadius = orig_radius*ratio;

% Call the function
shrunkPolytope = fcn_MapGen_polytopeShrinkToRadius(shrinker, newRadius, (-1));

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
% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',figNum);
figure(figNum);
close(figNum);

shrinker = fcn_INTERNAL_loadExampleData;

orig_radius = shrinker.max_radius;
ratio = 0.5;
newRadius = orig_radius*ratio;

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    shrunkPolytope = fcn_MapGen_polytopeShrinkToRadius(shrinker, newRadius, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    shrunkPolytope = fcn_MapGen_polytopeShrinkToRadius(shrinker, newRadius, (-1));
end
fast_method = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));

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
assert(~any(figHandles==figNum));


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
function shrinker = fcn_INTERNAL_loadExampleData

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
end % Ends fcn_INTERNAL_loadExampleData