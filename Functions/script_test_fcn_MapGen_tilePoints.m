% script_test_fcn_MapGen_tilePoints
% Tests function: fcn_MapGen_tilePoints

% REVISION HISTORY:
% 2023_02_24
% -- first written by S. Brennan
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

%% DEMO case: repeating set of numbers with tileDepth of 1
fig_num = 10001;
titleString = sprintf('DEMO case: repeating set of numbers with tileDepth of 1');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

Npoints = 20;
tileDepth = 1;
AABB = [0 0 2 2];
inputPoints = 2*rand(Npoints,2);

% Call the function
[tiledPoints] = fcn_MapGen_tilePoints(inputPoints, tileDepth, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(tiledPoints));

% Check variable sizes
Multiplier = (tileDepth*2+1)^2;
assert(length(tiledPoints(1,:))==2); % Is it 2 columns?
assert(length(tiledPoints(:,1))==Npoints*Multiplier); % Does it have right number of points?

% Check variable values
% (do this later)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));
 
%% DEMO case: points outside AABB - will tile them but they are shifted
fig_num = 10002;
titleString = sprintf('DEMO case: points outside AABB - will tile them but they are shifted');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

Npoints = 200;
tileDepth = 1;
AABB = [0 0 2 2];
inputPoints = 0.5*randn(Npoints,2)+[7,0] + 1;

% Call the function
[tiledPoints] = fcn_MapGen_tilePoints(inputPoints, tileDepth, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(tiledPoints));

% Check variable sizes
Multiplier = (tileDepth*2+1)^2;
assert(length(tiledPoints(1,:))==2); % Is it 2 columns?
assert(length(tiledPoints(:,1))==Npoints*Multiplier); % Does it have right number of points?

% Check variable values
% (do this later)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% DEMO case: AABB bigger than points pads points
fig_num = 10003;
titleString = sprintf('DEMO case: AABB bigger than points pads points');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

Npoints = 200;
tileDepth = 1;
AABB = [0 0 5 5];
inputPoints = 0.5*randn(Npoints,2)+1;

% Call the function
[tiledPoints] = fcn_MapGen_tilePoints(inputPoints, tileDepth, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(tiledPoints));

% Check variable sizes
Multiplier = (tileDepth*2+1)^2;
assert(length(tiledPoints(1,:))==2); % Is it 2 columns?
assert(length(tiledPoints(:,1))==Npoints*Multiplier); % Does it have right number of points?

% Check variable values
% (do this later)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% DEMO case: AABB shifted will still tile points
fig_num = 10004;
titleString = sprintf('DEMO case: AABB shifted will still tile points');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

Npoints = 200;
tileDepth = 1;
AABB = [0 0 2 2]- 7;
inputPoints = 0.5*randn(Npoints,2) + 1;

% Call the function
[tiledPoints] = fcn_MapGen_tilePoints(inputPoints, tileDepth, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(tiledPoints));

% Check variable sizes
Multiplier = (tileDepth*2+1)^2;
assert(length(tiledPoints(1,:))==2); % Is it 2 columns?
assert(length(tiledPoints(:,1))==Npoints*Multiplier); % Does it have right number of points?

% Check variable values
% (do this later)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% DEMO case: repeating set of numbers with tileDepth of 2
fig_num = 10005;
titleString = sprintf('DEMO case: repeating set of numbers with tileDepth of 1');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

Npoints = 20;
tileDepth = 2;
AABB = [0 0 2 2];
inputPoints = 2*rand(Npoints,2);

% Call the function
[tiledPoints] = fcn_MapGen_tilePoints(inputPoints, tileDepth, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(tiledPoints));

% Check variable sizes
Multiplier = (tileDepth*2+1)^2;
assert(length(tiledPoints(1,:))==2); % Is it 2 columns?
assert(length(tiledPoints(:,1))==Npoints*Multiplier); % Does it have right number of points?

% Check variable values
% (do this later)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% DEMO case: repeating set of numbers with tileDepth of 5
fig_num = 10006;
titleString = sprintf('DEMO case: repeating set of numbers with tileDepth of 5');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

Npoints = 20;
tileDepth = 5;
AABB = [0 0 2 2];
inputPoints = 2*rand(Npoints,2);

% Call the function
[tiledPoints] = fcn_MapGen_tilePoints(inputPoints, tileDepth, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(tiledPoints));

% Check variable sizes
Multiplier = (tileDepth*2+1)^2;
assert(length(tiledPoints(1,:))==2); % Is it 2 columns?
assert(length(tiledPoints(:,1))==Npoints*Multiplier); % Does it have right number of points?

% Check variable values
% (do this later)

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

Npoints = 20;
tileDepth = 1;
AABB = [0 0 2 2];
inputPoints = 2*rand(Npoints,2);

% Call the function
[tiledPoints] = fcn_MapGen_tilePoints(inputPoints, tileDepth, AABB, ([]));

% Check variable types
assert(isnumeric(tiledPoints));

% Check variable sizes
Multiplier = (tileDepth*2+1)^2;
assert(length(tiledPoints(1,:))==2); % Is it 2 columns?
assert(length(tiledPoints(:,1))==Npoints*Multiplier); % Does it have right number of points?

% Check variable values
% (do this later)

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

Npoints = 20;
tileDepth = 1;
AABB = [0 0 2 2];
inputPoints = 2*rand(Npoints,2);

% Call the function
[tiledPoints] = fcn_MapGen_tilePoints(inputPoints, tileDepth, AABB, (-1));

% Check variable types
assert(isnumeric(tiledPoints));

% Check variable sizes
Multiplier = (tileDepth*2+1)^2;
assert(length(tiledPoints(1,:))==2); % Is it 2 columns?
assert(length(tiledPoints(:,1))==Npoints*Multiplier); % Does it have right number of points?

% Check variable values
% (do this later)

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

Npoints = 20;
tileDepth = 1;
AABB = [0 0 2 2];
inputPoints = 2*rand(Npoints,2);

Niterations = 10;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [tiledPoints] = fcn_MapGen_tilePoints(inputPoints, tileDepth, AABB, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [tiledPoints] = fcn_MapGen_tilePoints(inputPoints, tileDepth, AABB, (-1));
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

if 1==0 % Fail cases
    %% Input incorrect, wrong number of argument
    inputPoints = rand(10,1);
    tileDepth = 2;
    fcn_MapGen_tilePoints(inputPoints,tileDepth);

    %% Input incorrect, wrong number of arguments
    
    inputPoints = rand(10,1);
    tileDepth = 2;
    AABB = [0 0 1 1];
    fcn_MapGen_tilePoints(inputPoints,tileDepth,AABB,3,4);

    %% Input incorrect for inputPoints, only 1 column
    
    inputPoints = rand(10,1);
    tileDepth = 2;
    AABB = [0 0 1 1];
    fcn_MapGen_tilePoints(inputPoints,tileDepth,AABB);
    %% Input incorrect for inputPoints, 3 column
    
    inputPoints = rand(10,1);
    tileDepth = 2;
    AABB = [0 0 1 1];
    fcn_MapGen_tilePoints(inputPoints,tileDepth,AABB);
    %% Input incorrect for tileDepth, not strictly positive
    
    inputPoints = rand(10,2);
    tileDepth = 0;
    AABB = [0 0 1 1];
    fcn_MapGen_tilePoints(inputPoints,tileDepth,AABB);
    %% Input incorrect for tileDepth, not an integer
    
    inputPoints = rand(10,2);
    tileDepth = 1.2;
    AABB = [0 0 1 1];
    fcn_MapGen_tilePoints(inputPoints,tileDepth,AABB);
    %% Input incorrect for AABB, not an 4x1
    
    inputPoints = rand(10,2);
    tileDepth = 1;
    AABB = [0 0 1];
    fcn_MapGen_tilePoints(inputPoints,tileDepth,AABB);
    %% Input incorrect for AABB, not an 4x1
    
    inputPoints = rand(10,2);
    tileDepth = 1;
    AABB = [0 0 1 1 1];
    fcn_MapGen_tilePoints(inputPoints,tileDepth,AABB);
    %% Input incorrect for AABB, not an 4x1
    
    inputPoints = rand(10,2);
    tileDepth = 1;
    AABB = [0 0 1 1; 0 0 0 0];
    fcn_MapGen_tilePoints(inputPoints,tileDepth,AABB);
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