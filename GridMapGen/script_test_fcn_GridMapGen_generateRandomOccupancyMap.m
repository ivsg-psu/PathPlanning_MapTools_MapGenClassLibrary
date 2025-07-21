% script_test_fcn_GridMapGen_generateRandomOccupancyMap
% Tests: fcn_GridMapGen_generateRandomOccupancyMap

%
% REVISION HISTORY:
%
% 2008_10_18 by Sean Brennan
% -- first write of function
% 2025_07_17 by Sean Brennan
% -- imported into MapGen library with updates to formatting

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
% DEMO figures start with 1

close all;
fprintf(1,'Figure: 1XXXXXX: DEMO cases\n');

%% DEMO case: occupancyMap using all defaults blank
fig_num = 10001;
titleString = sprintf('DEMO case: occupancyMap using all defaults blank');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set input arguments
mapSize = [];
occupancyRatio = [];
dilationLevel = [];
seedMap = [];

% Call the function
[occupancyMatrix, randomMatrixDilated, optimizedThreshold] = ...
    fcn_GridMapGen_generateRandomOccupancyMap(...
    'mapSize', (mapSize), ...
    'occupancyRatio',(occupancyRatio), ...
    'dilationLevel',(dilationLevel), ...
    'seedMap', (seedMap), ...
    'figNum',(fig_num));

sgtitle(titleString, 'Interpreter','none','Fontsize',10);

fcn_INTERNAL_checkData(occupancyMatrix, randomMatrixDilated, optimizedThreshold, mapSize, occupancyRatio);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% DEMO case: occupancyMap from no input arguments (no figure opens)
fig_num = 10002;
titleString = sprintf('DEMO case: occupancyMap from no input arguments (no figure opens)');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); close(fig_num);

% Set input arguments
% mapSize = [];
% mColumns = [];
% occupancyRatio = [];
% dilationLevel = [];
% seedMap = [];

% Call the function
[occupancyMatrix, randomMatrixDilated, optimizedThreshold]  = fcn_GridMapGen_generateRandomOccupancyMap;
% (...
%     (mapSize), ...
%     (mColumns), ...
%     (occupancyRatio), ...
%     (dilationLevel), ...
%     (seedMap), ...
%     (fig_num));

fcn_INTERNAL_checkData(occupancyMatrix, randomMatrixDilated, optimizedThreshold, mapSize, occupancyRatio);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% DEMO case: huge smoothing (dilationLevel = 40)
fig_num = 10003;
titleString = sprintf('DEMO case: huge smoothing (dilationLevel = 40)');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set input arguments
mapSize = [];
occupancyRatio = [];
dilationLevel = 40;
seedMap = [];

% Call the function
[occupancyMatrix, randomMatrixDilated, optimizedThreshold] = ...
    fcn_GridMapGen_generateRandomOccupancyMap(...
    'mapSize', (mapSize), ...
    'occupancyRatio',(occupancyRatio), ...
    'dilationLevel',(dilationLevel), ...
    'seedMap', (seedMap), ...
    'figNum',(fig_num));

sgtitle(titleString, 'Interpreter','none','Fontsize',10);

fcn_INTERNAL_checkData(occupancyMatrix, randomMatrixDilated, optimizedThreshold, mapSize, occupancyRatio);

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
% TEST figures start with 2

close all;
fprintf(1,'Figure: 2XXXXXX: TEST mode cases\n');

%% TEST case: using all inputs specified to defaults
fig_num = 20001;
titleString = sprintf('TEST case: occupancyMap using all inputs specified to defaults');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set input arguments
mapSize = [100 100];
occupancyRatio = 0.2;
dilationLevel = 2;
seedMap = 234;
leftDilationMultiplier = [];
rightDilationMultiplier = [];
threshold = [];
flagSkipThresholdOptimization = 0;

% Call the function
[occupancyMatrix, randomMatrixDilated, optimizedThreshold, ~, ~]  = ...
    fcn_GridMapGen_generateRandomOccupancyMap(...
    'mapSize', (mapSize),... % [nRows mCols])
    'occupancyRatio',(occupancyRatio),... % [1x1] value between 0 and 1
    'dilationLevel',(dilationLevel),.... % [1x1] strictly positive int
    'seedMap', (seedMap),... % [1x1] integer to be a random seed or NxM matrix of random numbers
    'leftDilationMultiplier', (leftDilationMultiplier),... %  [nRows nRows], ...
    'rightDilationMultiplier', (rightDilationMultiplier),... % [mCols mCols], ...
    'thresholdForced', (threshold), ... % [1x1] scalar
    'flagSkipThresholdOptimization',(flagSkipThresholdOptimization),...% [1x1] scalar
    'figNum',(fig_num));

sgtitle(titleString, 'Interpreter','none','Fontsize',10);

fcn_INTERNAL_checkData(occupancyMatrix, randomMatrixDilated, optimizedThreshold, mapSize, occupancyRatio);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: mapSize specified, producing 50x50 map
fig_num = 20002;
titleString = sprintf('TEST case: occupancyMap with mapSize specified, mColumns empty so inherets mapSize value producing 50x50 map');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set input arguments
mapSize = [50 50];
occupancyRatio = 0.2;
dilationLevel = 2;
seedMap = 234;

% Call the function
[occupancyMatrix, randomMatrixDilated, optimizedThreshold] = ...
    fcn_GridMapGen_generateRandomOccupancyMap(...
    'mapSize', (mapSize), ...
    'occupancyRatio',(occupancyRatio), ...
    'dilationLevel',(dilationLevel), ...
    'seedMap', (seedMap), ...
    'figNum',(fig_num));

sgtitle(titleString, 'Interpreter','none','Fontsize',10);

fcn_INTERNAL_checkData(occupancyMatrix, randomMatrixDilated, optimizedThreshold, mapSize, occupancyRatio);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: mapSize and mColumns specified differently, producing 50x100 map
fig_num = 20003;
titleString = sprintf('TEST case: mapSize and mColumns specified differently, producing 50x100 map');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set input arguments
mapSize = [50 100];
occupancyRatio = 0.2;
dilationLevel = 2;
seedMap = 234;

% Call the function
[occupancyMatrix, randomMatrixDilated, optimizedThreshold] = ...
    fcn_GridMapGen_generateRandomOccupancyMap(...
    'mapSize', (mapSize), ...
    'occupancyRatio',(occupancyRatio), ...
    'dilationLevel',(dilationLevel), ...
    'seedMap', (seedMap), ...
    'figNum',(fig_num));

sgtitle(titleString, 'Interpreter','none','Fontsize',10);

fcn_INTERNAL_checkData(occupancyMatrix, randomMatrixDilated, optimizedThreshold, mapSize, occupancyRatio);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: mapSize not specified, but mColumns specified differently, producing 100x50 map
fig_num = 20004;
titleString = sprintf('TEST case: mapSize and mColumns specified differently, producing 100x50 map');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set input arguments
mapSize = [];
occupancyRatio = 0.2;
dilationLevel = 2;
seedMap = 234;

% Call the function
[occupancyMatrix, randomMatrixDilated, optimizedThreshold] = ...
    fcn_GridMapGen_generateRandomOccupancyMap(...
    'mapSize', (mapSize), ...
    'occupancyRatio',(occupancyRatio), ...
    'dilationLevel',(dilationLevel), ...
    'seedMap', (seedMap), ...
    'figNum',(fig_num));

sgtitle(titleString, 'Interpreter','none','Fontsize',10);

fcn_INTERNAL_checkData(occupancyMatrix, randomMatrixDilated, optimizedThreshold, mapSize, occupancyRatio);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: occupancyRatio increased to 0.4
fig_num = 20005;
titleString = sprintf('TEST case: occupancyRatio increased to 0.4');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set input arguments
mapSize = [];
occupancyRatio = 0.4;
dilationLevel = 2;
seedMap = 234;

% Call the function
[occupancyMatrix, randomMatrixDilated, optimizedThreshold] = ...
    fcn_GridMapGen_generateRandomOccupancyMap(...
    'mapSize', (mapSize), ...
    'occupancyRatio',(occupancyRatio), ...
    'dilationLevel',(dilationLevel), ...
    'seedMap', (seedMap), ...
    'figNum',(fig_num));

sgtitle(titleString, 'Interpreter','none','Fontsize',10);

fcn_INTERNAL_checkData(occupancyMatrix, randomMatrixDilated, optimizedThreshold, mapSize, occupancyRatio);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: seedMap given as a seed gives repeatedly same map
fig_num = 20006;
titleString = sprintf('TEST case: seedMap given as a seed gives repeatedly same map');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set input arguments
mapSize = [];
occupancyRatio = [];
dilationLevel = 2;
seedMap = 101;

% Call the function
[occupancyMatrix, randomMatrixDilated, optimizedThreshold] = ...
    fcn_GridMapGen_generateRandomOccupancyMap(...
    'mapSize', (mapSize), ...
    'occupancyRatio',(occupancyRatio), ...
    'dilationLevel',(dilationLevel), ...
    'seedMap', (seedMap), ...
    'figNum',(fig_num));

sgtitle(titleString, 'Interpreter','none','Fontsize',10);

fcn_INTERNAL_checkData(occupancyMatrix, randomMatrixDilated, optimizedThreshold, mapSize, occupancyRatio);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: seedMap given as a matrix gives repeatedly same map
fig_num = 20007;
titleString = sprintf('TEST case: seedMap given as a matrix gives repeatedly same map');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set input arguments
mapSize = [];
occupancyRatio = [];
dilationLevel = 2;
seedMap = rand(30,50);

% Call the function
[occupancyMatrix, randomMatrixDilated, optimizedThreshold] = ...
    fcn_GridMapGen_generateRandomOccupancyMap(...
    'mapSize', (mapSize), ...
    'occupancyRatio',(occupancyRatio), ...
    'dilationLevel',(dilationLevel), ...
    'seedMap', (seedMap), ...
    'figNum',(fig_num));

sgtitle(titleString, 'Interpreter','none','Fontsize',10);

fcn_INTERNAL_checkData(occupancyMatrix, randomMatrixDilated, optimizedThreshold, size(seedMap), occupancyRatio);

% Call the function again, same seed
[occupancyMatrix2, randomMatrixDilated, optimizedThreshold] = ...
    fcn_GridMapGen_generateRandomOccupancyMap(...
    'mapSize', (mapSize), ...
    'occupancyRatio',(occupancyRatio), ...
    'dilationLevel',(dilationLevel), ...
    'seedMap', (seedMap), ...
    'figNum',(fig_num));

sgtitle(titleString, 'Interpreter','none','Fontsize',10);

fcn_INTERNAL_checkData(occupancyMatrix2, randomMatrixDilated, optimizedThreshold, size(seedMap), occupancyRatio);

assert(isequal(occupancyMatrix,occupancyMatrix2));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



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

% Set input arguments
mapSize = [100 100];
occupancyRatio = 0.2;
dilationLevel = 2;
seedMap = 234;

% Call the function
[occupancyMatrix, randomMatrixDilated, optimizedThreshold] = ...
    fcn_GridMapGen_generateRandomOccupancyMap(...
    'mapSize', (mapSize), ...
    'occupancyRatio',(occupancyRatio), ...
    'dilationLevel',(dilationLevel), ...
    'seedMap', (seedMap), ...
    'figNum',([]));

fcn_INTERNAL_checkData(occupancyMatrix, randomMatrixDilated, optimizedThreshold, mapSize, occupancyRatio);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Set input arguments
mapSize = [100 100];
occupancyRatio = 0.2;
dilationLevel = 2;
seedMap = 234;

% Call the function
[occupancyMatrix, randomMatrixDilated, optimizedThreshold] = ...
    fcn_GridMapGen_generateRandomOccupancyMap(...
    'mapSize', (mapSize), ...
    'occupancyRatio',(occupancyRatio), ...
    'dilationLevel',(dilationLevel), ...
    'seedMap', (seedMap), ...
    'figNum',(-1));

fcn_INTERNAL_checkData(occupancyMatrix, randomMatrixDilated, optimizedThreshold, mapSize, occupancyRatio);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Compare speeds of fast / normal and pre-calculation versus post-calculation 
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

% Set input arguments
mapSize = [100 100];
occupancyRatio = 0.2;
dilationLevel = 40;
seedMap = 234;
leftDilationMultiplier = [];
rightDilationMultiplier = [];
threshold = [];
flagSkipThresholdOptimization = 0;

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [occupancyMatrix, randomMatrixDilated, optimizedThreshold] = ...
        fcn_GridMapGen_generateRandomOccupancyMap(...
        'mapSize', (mapSize), ...
        'occupancyRatio',(occupancyRatio), ...
        'dilationLevel',(dilationLevel), ...
        'seedMap', (seedMap), ...
        'figNum',([]));
end
slow_method = toc;


% Do calculation without pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [occupancyMatrix2, randomMatrixDilated, optimizedThreshold, leftDilationMultiplier, rightDilationMultiplier] = ...
        fcn_GridMapGen_generateRandomOccupancyMap(...
        'mapSize', (mapSize), ...
        'occupancyRatio',(occupancyRatio), ...
        'dilationLevel',(dilationLevel), ...
        'seedMap', (seedMap), ...
        'figNum',(-1));
end
fast_method = toc;
assert(isequal(occupancyMatrix,occupancyMatrix2));


% Do calculation with dilation pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function, using pre-calculated left and right matrices (from
    % earlier) for the dilation
    [occupancyMatrix3, randomMatrixDilated, optimizedThreshold, ~, ~]  = ...
        fcn_GridMapGen_generateRandomOccupancyMap(...
        'mapSize', (mapSize),... % [nRows mCols])
        'occupancyRatio',(occupancyRatio),... % [1x1] value between 0 and 1
        'dilationLevel',(dilationLevel),.... % [1x1] strictly positive int
        'seedMap', (seedMap),... % [1x1] integer to be a random seed or NxM matrix of random numbers
        'leftDilationMultiplier', (leftDilationMultiplier),... %  [nRows nRows], ...
        'rightDilationMultiplier', (rightDilationMultiplier),... % [mCols mCols], ...
        'threshold', (threshold), ... % [1x1] scalar
        'flagSkipThresholdOptimization',(flagSkipThresholdOptimization),...% [1x1] scalar
        'figNum',(-1));
end
very_fast_method = toc;
assert(isequal(occupancyMatrix,occupancyMatrix3));


% Do calculation with dilation pre-calculation, fixed threshold, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function, using pre-calculated left and right matrices (from
    % earlier) for the dilation
    [occupancyMatrix3, randomMatrixDilated, ~, ~, ~]  = ...
        fcn_GridMapGen_generateRandomOccupancyMap(...
        'mapSize', (mapSize),... % [nRows mCols])
        'occupancyRatio',(occupancyRatio),... % [1x1] value between 0 and 1
        'dilationLevel',(dilationLevel),.... % [1x1] strictly positive int
        'seedMap', (seedMap),... % [1x1] integer to be a random seed or NxM matrix of random numbers
        'leftDilationMultiplier', (leftDilationMultiplier),... %  [nRows nRows], ...
        'rightDilationMultiplier', (rightDilationMultiplier),... % [mCols mCols], ...
        'thresholdForced', (optimizedThreshold), ... % [1x1] scalar
        'flagSkipThresholdOptimization',(1),...% [1x1] scalar
        'figNum',(-1));
end
very_very_fast_method = toc;
assert(isequal(occupancyMatrix,occupancyMatrix3));


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

% Plot results as bar chart
figure(373737);
clf;
hold on;

names = {'Normal mode','Fast mode','Fast mode, pre-calc dilation', 'Fast mode, pre-calc dilation, forced threshold'};
X = categorical(names);
X = reordercats(X,names); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method very_fast_method very_very_fast_method]*1000/Niterations;
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

%% fcn_INTERNAL_checkData
function fcn_INTERNAL_checkData(occupancyMatrix, randomMatrixDilated, optimizedThreshold, mapSize, occupancyRatio)
assert(islogical(occupancyMatrix));
assert(isnumeric(randomMatrixDilated));
assert(isnumeric(optimizedThreshold));

% Check variable sizes
if ~exist('mapSize','var') || isempty(mapSize)
    assert(isequal(size(occupancyMatrix,1),100));
else
    assert(isequal(size(occupancyMatrix),mapSize));
end
assert(isequal(size(randomMatrixDilated),size(occupancyMatrix)));
assert(isscalar(optimizedThreshold));

% Check variable values
percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix, -1);
if ~exist('occupancyRatio','var') || isempty(occupancyRatio)
    assert(abs(percentOccupied - 0.2)<0.1);
else
    assert(abs(percentOccupied - occupancyRatio)<0.1);
end

end