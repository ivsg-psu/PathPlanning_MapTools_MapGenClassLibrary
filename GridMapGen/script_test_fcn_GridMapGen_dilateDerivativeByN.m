% script_test_fcn_GridMapGen_dilateDerivativeByN
% Tests: fcn_GridMapGen_dilateDerivativeByN

%
% REVISION HISTORY:
%
% 2025_07_21 by Sean Brennan
% -- first write of function

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

%% DEMO case: random 100x100 matrix, 10% occupied, dilation of 1
fig_num = 10001;
titleString = sprintf('DEMO case: random 100x100 matrix, 10%% occupied, dilation of 1');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;


% Produce a smooth, random map
[~, inputMatrix]  = fcn_GridMapGen_generateRandomOccupancyMap(...
    'dilationLevel',60,.... % [1x1] strictly positive int
    'figNum',-1);

derivativeOrder = 1;

% Call the function
[rowDerivative, colDerivative, leftDilationMultiplier, rightDilationMultiplier] = ...
    fcn_GridMapGen_dilateDerivativeByN(inputMatrix, derivativeOrder, ...
    ([]), ([]), (fig_num));

sgtitle(titleString, 'Interpreter','none','FontSize',12)

% Check variable types
assert(isnumeric(rowDerivative));
assert(isnumeric(colDerivative));

assert(isnumeric(leftDilationMultiplier));
assert(isnumeric(leftDilationMultiplier));

% Check variable sizes
nRows = size(inputMatrix,1);
mCols = size(inputMatrix,2);
assert(isequal(size(rowDerivative), [nRows mCols]));
assert(isequal(size(colDerivative), [nRows mCols]));

assert(isequal(size(leftDilationMultiplier),[nRows nRows]));
assert(isequal(size(rightDilationMultiplier),[mCols mCols]));

% Check variable values
% assert(all(all(rowDerivative==2)));
% assert(all(all(colDerivative==0)));

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

%% TEST case: 10x10 matrix, derivate in cols exactly 1
fig_num = 20001;
titleString = sprintf('TEST case: 10x10 matrix, derivate in cols exactly 1');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;


inputMatrix = repmat((1:10),10,1);

derivativeOrder = 1;

% Call the function
[rowDerivative, colDerivative, leftDilationMultiplier, rightDilationMultiplier] = ...
    fcn_GridMapGen_dilateDerivativeByN(inputMatrix, derivativeOrder, ...
    ([]), ([]), (fig_num));

sgtitle(titleString, 'Interpreter','none','FontSize',12)

% Check variable types
assert(isnumeric(rowDerivative));
assert(isnumeric(colDerivative));

assert(isnumeric(leftDilationMultiplier));
assert(isnumeric(leftDilationMultiplier));

% Check variable sizes
nRows = size(inputMatrix,1);
mCols = size(inputMatrix,2);
assert(isequal(size(rowDerivative), [nRows mCols]));
assert(isequal(size(colDerivative), [nRows mCols]));

assert(isequal(size(leftDilationMultiplier),[nRows nRows]));
assert(isequal(size(rightDilationMultiplier),[mCols mCols]));

% Check variable values
assert(all(all(rowDerivative==0)));
assert(all(all(colDerivative==1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: 10x10 matrix, derivate in cols exactly 2
fig_num = 20002;
titleString = sprintf('TEST case: 10x10 matrix, derivate in cols exactly 2');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;


inputMatrix = repmat((2:2:20),10,1);

derivativeOrder = 1;

% Call the function
[rowDerivative, colDerivative, leftDilationMultiplier, rightDilationMultiplier] = ...
    fcn_GridMapGen_dilateDerivativeByN(inputMatrix, derivativeOrder, ...
    ([]), ([]), (fig_num));

sgtitle(titleString, 'Interpreter','none','FontSize',12)

% Check variable types
assert(isnumeric(rowDerivative));
assert(isnumeric(colDerivative));

assert(isnumeric(leftDilationMultiplier));
assert(isnumeric(leftDilationMultiplier));

% Check variable sizes
nRows = size(inputMatrix,1);
mCols = size(inputMatrix,2);
assert(isequal(size(rowDerivative), [nRows mCols]));
assert(isequal(size(colDerivative), [nRows mCols]));

assert(isequal(size(leftDilationMultiplier),[nRows nRows]));
assert(isequal(size(rightDilationMultiplier),[mCols mCols]));

% Check variable values
assert(all(all(rowDerivative==0)));
assert(all(all(colDerivative==2)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: 10x10 matrix, derivate in rows exactly 1
fig_num = 20003;
titleString = sprintf('TEST case: 10x10 matrix, derivate in rows exactly 1');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;


inputMatrix = repmat((1:10)',1,10);

derivativeOrder = 1;

% Call the function
[rowDerivative, colDerivative, leftDilationMultiplier, rightDilationMultiplier] = ...
    fcn_GridMapGen_dilateDerivativeByN(inputMatrix, derivativeOrder, ...
    ([]), ([]), (fig_num));

sgtitle(titleString, 'Interpreter','none','FontSize',12)

% Check variable types
assert(isnumeric(rowDerivative));
assert(isnumeric(colDerivative));

assert(isnumeric(leftDilationMultiplier));
assert(isnumeric(leftDilationMultiplier));

% Check variable sizes
nRows = size(inputMatrix,1);
mCols = size(inputMatrix,2);
assert(isequal(size(rowDerivative), [nRows mCols]));
assert(isequal(size(colDerivative), [nRows mCols]));

assert(isequal(size(leftDilationMultiplier),[nRows nRows]));
assert(isequal(size(rightDilationMultiplier),[mCols mCols]));

% Check variable values
assert(all(all(rowDerivative==1)));
assert(all(all(colDerivative==0)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: 10x10 matrix, derivate in rows exactly 2
fig_num = 20004;
titleString = sprintf('TEST case: 10x10 matrix, derivate in rows exactly 2');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;


inputMatrix = repmat((2:2:20)',1,10);

derivativeOrder = 1;

% Call the function
[rowDerivative, colDerivative, leftDilationMultiplier, rightDilationMultiplier] = ...
    fcn_GridMapGen_dilateDerivativeByN(inputMatrix, derivativeOrder, ...
    ([]), ([]), (fig_num));

sgtitle(titleString, 'Interpreter','none','FontSize',12)

% Check variable types
assert(isnumeric(rowDerivative));
assert(isnumeric(colDerivative));

assert(isnumeric(leftDilationMultiplier));
assert(isnumeric(leftDilationMultiplier));

% Check variable sizes
nRows = size(inputMatrix,1);
mCols = size(inputMatrix,2);
assert(isequal(size(rowDerivative), [nRows mCols]));
assert(isequal(size(colDerivative), [nRows mCols]));

assert(isequal(size(leftDilationMultiplier),[nRows nRows]));
assert(isequal(size(rightDilationMultiplier),[mCols mCols]));

% Check variable values
assert(all(all(rowDerivative==2)));
assert(all(all(colDerivative==0)));

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

% Produce a smooth, random map
[~, inputMatrix]  = fcn_GridMapGen_generateRandomOccupancyMap(...
    'dilationLevel',60,.... % [1x1] strictly positive int
    'figNum',-1);

derivativeOrder = 1;

% Call the function
[rowDerivative, colDerivative, leftDilationMultiplier, rightDilationMultiplier] = ...
    fcn_GridMapGen_dilateDerivativeByN(inputMatrix, derivativeOrder, ...
    ([]), ([]), ([]));

% Check variable types
assert(isnumeric(rowDerivative));
assert(isnumeric(colDerivative));

assert(isnumeric(leftDilationMultiplier));
assert(isnumeric(leftDilationMultiplier));

% Check variable sizes
nRows = size(inputMatrix,1);
mCols = size(inputMatrix,2);
assert(isequal(size(rowDerivative), [nRows mCols]));
assert(isequal(size(colDerivative), [nRows mCols]));

assert(isequal(size(leftDilationMultiplier),[nRows nRows]));
assert(isequal(size(rightDilationMultiplier),[mCols mCols]));

% Check variable values
% assert(all(all(rowDerivative==2)));
% assert(all(all(colDerivative==0)));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Produce a smooth, random map
[~, inputMatrix]  = fcn_GridMapGen_generateRandomOccupancyMap(...
    'dilationLevel',60,.... % [1x1] strictly positive int
    'figNum',-1);

derivativeOrder = 1;

% Call the function
[rowDerivative, colDerivative, leftDilationMultiplier, rightDilationMultiplier] = ...
    fcn_GridMapGen_dilateDerivativeByN(inputMatrix, derivativeOrder, ...
    ([]), ([]), (-1));

% Check variable types
assert(isnumeric(rowDerivative));
assert(isnumeric(colDerivative));

assert(isnumeric(leftDilationMultiplier));
assert(isnumeric(leftDilationMultiplier));

% Check variable sizes
nRows = size(inputMatrix,1);
mCols = size(inputMatrix,2);
assert(isequal(size(rowDerivative), [nRows mCols]));
assert(isequal(size(colDerivative), [nRows mCols]));

assert(isequal(size(leftDilationMultiplier),[nRows nRows]));
assert(isequal(size(rightDilationMultiplier),[mCols mCols]));

% Check variable values
% assert(all(all(rowDerivative==2)));
% assert(all(all(colDerivative==0)));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% FAST mode comparisons in normal mode
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons in normal mode\n',fig_num);
figure(fig_num);
close(fig_num);

% Produce a smooth, random map
[~, inputMatrix]  = fcn_GridMapGen_generateRandomOccupancyMap(...
    'dilationLevel',60,.... % [1x1] strictly positive int
    'figNum',-1);

derivativeOrder = 1;
     
Niterations = 10;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [rowDerivative, colDerivative, leftDilationMultiplier, rightDilationMultiplier] = ...
        fcn_GridMapGen_dilateDerivativeByN(inputMatrix, derivativeOrder, ...
        ([]), ([]), ([]));
end
slow_method = toc;

% Do calculation without pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [rowDerivative, colDerivative, leftDilationMultiplier, rightDilationMultiplier] = ...
        fcn_GridMapGen_dilateDerivativeByN(inputMatrix, derivativeOrder, ...
        ([]), ([]), (-1));
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