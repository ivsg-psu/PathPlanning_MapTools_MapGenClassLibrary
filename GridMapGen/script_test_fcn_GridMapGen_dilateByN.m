% script_test_fcn_GridMapGen_dilateByN
% Tests: fcn_GridMapGen_dilateByN

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

%% DEMO case: random 100x100 matrix, 10% occupied, dilation of 1
fig_num = 10001;
titleString = sprintf('DEMO case: random 100x100 matrix, 10%% occupied, dilation of 1');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

occupancyMatrix = rand(100,100)<0.1; % about 10 percent occupied
% percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix, (28282));
dilationLevel = 1;

% Call the function
dilatedMatrix = fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(dilatedMatrix));

% Check variable sizes
assert(isequal(size(dilatedMatrix),size(occupancyMatrix)));

% Check variable values
assert(all(all(dilatedMatrix>=occupancyMatrix)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% DEMO case: random 100x100 matrix, 2% occupied, dilation of 2
fig_num = 10001;
titleString = sprintf('DEMO case: random 100x100 matrix, 2%% occupied, dilation of 2');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

occupancyMatrix = rand(100,100)<0.02; % about 2 percent occupied
% percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix, (28282));
dilationLevel = 2;

% Call the function
dilatedMatrix = fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(dilatedMatrix));

% Check variable sizes
assert(isequal(size(dilatedMatrix),size(occupancyMatrix)));

% Check variable values
assert(all(all(dilatedMatrix>=occupancyMatrix)));

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

%% TEST case: 10x10 matrix, exactly 50% occupied by column, converted to exactly 60% by 1-step dilation
fig_num = 20001;
titleString = sprintf('TEST case: 10x10 matrix, exactly 50%% occupied by column, converted to exactly 60%% by 1-step dilation');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

occupancyMatrix = [zeros(10,5) ones(10,5)];
dilationLevel = 1;

% Call the function
dilatedMatrix = fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(dilatedMatrix));

% Check variable sizes
assert(isequal(size(dilatedMatrix),size(occupancyMatrix)));

% Check variable values
assert(all(all(dilatedMatrix>=occupancyMatrix)));
percentOccupied = fcn_GridMapGen_dilateOccupancyStats(dilatedMatrix,-1);
assert(isequal(round(percentOccupied,6),0.6));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: 10x10 matrix, exactly 50% occupied by column, converted to exactly 70% by 2-step dilation
fig_num = 20001;
titleString = sprintf('TEST case: 10x10 matrix, exactly 50%% occupied by column, converted to exactly 70%% by 2-step dilation');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

occupancyMatrix = [zeros(10,5) ones(10,5)];
dilationLevel = 2;

% Call the function
dilatedMatrix = fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(dilatedMatrix));

% Check variable sizes
assert(isequal(size(dilatedMatrix),size(occupancyMatrix)));

% Check variable values
assert(all(all(dilatedMatrix>=occupancyMatrix)));
percentOccupied = fcn_GridMapGen_dilateOccupancyStats(dilatedMatrix,-1);
assert(isequal(round(percentOccupied,6),0.7));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% TEST case: 10x10 matrix, exactly 50% occupied by row, converted to exactly 70% by 2-step dilation
fig_num = 20001;
titleString = sprintf('TEST case: 10x10 matrix, exactly 50%% occupied by row, converted to exactly 70%% by 2-step dilation');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

occupancyMatrix = [zeros(5,10); ones(5,10)];
dilationLevel = 2;

% Call the function
dilatedMatrix = fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(dilatedMatrix));

% Check variable sizes
assert(isequal(size(dilatedMatrix),size(occupancyMatrix)));

% Check variable values
assert(all(all(dilatedMatrix>=occupancyMatrix)));
percentOccupied = fcn_GridMapGen_dilateOccupancyStats(dilatedMatrix,-1);
assert(isequal(round(percentOccupied,6),0.7));

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

occupancyMatrix = rand(100,100)<0.1; % about 10 percent occupied
% percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix, (28282));
dilationLevel = 1;

% Call the function
dilatedMatrix = fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, ([]));

% Check variable types
assert(isnumeric(dilatedMatrix));

% Check variable sizes
assert(isequal(size(dilatedMatrix),size(occupancyMatrix)));

% Check variable values
assert(all(all(dilatedMatrix>=occupancyMatrix)));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

occupancyMatrix = rand(100,100)<0.1; % about 10 percent occupied
% percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix, (28282));
dilationLevel = 1;

% Call the function
dilatedMatrix = fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, (-1));

% Check variable types
assert(isnumeric(dilatedMatrix));

% Check variable sizes
assert(isequal(size(dilatedMatrix),size(occupancyMatrix)));

% Check variable values
assert(all(all(dilatedMatrix>=occupancyMatrix)));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

occupancyMatrix = rand(100,100)<0.1; % about 10 percent occupied
% percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix, (28282));
dilationLevel = 1;
     
Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    dilatedMatrix = fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    dilatedMatrix = fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, (-1));
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