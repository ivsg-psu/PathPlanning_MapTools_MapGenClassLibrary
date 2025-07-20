% script_demo_dilateCompareDilateByNSpeeds.m
% Example script to compare the timing between the different dilation
% functions

% REVISION HISTORY
% 2008 - S. Brennan
% -- first draft for ScriptTest_TimeCompareDialationFunctions
% 2025_07_20 - By S. Brennan
% -- Modified for new GridMapGen functions

%clear;
clc;

%% Create random map

nRows=200;
mColumns=300;
occupancyRatio = 0.15;
dilationLevel = 60;
seedMap = 999;
map = fcn_GridMapGen_generateRandomOccupancyMap(...
    (nRows), ...
    (mColumns), ...
    (occupancyRatio), ...
    (dilationLevel), ...
    (seedMap), ...
    (3737));


% Call the function
fig_num = 1000;
figure(fig_num); clf;
fcn_GridMapGen_dilateByN(map, dilationLevel, [], [], (fig_num));

fig_num = 2000;
figure(fig_num); clf;
fcn_GridMapGen_dilateOccupancyByN(map, dilationLevel, ([]), (fig_num));

%% FAST mode comparisons with regular vs precalculated modes
fig_num = 80001;
fprintf(1,'Figure: %.0f: FAST mode comparisons with regular vs precalculated modes\n',fig_num);
figure(fig_num);
close(fig_num);

occupancyMatrix = rand(100,80)<0.1; % about 10 percent occupied
% percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix, (28282));
dilationLevel = 1;
     
Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [dilatedMatrix, leftDilationMultiplier, rightDilationMultiplier] = ...
    fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, ...
    ([]), ([]), ([]));
end
slow_method = toc;

% Do calculation without pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [dilatedMatrix, leftDilationMultiplier, rightDilationMultiplier] = ...
        fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, ...
        ([]), ([]), (-1));
end
fast_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [dilatedMatrix, leftDilationMultiplier, rightDilationMultiplier] = ...
        fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, ...
        (leftDilationMultiplier), (rightDilationMultiplier), (-1));
end
very_fast_method = toc;

% Do calculation without pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [dilatedMatrix, dilationIndexShift] = ...
        fcn_GridMapGen_dilateOccupancyByN(occupancyMatrix, dilationLevel, ...
        ([]), (-1));
end
fast_method_Occupancy = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [dilatedMatrix, dilationIndexShift] = ...
        fcn_GridMapGen_dilateOccupancyByN(occupancyMatrix, dilationLevel, ...
        (dilationIndexShift), (-1));
end
very_fast_method_Occupancy = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

% Plot results as bar chart
figure(373737);
clf;
hold on;

names = {'Dialate Normal mode','Dialate Fast mode','Dialate Fast+Precalc mode','Occupancy Fast mode','Occupancy fast+Precalc'};
X = categorical(names);
X = reordercats(X,names); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method very_fast_method fast_method_Occupancy very_fast_method_Occupancy]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')

% Plot results as bar chart
figure(23434);
clf;
hold on;

names = {'Dialate Fast+Precalc mode','Occupancy fast+Precalc'};
X = categorical(names);
X = reordercats(X,names); % Forces bars to appear in this exact order, not alphabetized
Y = [very_fast_method very_fast_method_Occupancy]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')

