% script_demo_dilateCompareDilateByNSpeeds.m
% Example script to compare the timing between the different dilation
% functions

% REVISION HISTORY
% 2025_07_20 - S. Brennan
% -- first draft testing random variation propogation



%% Create random map
% Set input arguments
nRows = 100;
mColumns = 100;
mapSize = [nRows mColumns];
Nsteps = 10;
occupancyRatio = 0.2;
dilationLevel = 200;
seedMap = rand(nRows,mColumns);
Nseeds = numel(seedMap);
leftDilationMultiplier = [];
rightDilationMultiplier = [];
optimizedThreshold = [];

% Call the function once to initialize settings for upcoming calls
[occupancyMatrix, randomMatrixDilated, forcedThreshold, leftDilationMultiplier, rightDilationMultiplier] = ...
    fcn_GridMapGen_generateRandomOccupancyMap(...
    'mapSize', (mapSize),... % [nRows mCols])
    'occupancyRatio',(occupancyRatio),... % [1x1] value between 0 and 1
    'dilationLevel',(dilationLevel),.... % [1x1] strictly positive int
    'seedMap', (seedMap),... % [1x1] integer to be a random seed or NxM matrix of random numbers
    'leftDilationMultiplier', (leftDilationMultiplier),... %  [nRows nRows], ...
    'rightDilationMultiplier', (rightDilationMultiplier),... % [mCols mCols], ...
    'thresholdForced', (optimizedThreshold), ... % [1x1] scalar
    'flagSkipThresholdOptimization',(0),...% [1x1] scalar
    'figNum',(-1));

colorMin = min(randomMatrixDilated,[],"all");
colorMax = max(randomMatrixDilated,[],"all");

fig_num = 1111;
figure(fig_num); clf;
numColors = 256;
cmap = turbo(numColors);
colormap(cmap);

h_fig = figure(fig_num);
set(h_fig,'Name','animatedRandom','NumberTitle','off'); %, 'Position',[684 85 592 317]);

for ith_step = 1:50 %:Nsteps
    
    % % Resample top Nrand values
    % % Increasing this number causes objects to "disappear" more as they
    % % progress in time
    % Nrand = 20;
    % seedVector = reshape(seedMap,[],1);
    % [~,sortedRandInd] = sort(seedVector,'descend');
    % seedVector(sortedRandInd(1:Nrand),1) = rand(Nrand,1);
    % seedMap = reshape(seedVector,nRows,mColumns);

    %%%%%%%%%%%
    % Change the map slightly to "evolve" the seeds
    % Resample Nrand values
    Nrand = 20;
    randomThreshold = Nrand/Nseeds;
    randomChange = rand(nRows,mColumns);
    indicesToChange = find(randomChange<randomThreshold);
    seedMap(indicesToChange) = rand(length(indicesToChange),1);

    %%%%%%%%%%%
    % Randomly walk sideways
    movementSideways = 0.6; %.5; %2.3;
    percentageSideways = mod(movementSideways,1); % A value between 0 and 1
    columnsSideways = floor(movementSideways);
    
    % Move the percentage
    if percentageSideways>0
        % Do not walk last columns, and refill first columns
        randomChange = rand(nRows,mColumns);
        indicesChange = find(randomChange<percentageSideways);
        indicesChange = indicesChange(indicesChange<(nRows*(mColumns-1)));
        seedMap(indicesChange+nRows) = seedMap(indicesChange);
        firstColumnChanged = find(indicesChange<=nRows);
        seedMap(firstColumnChanged) = rand(length(firstColumnChanged),1);
    end

    % Move the columns
    if columnsSideways>0
        randomChange = rand(nRows,columnsSideways);
        seedMap = [randomChange seedMap(:,1:(mColumns-columnsSideways))];
    end


    %%%%%%%%%%%
    % Update the map based on changed seed
    % Call the function once to initialize settings for upcoming calls
    [~, randomMatrixDilated, updatedThreshold, leftDilationMultiplier, rightDilationMultiplier] = ...
        fcn_GridMapGen_generateRandomOccupancyMap(...
        'seedMap', (seedMap),... % [1x1] integer to be a random seed or NxM matrix of random numbers
        'leftDilationMultiplier', (leftDilationMultiplier),... %  [nRows nRows], ...
        'rightDilationMultiplier', (rightDilationMultiplier),... % [mCols mCols], ...
        'thresholdForced', (forcedThreshold), ... % [1x1] scalar
        'flagSkipThresholdOptimization',(0),...% [1x1] scalar
        'figNum',(-1));
    % 'mapSize', (mapSize),... % [nRows mCols])
    % 'occupancyRatio',(occupancyRatio),... % [1x1] value between 0 and 1
    forcedThreshold = 0.9*forcedThreshold + 0.1*updatedThreshold;

    %%%%
    % Rescale the colors to integer values

    % Filter the color minimum and max values
    keep = 0.8;
    colorMin = (keep)*colorMin + (1-keep)*min(randomMatrixDilated,[],"all");
    colorMax = (keep)*colorMax + (1-keep)*max(randomMatrixDilated,[],"all");

    rescaledToColorIntegers = fcn_INTERNAL_rescaleToColors(randomMatrixDilated, colorMin, colorMax, numColors);
    clf;
    image(rescaledToColorIntegers);
    hold on;
    [M,c] = contour(randomMatrixDilated,50,'k-','Linewidth',0.2);

    if 1==0
        % Trying to use a different color map for the contours
        %% See:
        % https://www.mathworks.com/matlabcentral/answers/194554-how-can-i-use-and-display-two-different-colormaps-on-the-same-figure
        colorizedImage = floor(rescale(randomMatrixDilated,1,numColors));
        clf;
        ax1 = axes;
        image(colorizedImage);
        % colormap(ax1, cmap);

        ax2 = axes;
        [M,c] = contour(randomMatrixDilated,50,'Linewidth',0.2);

        %%Link them together
        linkaxes([ax1,ax2])

        %%Hide the top axes
        ax2.Visible = 'off';
        ax2.XTick = [];
        ax2.YTick = [];

        %%Give each one its own colormap
        colormap(ax1,cmap)
        colormap(ax2,'pink')

    end
    pause(0.01);
end

%% fcn_INTERNAL_rescaleToColors
function rescaledToColorIntegers = fcn_INTERNAL_rescaleToColors(randomMatrixDilated, colorMin, colorMax, numColors)

colorInterval = (colorMax-colorMin)/(numColors-1); % The interval represents the "jump" between different colors
offsetRemovedMatrix = randomMatrixDilated - colorMin; % Remove the offset
countMatrix = floor(offsetRemovedMatrix./colorInterval)+1;
rescaledToColorIntegers = min(max(countMatrix,1),numColors);

end % Ends fcn_INTERNAL_rescaleToColors