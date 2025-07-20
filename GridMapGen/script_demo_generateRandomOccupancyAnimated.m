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
Nsteps = 10;
occupancyRatio = 0.2;
dilationLevel = 40;
seedMap = rand(nRows,mColumns);
Nseeds = numel(seedMap);

fig_num = 1111;
for ith_step = 1:300 %:Nsteps
    
    % Resample top Nrand values
    Nrand = 20;
    seedVector = reshape(seedMap,[],1);
    [~,sortedRandInd] = sort(seedVector,'descend');
    seedVector(sortedRandInd(1:Nrand),1) = rand(Nrand,1);
    seedMap = reshape(seedVector,nRows,mColumns);

    % Resample Nrand values
    Nrand = 200;
    randomThreshold = Nrand/Nseeds;
    randomChange = rand(nRows,mColumns);
    indicesToChange = find(randomChange<randomThreshold);
    seedMap(indicesToChange) = rand(length(indicesToChange),1);

    % Randomly walk sideways
    percentageSideways = 1;
    randomChange = rand(nRows,mColumns);
    indicesChange = find(randomChange<percentageSideways);
    % Do not walk last column
    indicesChange = indicesChange(indicesChange<(nRows*(mColumns-1)));
    seedMap(indicesChange+nRows) = seedMap(indicesChange);
    firstColumnChanged = find(indicesChange<=nRows);
    seedMap(firstColumnChanged) = rand(length(firstColumnChanged),1);

    

    % Call the function, plotting results in fig_num
    fcn_GridMapGen_generateRandomOccupancyMap(...
        (nRows), ...
        (mColumns), ...
        (occupancyRatio), ...
        (dilationLevel), ...
        (seedMap),... %(1.0*(seedMap<0.2)), ...
        (fig_num));

    pause(0.1);
end
