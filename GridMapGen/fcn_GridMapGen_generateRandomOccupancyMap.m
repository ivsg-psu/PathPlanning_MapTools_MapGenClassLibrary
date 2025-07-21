function [occupancyMatrix, randomMatrixDilated, optimizedThreshold, leftDilationMultiplier, rightDilationMultiplier]  = ...
    fcn_GridMapGen_generateRandomOccupancyMap(varargin)
% fcn_GridMapGen_generateRandomOccupancyMap  generates a random occupancy map
% 
% an occupancy map is generated that consists of 0's where there is open
% space, 1 otherwise. The method to generate the random map is to dilate a
% random matrix, then iterate through different threshold ranges until one
% is found that is close to the desired occupancyRatio. The threshold is
% chosen such that the occupancyMatrix is set to 1 for values of the
% randomMatrixDilated larger than the threshold.
%
% FORMAT:
% 
%     [occupancyMatrix, randomMatrixDilated, optimizedThreshold, leftDilationMultiplier, rightDilationMultiplier]  = ...
%       fcn_GridMapGen_generateRandomOccupancyMap(...
%       'mapSize', (mapSize),... % [nRows mCols])
%       'occupancyRatio',(occupancyRatio),... % [1x1] value between 0 and 1
%       'dilationLevel',(dilationLevel),.... % [1x1] strictly positive int
%       'seedMap', (seedMap),... % [1x1] integer to be a random seed or NxM matrix of random numbers
%       'leftDilationMultiplier', (leftDilationMultiplier),... %  [nRows nRows], ...
%       'rightDilationMultiplier', (rightDilationMultiplier),... % [mCols mCols], ...
%       'thresholdforced', (thresholdforced), ... % [1x1] scalar
%       'flagSkipThresholdOptimization',(flagSkipThresholdOptimization),...% [1x1] scalar
%       'figNum',(fig_num));
% 
% INPUTS:
% 
%     'mapSize', [nRows mCols]: lets user specify the number of rows in the
%     output matrix. Defaults to 100x100 if no input given.
%
%     'occupancyRatio',(1x1 value between 0 and 1): allows the user to specify the
%     ratio of occupied elements to total elements. For example, if R =
%     0.5, then 50% of the elements will be occupied. The default value is
%     0.2, which is approximately 20% occupied.
%
%     'dilationLevel',(1x1 integer): allows the user to specify the
%     dilation level (integer), which in turn affects the smoothness.
%     Larger numbers cause more smoothness. Dilation levels of 0 are the
%     same as the random plot. A value of 1, 2, or 3 usually gives good
%     maps, though larger values can be used. The default is 2.
%
%     'seedMap', (1x1 random seed or NxM matrix): allows the user to
%     specify a seed map as an N-by-M matrix, or the user can enter an
%     integer which is then used as a seed value for the random number
%     generator. This is especially useful to compare the same map data at
%     different resolutions or levels of smoothing. Note: if a seedMap is
%     given as a matrix, the size of the seed map overwrites the values of
%     N and M in the first two arguments.
%
%    'leftDilationMultiplier',  [nRows nRows]; See fcn_GridMapGen_dilateByN
%    'rightDilationMultiplier', [mCols mCols]; See fcn_GridMapGen_dilateByN
%
%    'thresholdForced', (1x1 scalar): allows user to give a user-defined
%    threshold for occupancy calculation, instead of using the optimized
%    value. This is rarely used, but useful for animations of occupancy
%    where smoothness is desired between maps that are similar but not
%    exactly identical. Default is empty, to use the optimized value.
%
%    'flagSkipThresholdOptimization', (1x scalar) where, if set to 1, skips
%    the optimization step and uses only the thresholdForced. Note:
%    thresholdForced MUST be non-empty for this flag to work. Default is to
%    set to 0, to cause an optimization of the threshold. This flag is
%    rarely used, but can significantly speed up calculations of animations
%    when the random maps are very similar to each other.
%
%    'figNum',(1x1 integer): allows user to specify a figure number to plot
%    results. If set to -1, skips any input checking or debugging, no
%    figures will be generated, and sets up code to maximize speed. As
%    well, if given, this forces the variable types to be displayed as
%    output and as well makes the input check process verbose.
% 
% OUTPUTS:
% 
%     occupancyMatrix: an NxM matrix filled with random occupancy as
%     logical values, or if user gives numeric seed, numeric values of 0 or
%     1.
%
%     randomMatrixDilated: the NxM matrix of scalar values that were
%     thresholded to generate the occupancy matrix
%
%     optimizedThreshold: the 1x1 threshold that, when applied to the
%     randomMatrixDialated, produces occupancyRatio percentage of values
%     greater than the threshold
% 
%     leftDilationMultiplier, the outputs from the
%     dilation function (see fcn_GridMapGen_dilateByN)
% 
%     rightDilationMultiplier: the outputs from the
%     dilation function (see fcn_GridMapGen_dilateByN)
% 
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_GridMapGen_dilateOccupancyStats
%     fcn_GridMapGen_dilateByN
% 
% EXAMPLES:
%
% See the script: script_test_fcn_GridMapGen_generateRandomOccupancyMap
% for a full test suite.
% 
% This function was written on 2008_10_18 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

% 
% REVISION HISTORY:
% 
% 2008_10_18 by Sean Brennan
% -- first write of function
% 2025_07_17 by Sean Brennan
% -- imported into MapGen library with updates to formatting
% 2025_07_20 by Sean Brennan
% -- improved input handling and input options. Switched inputs to matched
%    % pair format as the input list was getting way too long.

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 18; % The largest Number of argument inputs to the function
flag_max_speed = 0;
figNumIndex = find(strcmp(varargin,'figNum'));
if (~isempty(figNumIndex) && isequal(varargin{figNumIndex+1},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS");
    MATLABFLAG_MAPGEN_FLAG_DO_DEBUG = getenv("MATLABFLAG_MAPGEN_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_MAPGEN_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_MAPGEN_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; 
end


%% check input arguments?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (0==flag_max_speed)
    if 1 == flag_check_inputs

        % Are there the right number of inputs?
        narginchk(0,MAX_NARGIN);

        % Are the input arguments paired?
        if mod(nargin,2)==1
            error('Unmatched input argument detected. There should be an even number of input arguments');
        end
        
    end
end

% Fill in defaults
mapSize = [100 100];
occupancyRatio = 0.2; % use default, which is approximately 20%
dilationLevel = 2; % use default, 2, which gives fairly smooth results
seedMap = []; % Default is no seed map
flag_GenerateNewSeedMap = 1;
flag_do_plots = 0; % Default is to NOT show plots
leftDilationMultiplierInput = [];
rightDilationMultiplierInput = [];
thresholdForced = [];
flagSkipThresholdOptimization = 0;

for ith_input = 1:(nargin/2)
    inputName = varargin{2*ith_input-1};
    inputValue = varargin{2*ith_input};
    if ~isempty(inputValue)
        switch lower(inputName)
            case {'mapsize'}
                mapSize = inputValue;
                if flag_check_inputs
                    % Check the nRows input
                    fcn_DebugTools_checkInputsToFunctions(mapSize, 'strictlypositive_2column_of_integers',[1 1]);
                end

            case {'occupancyratio'}                
                occupancyRatio = inputValue;
                if flag_check_inputs
                    % Check the occupancyRatio input
                    fcn_DebugTools_checkInputsToFunctions(occupancyRatio, 'strictlypositive_1column_of_numbers',[1 1]);
                end

            case {'dilationlevel'}
                dilationLevel = inputValue;
                if flag_check_inputs
                    % Check the dilationLevel input
                    fcn_DebugTools_checkInputsToFunctions(dilationLevel, 'strictlypositive_1column_of_integers',[1 1]);
                end

            case {'seedmap'}
                seedMap = inputValue;
                if isscalar(seedMap)
                    if flag_check_inputs
                        % Check the scalar seedMap input. Must be positive 1x1
                        % integer
                        fcn_DebugTools_checkInputsToFunctions(seedMap, 'positive_1column_of_numbers',[1 1]);
                    end
                    rng(seedMap);
                else
                    if flag_check_inputs
                        % Check the matrix seedMap input. Must be positive
                        % 2+ by 2+  matrix
                        fcn_DebugTools_checkInputsToFunctions(seedMap, 'positive_2ormorecolumn_of_numbers',[2 3]);
                    end
                    flag_GenerateNewSeedMap = 0;
                    mapSize = size(seedMap);
                end

            case {'leftdilationmultiplier'}
                leftDilationMultiplierInput = inputValue;
            case {'rightdilationmultiplier'}
                rightDilationMultiplierInput = inputValue;
            case {'thresholdforced'}
                thresholdForced = inputValue;
            case {'flagskipthresholdoptimization'}
                flagSkipThresholdOptimization = inputValue;
            case {'fignum'}
                if (0==flag_max_speed)
                    fig_num = inputValue;
                    figure(fig_num);
                    flag_do_plots = 1;
                end
            otherwise
                error('Unrecognized argument input string: %s',inputName);
        end
    end
end



%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%See: http://patorjk.com/software/taag/#p=display&f=Big&t=Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§
nRows = mapSize(1,1); 
mCols = mapSize(1,2); 

if 1==flag_GenerateNewSeedMap
    seedMap = rand(nRows,mCols);
end

% Convert random values or seed input map into occupancyMap, which has only
% 1's or 0's.
initialOccupancyMap = seedMap> (1-occupancyRatio);

% Display occupancyMap (for debugging)?
if flag_do_debug
    h_fig = figure(debug_fig_num);
    set(h_fig,'Name','Original Random','NumberTitle','off', 'Position',[88 326 592 317])
    image(initialOccupancyMap+1);
    colormap([1 1 1;0 0 0])
    percentOccupied = fcn_GridMapGen_dilateOccupancyStats(initialOccupancyMap, -1);
    title(sprintf('Percent occupied at dialation %d, r of %.2f: %.4f\n',dilationLevel, occupancyRatio,percentOccupied));
end

% Dilate on original random image with dialation function
% randomMatrixDilated = fcn_GridMapGen_dilateByN(seedMap, dilationLevel, [], [], -1);
[randomMatrixDilated, leftDilationMultiplier, rightDilationMultiplier] = ...
    fcn_GridMapGen_dilateByN(seedMap, dilationLevel, ...
    (leftDilationMultiplierInput), (rightDilationMultiplierInput), (-1));

% Display randomMatrixDilated (for debugging)?
if flag_do_debug
    h_fig = figure(debug_fig_num+1);
    clf;

    set(h_fig,'Name','randomMatrixDilated','NumberTitle','off', 'Position',[684 326 592 317]);

    numColors = 256;
    cmap = turbo(numColors);
    colorizedImage = floor(rescale(randomMatrixDilated,1,numColors));
    image(colorizedImage);
    colormap(cmap);

    percentOccupied = fcn_GridMapGen_dilateOccupancyStats(randomMatrixDilated, -1);
    title(sprintf('Percent occupied at dialation %d, r of %.2f: %.4f\n',dilationLevel, occupancyRatio,percentOccupied));
end

% Did the user shut off the threshold optimization?
if 1==flagSkipThresholdOptimization  && ~isempty(thresholdForced)
    optimizedThreshold = nan;
else
    % Iterate to find the right value of r that preserves ratio
    low = min(min(randomMatrixDilated));
    high = max(max(randomMatrixDilated));

    % Find 100 levels between high and low value
    levels = low:(high-low)/100:high;

    % Loop through levels, saving results into percentages array
    percentages = zeros(length(levels),1);
    for i=1:length(levels)
        temp_occupancy = randomMatrixDilated>levels(i);
        percentages(i) = fcn_GridMapGen_dilateOccupancyStats(temp_occupancy,-1);
    end

    % Find the index that matches the desired occupancy ratio
    indexMatchingDesiredRatio = find(percentages<=occupancyRatio,1);
    optimizedThreshold = levels(indexMatchingDesiredRatio);
end

if flag_do_debug
    fprintf(1,'Current threshold: %.2f\n', optimizedThreshold);
end

% Did user force the threshold?
if isempty(thresholdForced)
    threshold = optimizedThreshold;
else
    threshold = thresholdForced;
end
    
%% Save final result
occupancyMatrix = randomMatrixDilated>threshold;

% Display occupancyMatrix (for debugging)?
if flag_do_debug
    h_fig = figure(debug_fig_num+2);
    clf;
    set(h_fig,'Name','OccupancyMatrixWithSeeds','NumberTitle','off', 'Position',[88 85 592 317]);
    
    colorValues = 1 + 2*(initialOccupancyMap>0) + (occupancyMatrix>0);
    image(colorValues); % Will call color map's 1st row, 2nd row, and 3rd row respectively
    colormap([1 1 1;0 0 0; 1 0 0])

    percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix, -1);
    title(sprintf('Percent occupied at dialation %d, r of %.2f: %.4f\n',dilationLevel, occupancyRatio,percentOccupied));
end

%§
%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_do_plots
    h_fig = figure(fig_num);
    set(h_fig,'Name','generateRandomOccupancyMap','NumberTitle','off'); %, 'Position',[684 85 592 317]);

    subplot(1,2,1);
    numColors = 256;
    cmap = turbo(numColors);
    colorizedImage = floor(rescale(randomMatrixDilated,1,numColors));
    image(colorizedImage);
    colormap(gca, cmap);

    %percentOccupied = fcn_GridMapGen_dilateOccupancyStats(randomMatrixDilated, -1);
    title('randomMatrixDilated','FontSize',10);

    subplot(1,2,2);
    image(occupancyMatrix+1);
    colormap(gca, [1 1 1;0 0 0])
    
    percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix, -1);
    text(20,80,sprintf([...
        'Percent occupied: %.1f %%\n' ...
        '\t dialation: %.0d,\n' ...
        '\t target occupancyRatio: of %.0f%%'],percentOccupied*100, dilationLevel, occupancyRatio*100),...
        'BackgroundColor',[1 1 1],'FontSize',10);
    title('occupancyMatrix','FontSize',10);
    
end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends the function

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




