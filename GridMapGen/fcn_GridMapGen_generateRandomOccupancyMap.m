function occupancyMatrix = fcn_GridMapGen_generateRandomOccupancyMap(varargin)
% fcn_GridMapGen_generateRandomOccupancyMap  generates a random occupancy map
% 
% an occupancy map is generated that consists of 0's where there is open
% space, 1 otherwise
%
% FORMAT:
% 
%     occupancyMatrix = fcn_GridMapGen_generateRandomOccupancyMap(...
%       (nRows), ...
%       (mColumns), ...
%       (occupancyRatio), ...
%       (dilationLevel), ...
%       (seedMap), ...
%       (fig_num))
% 
% INPUTS:
% 
%     (optional inputs)
%     nRows: the number of rows in the output matrix. Defaults to 100 if no
%     input given.
%
%     mColumns: the number of columns in the output matrix. Defaults to
%     nRows if no input given.
%
%     occupancyRatio: allows the user to specify the ratio of occupied
%     elements to total elements. For example, if R = 0.5, then 50% of the
%     elements will be occupied. The default value is 0.2, which is
%     approximately 20% occupied.
%
%     dilationLevel: allows the user to specify the dilation level
%     (integer), which in turn affects the smoothness. Larger numbers cause
%     more smoothness. Dilation levels of 0 are the same as the random
%     plot. A value of 1, 2, or 3 usually gives good maps, though larger
%     values can be used. The default is 2.
%
%     seedMap: allows the user to specify a seed map as an N-by-M matrix,
%     or the user can enter an integer which is then used as a seed value
%     for the random number generator. This is especially useful to compare
%     the same map data at different resolutions or levels of smoothing.
%     Note: if a seedMap is given as a matrix, the size of the seed map
%     overwrites the values of N and M in the first two arguments.
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
% 
% OUTPUTS:
% 
%     occupancyMatrix: an NxM matrix filled with random occupancy
% 
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
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

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 6; % The largest Number of argument inputs to the function
flag_max_speed = 0;
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
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

        % % Check the occupancyMatrix input, must be [2+ 2+] in size
        % fcn_DebugTools_checkInputsToFunctions(occupancyMatrix*1.0, 'positive_2orMorecolumn_of_numbers',[2 3]);
        % 
        % % Check the dilationLevel input, must be [1 1] integer
        % fcn_DebugTools_checkInputsToFunctions(dilationLevel, 'positive_1column_of_integers',[1 1]);
        
    end
end

% check variable argument nRows
nRows = 100; % use default nRows
if 1 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        nRows = temp;
        if flag_check_inputs
            % Check the nRows input
            fcn_DebugTools_checkInputsToFunctions(nRows, 'strictlypositive_1column_of_integers',[1 1]);
        end
    end
end

% check variable argument mColumns
mColumns = nRows; % use default, which is to match nRows
if 2 <= nargin
    temp = varargin{2};
    if ~isempty(temp)
        mColumns = temp;
        if flag_check_inputs
            % Check the mColumns input
            fcn_DebugTools_checkInputsToFunctions(mColumns, 'strictlypositive_1column_of_integers',[1 1]);
        end
    end
end

% check variable argument occupancyRatio
occupancyRatio = 0.2; % use default, which is approximately 20%
if 3 <= nargin
    temp = varargin{3};
    if ~isempty(temp)
        occupancyRatio = temp;
        if flag_check_inputs
            % Check the occupancyRatio input
            fcn_DebugTools_checkInputsToFunctions(occupancyRatio, 'strictlypositive_1column_of_numbers',[1 1]);
        end
    end
end

% check variable argument dilationLevel
dilationLevel = 2; % use default, 2, which gives fairly smooth results
if 4 <= nargin
    temp = varargin{4};
    if ~isempty(temp)
        dilationLevel = temp;
        if flag_check_inputs
            % Check the dilationLevel input
            fcn_DebugTools_checkInputsToFunctions(dilationLevel, 'strictlypositive_1column_of_integers',[1 1]);
        end
    end
end

% check variable argument dilationLevel
seedMap = []; % Default is no seed map
flag_GenerateNewSeedMap = 1;
if 5 <= nargin
    temp = varargin{5};
    if ~isempty(temp)
        seedMap = temp;
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
            nRows = size(seedMap,1);
            mColumns = size(seedMap,2);
        end
    end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
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

if 1==flag_GenerateNewSeedMap
    seedMap = rand(nRows,mColumns);
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
randomMap_dilated = fcn_GridMapGen_dilateByN(seedMap, dilationLevel, [], [], -1);


% Display randomMap_dilated (for debugging)?
if flag_do_debug
    h_fig = figure(debug_fig_num+1);
    clf;

    set(h_fig,'Name','randomMap_dilated','NumberTitle','off', 'Position',[684 326 592 317]);

    numColors = 256;
    cmap = turbo(numColors);
    colorizedImage = floor(rescale(randomMap_dilated,1,numColors));
    image(colorizedImage);
    colormap(cmap);

    percentOccupied = fcn_GridMapGen_dilateOccupancyStats(randomMap_dilated, -1);
    title(sprintf('Percent occupied at dialation %d, r of %.2f: %.4f\n',dilationLevel, occupancyRatio,percentOccupied));
end

% Iterate to find the right value of r that preserves ratio
low = min(min(randomMap_dilated));
high = max(max(randomMap_dilated));

% Find 100 levels between high and low value
levels = low:(high-low)/100:high;

% Loop through levels, saving results into percentages array
percentages = zeros(length(levels),1);
for i=1:length(levels)
    temp_occupancy = randomMap_dilated>levels(i);
    percentages(i) = fcn_GridMapGen_dilateOccupancyStats(temp_occupancy,-1);
end

% Find the index that matches the desired occupancy ratio
indexMatchingDesiredRatio = find(percentages<=occupancyRatio,1);
goodThresholdToUse = levels(indexMatchingDesiredRatio);


%% Save final result
occupancyMatrix = randomMap_dilated>goodThresholdToUse;

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
    set(h_fig,'Name','occupancyMatrix','NumberTitle','off', 'Position',[684 85 592 317]);

    image(occupancyMatrix+1);
    colormap([1 1 1;0 0 0])
    
    percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix, -1);
    text(20,80,sprintf([...
        'Percent occupied: %.1f %%\n' ...
        '\t dialation: %.0d,\n' ...
        '\t target occupancyRatio: of %.0f%%'],percentOccupied*100, dilationLevel, occupancyRatio*100),...
        'BackgroundColor',[1 1 1],'FontSize',10);
    
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




