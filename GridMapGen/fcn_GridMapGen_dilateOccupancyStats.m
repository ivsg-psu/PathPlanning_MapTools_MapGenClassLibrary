function percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix, varargin)
% fcn_GridMapGen_dilateOccupancyStats  calculates percentage occupancy,
% namely the percentage of pixels greater than or equal to 1 (occupied). 
% 
% FORMAT:
% 
%     percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix), (fig_num))
% 
% INPUTS:
% 
%     occupancyMatrix: N-by-M matrix
% 
%     (optional inputs)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
% 
% OUTPUTS:
% 
%     percentOccupied: a 1x1 scalar between 0 and 1 representing percentage
%     of elements in the occupancyMatrix that are greater than or equal to
%     1
% 
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
% 
% EXAMPLES:
%
%      % Example: Generates about 50% occupied
%      occupancyMatrix = rand(100,100);
%      percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix)
%
% See the script: script_test_fcn_GridMapGen_dilateOccupancyStats
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
MAX_NARGIN = 2; % The largest Number of argument inputs to the function
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
    debug_fig_num = 999978; %#ok<NASGU>
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
        narginchk(1,MAX_NARGIN);

        % Check the occupancyMatrix input, must be [2+ 2+] in size
        fcn_DebugTools_checkInputsToFunctions(occupancyMatrix*1.0, 'positive_2orMorecolumn_of_numbers',[2 3]);

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

% Both methods are about the same, in timing. See notes below
if 1==0
    % Takes 0.056462 milliseconds per run, on average of 100 repeats
    totalElements = numel(occupancyMatrix);
    occupiedElements = sum(occupancyMatrix>=1,'all');
    percentOccupied = occupiedElements/totalElements;
else
    % Takes 0.056518 milliseconds per run, on average of 100 repeats
    [n,m] = size(occupancyMatrix);
    percentOccupied = sum(sum(occupancyMatrix>=1))/(n*m);
end


% If all are 1's and 0's, can use this instead (faster)
% percentOccupied = sum(sum(matrix))/(n*m);

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
    fprintf(1,'Percent occupancy: %.2f %%\n',percentOccupied*100);

    figure(fig_num);

    % Display dialated on top of original by making testOccupancy
    % pixels have value of 3, dilated have values of 2, 1 otherwise
    % Then use a colormap to map 3 values to white, 2 to grey, 1 to
    % white
    image((occupancyMatrix>0) + 1);
    colormap([1 1 1; 0 0 0])
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




