function dilatedMatrix = fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, varargin)
% fcn_GridMapGen_dilateByN  dilates a matrix by N cells
% 
% given an N-by-M occupancyMatrix, returns an N-by-M dilatedMatrix which is
% occupancyMatrix dialated by n pixels. 
% 
% NOTE: this is NOT true dilation. For true dilation, see the image
% processing toolbox.
% 
% FORMAT:
% 
%     dilatedMatrix = fcn_GridMapGen_dilateByN(occupancyMatrix, n, (fig_num))
% 
% INPUTS:
% 
%     occupancyMatrix: N-by-M matrix
%   
%     dilationLevel: number of cells (or pixels) to dilate by
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
%     dilatedMatrix: a matrix where all elements of occupancyMatrix>=1 have
%     been dilated outward by n pixels
% 
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
% 
% EXAMPLES:
%
%     occupancyMatrix = rand(100,100)<0.1; % about 10 percent occupied
%     percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix*1.0, (28282));
%     dilationLevel = 1;
%     
%     % Call the function
%     dilatedMatrix = fcn_GridMapGen_dilateByN(occupancyMatrix, dilationLevel, (fig_num));
%
% See the script: script_test_fcn_GridMapGen_dilateByN
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
MAX_NARGIN = 3; % The largest Number of argument inputs to the function
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
        narginchk(2,MAX_NARGIN);

        % Check the occupancyMatrix input, must be [2+ 2+] in size
        fcn_DebugTools_checkInputsToFunctions(occupancyMatrix*1.0, 'positive_2orMorecolumn_of_numbers',[2 3]);

        % Check the dilationLevel input, must be [1 1] integer
        fcn_DebugTools_checkInputsToFunctions(dilationLevel, 'positive_1column_of_integers',[1 1]);
        
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
[n,m] = size(occupancyMatrix);
dilatedMatrix = occupancyMatrix;

multiplier1 = diag(ones(n,1),0) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
multiplier2 = diag(ones(m,1),0) + diag(ones(m-1,1),1) + diag(ones(m-1,1),-1);

for i=1:dilationLevel
    dilatedMatrix = multiplier1*dilatedMatrix*multiplier2;

    % The following code makes the uppper/lower/right/left walls all
    % "occupied"
    if 1==0
        dilatedMatrix(:,1) = max(max(dilatedMatrix));
        dilatedMatrix(:,end) = max(max(dilatedMatrix));
        dilatedMatrix(1,:) = max(max(dilatedMatrix));
        dilatedMatrix(end,:) = max(max(dilatedMatrix));
    end
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
    figure(fig_num);

    % Display dialated on top of original by making testOccupancy
    % pixels have value of 3, dilated have values of 2, 1 otherwise
    % Then use a colormap to map 3 values to white, 2 to grey, 1 to
    % white
    image(2*(occupancyMatrix>0) + 1.0*(dilatedMatrix>0) + 1);
    colormap([1 1 1;0.5 0.5 0.5;0 0 0])
    
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




