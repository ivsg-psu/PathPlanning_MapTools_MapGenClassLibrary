function h_plot = fcn_MapGen_plotPolytopes(polytopes, varargin)
%fcn_MapGen_plotPolytopes    plots XY data with user-defined formatting strings
% 
% FORMAT:
%
%      h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat), (fillFormat), (fig_num))
%
% INPUTS:  
%
%      polytopes: an array of structures containing polytope information.
%      
%      (OPTIONAL INPUTS)
%
%      plotFormat: one of the following:
%      
%          * a format string, e.g. 'b-', that dictates the plot style
%          * a [1x3] color vector specifying the RGB ratios from 0 to 1
%          * a structure whose subfields for the plot properties to change, for example:
%            plotFormat.LineWidth = 3;
%            plotFormat.MarkerSize = 10;
%            plotFormat.Color = [1 0.5 0.5];
%            A full list of properties can be found by examining the plot
%            handle, for example: h_plot = plot(1:10); get(h_plot)
%
%      fillFormat: a 1-by-5 vector to specify wether or not there is fill,
%      the color of fill, and the opacity of the fill [Y/N, R, G, B,
%      alpha]. Default is [0 0 0 0 0];
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      h_plot: the handle to the plotting result
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
% 
%       script_test_fcn_MapGen_plotPolytopes.m 
%  
%       for a full test suite.
%
% This function was written on 2021_06_06 by Sean Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history
% 2021_06_06 - S.Brennan
% -- revised from fcn_plot_polytopes
% -- added revisions to prep for MapGen library
% 2021_06_07 - S.Brennan
% -- updated examples in header
% -- added test script for function
% 2023_01_15 - S.Brennan
% -- uses narginchk now
% 2023_02_20 - S.Brennan
% -- checks if figure exists
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_07_16 - Sean Brennan
% -- imported plot structure into format from fcn_plotRoad_plotXY
% 2025_07_17 by Sean Brennan
% -- standardized Debugging and Input checks area, Inputs area
% -- made codes use MAX_NARGIN definition at top of code, narginchk
% -- made plotting flag_do_plots and code consistent across all functions

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 4; % The largest Number of argument inputs to the function
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

        % Check the polytopes input
        % fcn_DebugTools_checkInputsToFunctions(polytopes, 'polytopes');

    end
end

% Set plotting defaults
plotFormat = 'k';

% formatting_type type is an integer to save the type of formatting. 
% The numbers refer to 1: a string is given or 2: a color is given, or 3: a structure is given
formatting_type = 1;  

% Check to see if user specifies plotFormat?
if 2 <= nargin
    input = varargin{1};
    if ~isempty(input)
        plotFormat = input;
        if ischar(plotFormat) && length(plotFormat)<=4
            formatting_type = 1;
        elseif isnumeric(plotFormat)  % Numbers are a color style
            formatting_type = 2;
        elseif isstruct(plotFormat)  % Structures give properties
            formatting_type = 3;
        else
            warning('on','backtrace');
            warning('An unkown input format is detected - throwing an error.')
            error('Unknown plotFormat input detected')
        end
    end
end

% Set fillFormat defaults
fillFormat = [0 0 0 0 0]; % initialize empty values
if 2 <= nargin
    temp = varargin{2};
    if ~isempty(temp)
        fillFormat = temp;
    end
end

% Does user want to show the plots?
flag_do_plots = 1; % Default is to show plots
fig_num = []; % Empty by default
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
end

if isempty(fig_num)
    temp = gcf;
    fig_num = temp.Number;
end


%% Solve for the Maxs and Mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the output
h_plot = 0;

% Fill in the x and y data
polytope_plot_data_x = [];
polytope_plot_data_y = [];
for polys = 1:size(polytopes,2) % plot each polytope
    polytope_plot_data_x = [polytope_plot_data_x; polytopes(polys).vertices(:,1); nan]; %#ok<AGROW>
    polytope_plot_data_y = [polytope_plot_data_y; polytopes(polys).vertices(:,2); nan]; %#ok<AGROW>
end

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
    % check whether the figure already has data
    temp_h = figure(fig_num);
    flag_rescale_axis = 0; 
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 0; % Set to 1 to force rescaling
    end        

    hold on;
    axis equal

    xlabel('X [m]');
    ylabel('Y [m]');

    % make plots
    if formatting_type==1
        finalPlotFormat = fcn_DebugTools_extractPlotFormatFromString(plotFormat, (-1));
    elseif formatting_type==2
        finalPlotFormat.Color = plotFormat;
    elseif formatting_type==3        
        finalPlotFormat = plotFormat;
    else
        warning('on','backtrace');
        warning('An unkown input format is detected in the main code - throwing an error.')
        error('Unknown plot type')
    end

    % If plotting only one point, make sure point style is filled
    NplotPoints = length(polytope_plot_data_x(:,1));
    if NplotPoints==1
        if ~isfield(plotFormat,'Marker') || strcmp(plotFormat.Marker,'none')
            finalPlotFormat.Marker = '.';
            finalPlotFormat.LineStyle = 'none';
        end
    end


    % Do plotting
    % Plot polytope as filled object, using 'fill'
    if fillFormat(1) == 1
        for polys = 1:size(polytopes,2)
            filler = fill(polytopes(polys).vertices(:,1)',polytopes(polys).vertices(:,2)',fillFormat(2:4));
            filler.FaceAlpha = polytopes(polys).cost;
        end
    end


    % Plot polytope edges depending on line style
    h_plot = plot(polytope_plot_data_x,polytope_plot_data_y);    
    list_fieldNames = fieldnames(finalPlotFormat);
    for ith_field = 1:length(list_fieldNames)
        thisField = list_fieldNames{ith_field};
        h_plot.(thisField) = finalPlotFormat.(thisField);
    end

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

    axis equal
    axis tight

end % Ends check if plotting

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function

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



