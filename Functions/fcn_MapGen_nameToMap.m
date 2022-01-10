function [polytopes,fig]=fcn_MapGen_nameToMap(...
    map_name,...
    plot_flag,...
    disp_name,...
    varargin)
% fcn_MapGen_nameToMap generates a map based on
% map_name which specifies the map characteristics in a way that allows map
% to be recreated exactly from the name alone.
%
% FORMAT:
%
% [polytopes,fig]=fcn_MapGen_nameToMap(...
%     map_name,...
%     plot_flag,...
%     disp_name,...
%     (fig_num),...
%     (line_spec),...
%     (line_width),...
%     (color),...
%     (axis_limits),
%     (axis_style),...
%     (fill_info))
%
% INPUTS:
%
%     map_name: string with map characteristic numbers seperated by indicating
%       letters. See details below.
%
%     plot_flag: variable determining whether the map should be plotted (1=Yes)
%
%     disp_name: variable determining whether the map name is displayed on the plot
%       (1=Yes), where the name goes (x,y), what color ([r g b]), and the font
%       size (e.g. 12)
%
%       Example [Yes, origin, red, size 12] = [1, 0 0, 1 0 0, 12]
%
%       If not desired set the first value to anything but 1 and the rest will
%       be ignored
%       Example [No, NA, NA, NA] = [0, 1 2, 3 4 5, 6], or 0, or 999
%
%    ( OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results.
%
%     line_spec: a string, line specifications such as color, 'r', and line or
%     point type, '--'
%
%     line_width: width of the line to be plotted
%
%     color: a 1-by-3 vector to specify the RGB plot colors [Red Green Blue],
%     where 0 <= Red,Green,Blue <= 1
%
%     axis_limits: a 1-by-4 vector to specify the start and end of the x and y
%     axes, [xstart xend ystart yend
%
%     axis_style: controls the style of the axis, such as square or equal
%
%     fill_info: a 1-by-5 vector to specify wether or not there is fill,
%     the color of fill, and the opacity of the fill [Y/N, R, G, B, alpha]
%
%
% OUTPUTS:
%
%     polytopes: a 1-by-n seven field structure, where n <= number of
%       the map polytopes with fields:
%       vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%         the number of the individual polytope vertices
%       xv: a 1-by-m vector of vertice x-coordinates
%       yv: a 1-by-m vector of vertice y-coordinates
%       distances: a 1-by-m vector of perimeter distances from one point to the
%         next point, distances(i) = distance from vertices(i) to vertices(i+1)
%       mean: centroid xy coordinate of the polytope
%       area: area of the polytope
%       max_radius: distance from the mean to the farthest vertex
%
%     fig: variable containing the figure number of the plot, if plot_flag
%     is 1
%
% CHARACTERISTICS:
%
% %% List of Map Characteristics:
% %%% Generation:
%       HST: Halton set voronoi tiling
%           requires: minimum value and maximum value of Halton set
%           example: HST 1 1000
% %%% Trimming
%       SQT: Square trimming
%           requires: low x value, high x value, low y value, high y value
%           example: SQT 0 1 0 1
% %%% Shrinking
%       SMV: Shrink to mean radius with specified variance
%           requires: mean radius, standard deviation of radius, minimum
%           radius to shrink to, shrink seed
%           example: SMV 0.02 0.005 1e-6 1234
%
% DEPENDENCIES:
%
%      fcn_MapGen_checkInputsToFunctions
%      fcn_MapGen_haltonVoronoiTiling
%      fcn_MapGen_polytopeCropEdges
%      fcn_MapGen_polytopesShrinkToRadius
%      fcn_MapGen_plotPolytopes
%
% EXAMPLES:
%
% Basic Example:
%   map_name = "HST 1 100 SQT 0 1 0 1 SMV 0.01 0.001 1e-6 1111";
%   plot_flag = 1; disp_name = 0; fig_num = []; line_style = 'r-';
%   line_width = 2;
%   [polytopes,fig]=fcn_MapGen_nameToMap(map_name,plot_flag,disp_name,fig_num,line_style,line_width);
%
% Advanced Example
%   map_name = "HST 30 450 SQT 0 1 0 1 SMV 0.02 0.005 1e-6 1234";
%   plot_flag = 1; disp_name = [1, 0.05 -0.05, 0.5 0.5 0.5, 10];
%   fig_num = 999; line_style = '-'; line_width = 2; color = [0 0 1];
%   axis_limits = [0 1 -0.1 1]; axis_style = 'square';
%   fill_info = [1 1 0 1 0.5];
%   [polytopes,fig]=fcn_MapGen_nameToMap(map_name,plot_flag,disp_name,fig_num,line_style,line_width,color,axis_limits,axis_style,fill_info);
%
%
% See the script: script_test_fcn_MapGen_nameToMap
% for a full test suite.
%
% This function was written on 2020_07_02 by Seth Tau
% Comments added on 2021_02_23 by Seth Tau
% Questions or comments? sat5340@psu.edu


% Revision History:
% 2021-06-08 - S. Brennan
% -- revised function to prep for MapGen class

% TO DO:
% -- (none)

%% Debugging and Input checks
% set an environment variable on your machine with the getenv function...
% in the Matlab console.  Char array of '1' will be true and '0' will be false.
flag_check_inputs = getenv('ENV_FLAG_CHECK_INPUTS');  % '1' will check input args
flag_do_plot = getenv('ENV_FLAG_DO_PLOT'); % '1' will make plots
flag_do_debug = getenv('ENV_FLAG_DO_DEBUG'); % '1' will enable debugging

% if the char array has length 0, assume the env var isn't set and default to...
% dipslaying more information rather than potentially hiding an issue
if length(flag_check_inputs) == 0
    flag_check_inputs = '1';
end
if length(flag_do_plot) == 0
    flag_do_plot = '1';
end
if length(flag_do_debug) == 0
    flag_do_debug = '1';
end

% convert flag from char string to logical
flag_check_inputs = flag_check_inputs == '1';
flag_do_plot = flag_do_plot == '1';
flag_do_debug = flag_do_debug == '1';

if flag_do_debug
    fig_for_debug = 9993;
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
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

if flag_check_inputs
    % Are there the right number of inputs?
    %  not enough   plotting but not enough   too many
    if (nargin<3)||((nargin>3)&&(nargin<6))||(nargin>10)
        error('Incorrect number of arguments.')
    end

    % Check the map_name input
    if ~isstring(map_name) % input must be a string
        if ischar(map_name) % convert to string if a character
            map_name=convertCharsToStrings(map_name);
        else
            error('map_name must be a string.')
        end
    end

    % Check the plot_flag input
    fcn_MapGen_checkInputsToFunctions(...
        plot_flag, '1column_of_numbers',1);

end


% Does user want to show the plots?
if 4 <= nargin
    fig_num = varargin{1};
    if ~isempty(fig_num)
        figure(fig_num);
        flag_do_plot = 1;
    end

else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
        flag_do_plot = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Split map_name
split_name = split(map_name); % split name string at each space

%% Base Map Generation
if sum(split_name=="HST")>0 % Check for Halton Set Tiling (HST)
    HST_index = find(split_name=="HST"); % index of the HST string
    if length(HST_index)>1 % more than one instance of HST
        error('HST is repeated in map_name.')
    end
    % generate base map based on the values following HST
    halton_range = ...
        [str2double(split_name(HST_index+1)),str2double(split_name(HST_index+2))];
    base_polytopes = fcn_MapGen_haltonVoronoiTiling(halton_range);

% elseif % check for other generating methods

    % generate using other method
else % no map generation method specified
    error('No map generation method specified.')
end

%% Trim map if necessary
if sum(split_name=="SQT")>0 % check for square triming (SQT)
    SQT_index = find(split_name=="SQT"); % index of the SQT string
    if length(SQT_index)>1 % more than one instance of SQT
        error('SQT is repeated in map_name')
    end

    % trim the base polytopes based on the values following SQT
    bounding_box = ...
        [str2double(split_name(SQT_index+1)),str2double(split_name(SQT_index+3));
        str2double(split_name(SQT_index+2)),str2double(split_name(SQT_index+4))];

    trim_polytopes = fcn_MapGen_polytopeCropEdges(...
        base_polytopes,bounding_box);
% elseif % check for other trim methods
    % trim polytopes
else % no trim method specified
    trim_polytopes = base_polytopes; % no trimming so trim_polytopes is the same
end

%% Shrink map if necessary
if sum(split_name=="SMV")>0 % check for shrinking to mean and variance (SMV)
    SMV_index = find(split_name=="SMV"); % index of the SMV string
    if length(SMV_index)>1 % more than 1 instance of SMV
        error('SMV is repeated in map_name')
    end
    % shrink based on the values following SMV
    rng(str2double(split_name(SMV_index+4))) % set the rng to the shrink seed
    polytopes = fcn_MapGen_polytopesShrinkToRadius(...
        trim_polytopes,str2double(...
        split_name(SMV_index+1)),...
        str2double(split_name(SMV_index+2)),...
        str2double(split_name(SMV_index+3)));
% elseif % check for other shrink methods
    % trim polytopes
else % no trim method specified
    polytopes = trim_polytopes; % no shrinking so polytopes is the same
end

fig = []; % set value empty to return as default. Value is filled below if plotting is turned on

%% Plot results?
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

if flag_do_plot

    if length(varargin)>3 % detailed plotting
        [fig] = fcn_MapGen_plotPolytopes(polytopes,varargin{1},varargin{2},varargin{3},varargin{4:end});
    else % basic plotting
        [fig] = fcn_MapGen_plotPolytopes(polytopes,varargin{1},varargin{2},varargin{3});
    end

    % Show the name
    if disp_name(1) == 1 % add map_name to plot
        text(disp_name(2),disp_name(3),map_name,'color',disp_name(4:6),'FontSize',disp_name(7));
    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends the function







