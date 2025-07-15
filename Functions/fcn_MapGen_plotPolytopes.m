function [fig] = fcn_MapGen_plotPolytopes(polytopes,fig_num,line_spec,line_width,varargin)
% fcn_MapGen_plotPolytopes plots a set of polytopes given user-defined
% inputs
%
% [FIG_HANDLE]=fcn_MapGen_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH)
% returns:
% a figure with the polytopes plotted as specified
%
% with inputs:
% POLYTOPES: a 1-by-n seven field structure, where n <= number of polytopes
%   with fields:
%   vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%     the number of the individual polytope vertices
%   xv: a 1-by-m vector of vertice x-coordinates
%   yv: a 1-by-m vector of vertice y-coordinates
%   distances: a 1-by-m vector of perimeter distances from one point to the
%     next point, distances(i) = distance from vertices(i) to vertices(i+1)
%   mean: average xy coordinate of the polytope
%   area: area of the polytope
%   max_radius: distance from the mean to the farthest vertex
% FIG_NUM: figure number to plot the values on
% LINE_SPEC: a string, line specifications such as color, 'r', and line or
% point type, '--'
% LINE_WIDTH: width of the line to be plotted
% By default the axes will be determined by the plot function
%
% [FIG]=fcn_MapGen_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,COLOR)
% allows the user to specify the input:
% COLOR: a 1-by-3 vector to specify the RGB plot colors [Red Green Blue],
% where 0 <= Red,Green,Blue <= 1
%
% [FIG]=fcn_MapGen_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,AXIS_LIMITS)
% allows the user to specify the input:
% AXIS_LIMTS: a 1-by-4 vector to specify the start and end of the x and y
% axes, [xstart xend ystart yend]
%
% [FIG]=fcn_MapGen_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,AXIS_STYLE)
% allows the user to specify the input:
% AXIS_STYLE: controls the style of the axis, such as square or equal
%
% [FIG]=fcn_MapGen_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,FILL_INFO)
% allows the user to specify the input:
% FILL_INFO: a 1-by-5 vector to specify wether or not there is fill, the
% color of fill, and the opacity of the fill [Y/N, R, G, B, alpha]
%
% [FIG]=fcn_MapGen_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,COLOR,AXIS_LIMITS,AXIS_STYLE,FILL_INFO)
% allows the user to specify any combination of all four inputs in any
% order after LINE_WIDTH
%
% Examples:
%
%      % BASIC example
%      mapx = 1;
%      mapy = 1;
%      seedGeneratorNames = 'haltonset';
%      seedGeneratorRanges = [1 100];
%      AABBs = [0 0 1 1];
%      mapStretchs = [1 1];
%      [polytopes] = fcn_MapGen_voronoiTiling(...
%      seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
%      seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
%      (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
%      (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
%      (-1));
%      fig1 = fcn_MapGen_plotPolytopes(polytopes,[],'-',2,[0.5 0 0]);
%      fig2 = fcn_MapGen_plotPolytopes(polytopes,998,'-',2,[0 0 0.5],[0 mapx 0 mapy]);
%      fig3 = fcn_MapGen_plotPolytopes(polytopes,999,'-',2,[0 0.5 0],[0 mapx 0 mapy],'square');
%      fig4 = fcn_MapGen_plotPolytopes(polytopes,1000,'-',2,[0 0 0],[0 mapx 0 mapy],'square',[1 0 0 0 0.5]);
%      seedGeneratorNames = 'haltonset';
%      seedGeneratorRanges = [101 200];
%      AABBs = [0 0 1 1];
%      mapStretchs = [1 1];
%      [polytopes2] = fcn_MapGen_voronoiTiling(...
%     seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
%     seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
%     (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
%     (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
%     (-1));
%      fcn_MapGen_plotPolytopes(polytopes2,fig1,'r-',2)
%      fcn_MapGen_plotPolytopes(polytopes2,fig2,'b--',2,[0 mapx 0 mapy])
%      fcn_MapGen_plotPolytopes(polytopes2,fig3,'g-',3,[0 mapx 0 mapy],'square')
%      fcn_MapGen_plotPolytopes(polytopes2,fig4,'k-',3,[0 mapx 0 mapy],'square',[1 0 0 0 0.5])
%
% For examples, see: script_test_fcn_MapGen_plotPolytopes.m
%
% This function was written on 2018_12_10 by Seth Tau
% Questions or comments? sat5340@psu.edu
%


% REVISION HISTORY:
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

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin>=2 && isequal(fig_num,-1))
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
else
    debug_fig_num = []; %#ok<NASGU>
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
        narginchk(4,8);

        % % Check the polytopes input, make sure it is 'polytopes' type
        % fcn_DebugTools_checkInputsToFunctions(...
        %     polytopes, 'polytopes');
        % 
        % 
        % % Check the exp_dist input, make sure it is 'positive_column_of_numbers' type
        % fcn_DebugTools_checkInputsToFunctions(...
        %     exp_dist, 'positive_1column_of_numbers',1);

    end
end



%% open figures
if isempty(fig_num)
    fig = figure; % create new figure with next default index
else
    % check to see that the handle is an axis. If so, use it and don't just
    % go to a new figure
    if isgraphics(fig_num,'axes')
        axes(fig_num);
    else
        fig = figure(fig_num); % open specific figure
    end
end
hold on % allow multiple plot calls

%% determine color and axis values
plots = 1; % basic plot
color = []; axis_limits = []; axis_style = []; fill_info = [0 0 0 0 0]; % initialize empty values
if nargin > 4
    for arg = 1:nargin-4
        argument = varargin{arg};
        arg_size = length(argument);
        if ischar(argument)
            axis_style = argument; % axis style (ie square or equal)
        elseif arg_size == 0
            % Do nothing - user left it blank
        elseif arg_size == 3
            color = argument; % color to plot polytopes
            plots = 2; % specific color plot
        elseif arg_size == 4
            axis_limits = argument; % limits of x and y axes
        elseif arg_size == 5
            fill_info = argument; % all the fill information [Y/N R G B alpha]
        else
            warning('on','backtrace');
            warning('Invalid argument. Argument ignored.')
        end
    end
end


% if nargin == 4
%     plots = 1; % basic plot
% else % nargin > 4
%     if length(varargin{1}) == 3 % first optional argument is color
%         color = varargin{1}; % color to plot polytopes
%         plots = 2; % plot polytopes with special color
%         if nargin == 8 % color, limits, limit style, and fill info specified
%             axis_limits = varargin{2}; % limits of x and y axes
%             axis_style = varargin{3}; % axis style (ie square or equal)
%             fill_info = varargin{4}; % all the fill information [Y/N R G B alpha]
%         elseif nargin == 7 % color, limits, and limit style or fill info
%             axis_limits = varargin{2}; % limits of x and y axes
%             if ischar(varargin{3}) % limit style specified
%                 axis_style = varargin{3}; % axis style (ie square or equal)
%             else % fill info specified
%                 fill_info = varargin{3}; % all the fill information [Y/N R G B alpha]
%             end
%         elseif nargin == 6 % only color and limits specified
%             if length(varargin{2}
%             axis_limits = varargin{2}; % limits of x and y axes
%         end
%     elseif length(varargin{1}) == 4 % first optional argument is axis limits
%         plots = 1; % basic plot
%         axis_limits = varargin{1}; % limits of x and y axes
%         if nargin == 6 % only color and limits specified
%             axis_style = varargin{2}; % axis style (ie square or equal)
%         end
%     else % first optional argument is fill info
%
%     end
% end

%% plot polytopes
% Plot polytope as filled object, using 'fill'
if fill_info(1) == 1
    for polys = 1:size(polytopes,2)
        filler = fill(polytopes(polys).vertices(:,1)',polytopes(polys).vertices(:,2)',fill_info(2:4));
        filler.FaceAlpha = polytopes(polys).cost;
    end
end

% Fill in the x and y data
polytope_plot_data_x = [];
polytope_plot_data_y = [];
for polys = 1:size(polytopes,2) % plot each polytope
    polytope_plot_data_x = [polytope_plot_data_x; polytopes(polys).vertices(:,1); nan]; %#ok<AGROW>
    polytope_plot_data_y = [polytope_plot_data_y; polytopes(polys).vertices(:,2); nan]; %#ok<AGROW>
end

% Plot polytope edges depending on line style
if plots == 1 % basic plot
    plot(polytope_plot_data_x,polytope_plot_data_y,line_spec,'linewidth',line_width)
else
    plot(polytope_plot_data_x,polytope_plot_data_y,line_spec,'Color',color,'linewidth',line_width)
end


axis(axis_limits);
axis(axis_style);






