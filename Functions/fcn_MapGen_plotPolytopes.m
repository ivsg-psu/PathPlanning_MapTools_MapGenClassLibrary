function [fig] = fcn_MapGen_plotPolytopes(polytopes,fig_num,line_spec,line_width,varargin)
% fcn_MapGen_plotPolytopes plots a set of polytopes given user-defined
% inputs
%
% [FIG]=fcn_MapGen_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH)
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
%      halton_range = [1 100]; % Sets the indices to use from halton set
%      [polytopes] = fcn_MapGen_haltonVoronoiTiling(halton_range);
%      fig1 = fcn_MapGen_plotPolytopes(polytopes,[],'-',2,[0.5 0 0]);
%      fig2 = fcn_MapGen_plotPolytopes(polytopes,998,'-',2,[0 0 0.5],[0 mapx 0 mapy]);
%      fig3 = fcn_MapGen_plotPolytopes(polytopes,999,'-',2,[0 0.5 0],[0 mapx 0 mapy],'square');
%      fig4 = fcn_MapGen_plotPolytopes(polytopes,1000,'-',2,[0 0 0],[0 mapx 0 mapy],'square',[1 0 0 0 0.5]);
%      halton_range = [101 200]; % Sets the indices to use from halton set
%      [polytopes2] = fcn_MapGen_haltonVoronoiTiling(halton_range);
%      fcn_MapGen_plotPolytopes(polytopes2,fig1,'r-',2)
%      fcn_MapGen_plotPolytopes(polytopes2,fig2,'b--',2,[0 mapx 0 mapy])
%      fcn_MapGen_plotPolytopes(polytopes2,fig3,'g-',3,[0 mapx 0 mapy],'square')
%      fcn_MapGen_plotPolytopes(polytopes2,fig4,'k-',3,[0 mapx 0 mapy],'square',[1 0 0 0 0.5])
%
% For more examples, see: script_test_fcn_MapGen_plotPolytopes.m
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

%% chec input arguments
if nargin < 4 || nargin > 8
    error('Incorrect number of input arguments')
end

%% open figures
is_axes = 0;
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
        elseif arg_size == 3
            color = argument; % color to plot polytopes
            plots = 2; % specific color plot
        elseif arg_size == 4
            axis_limits = argument; % limits of x and y axes
        elseif arg_size == 5
            fill_info = argument; % all the fill information [Y/N R G B alpha]
        else
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
if fill_info(1) == 1
    for polys = 1:size(polytopes,2)
        filler = fill(polytopes(polys).vertices(:,1)',polytopes(polys).vertices(:,2)',fill_info(2:4));
        % use polytope shading darkness based on polytope traversal cost
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

% Plot depending on line style
if plots == 1 % basic plot
    plot(polytope_plot_data_x,polytope_plot_data_y,line_spec,'linewidth',line_width)
else
    plot(polytope_plot_data_x,polytope_plot_data_y,line_spec,'Color',color,'linewidth',line_width)
end


axis(axis_limits);
axis(axis_style);





%{

function [map_polytopes,all_pts,mu_rad_final,sigma_rad_final] = ...
    fcn_MapGen_polytopeMapGen(...
    low_pt,high_pt,...
    xlow,xhigh,ylow,yhigh,...
    des_radius,sigma_radius,min_rad,shrink_seed)
% fcn_MapGen_polytopeMapGen creates a polytope map generated by a Halton
% set polytopes, with given desired radius and standard deviation in
% radius.
%
% FORMAT:
%
% [map_polytopes,all_pts,mu_rad_final,sigma_rad_final] = ...
%     fcn_MapGen_polytopeMapGen(...
%     low_pt,high_pt,...
%     xlow,xhigh,ylow,yhigh,...
%     des_radius,sigma_radius,min_rad,shrink_seed)
%
% INPUTS:
%
%    Halton_range: 1 x 2 vector of [low high] range of Halton points to use
%    to generate the tiling
%
%    xlow:
%
%    xhigh:
%
%    ylow:
%
%    yhigh:
%
%    des_radius:
%
%    sigma_radius:
%
%    min_rad:
%
%    shrink_seed:
%
%    (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%     map_polytopes: an N x 2 matrix representing the [x y] vector of starting
%     points of the "walls", where N is # of wall segments
%
%     all_pts: an N x 2 matrix representing the [x y] vector of ending
%     points of the "walls", where N is # of wall segments
%
%     mu_rad_final:
%
%     sigma_rad_final
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_polytopeMapGen
% for a full test suite.
%
% This function was written on 2020_06_06 by S. Brennan by editing
% fcn_basic_polytope_map_generation written by S. Tau.
% Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
% 2021_06_06
% -- first written by S. Brennan.

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

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
%     % Are there the right number of inputs?
%     if nargin < 0 || nargin > 1
%         error('Incorrect number of input arguments')
%     end

    %     % Check the wall_start input
    %     fcn_vis_checkInputsToFunctions(...
    %         wall_start, '2column_of_numbers');
    %
    %     % Use wall_start to calculate the number of walls
    %     Nwalls = length(wall_start(:,1));
    %
    %     % Check the wall_end input
    %     fcn_vis_checkInputsToFunctions(...
    %         wall_end, '2column_of_numbers',Nwalls);
    %
    %     % Check the sensor_location input
    %     fcn_vis_checkInputsToFunctions(...
    %         sensor_location, '2column_of_numbers',[1 1]);
end

%
% % Does user want to show the plots?
% if 1 == nargin
%     fig_num = varargin{end};
%     figure(fig_num);
%     flag_do_plot = 1;
% else
%     if flag_do_debug
%         fig = figure;
%         fig_for_debug = fig.Number;
%         flag_do_plot = 1;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Halton_range = [low_pt, high_pt];

% generate Voronoi tiling from Halton points
tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(low_pt,high_pt);
% remove the edge polytope that extend past the high and low points
trim_polytopes = fcn_polytope_editing_remove_edge_polytopes(tiled_polytopes,xlow,xhigh,ylow,yhigh);
% shink the polytopes so that they are no longer tiled
rng(shrink_seed) % set the random number generator with the shrink seed
[map_polytopes,mu_rad_final,sigma_rad_final] = fcn_MapGen_polytopeShrinkToRadius(trim_polytopes,des_radius,sigma_radius,min_rad);

% gather data on all the points
point_tot = length([map_polytopes.xv]); % total number of vertices in the polytopes
beg_end = zeros(1,point_tot); % is the point the start/end of an obstacle
curpt = 0;
for poly = 1:size(map_polytopes,2) % check each polytope
    verts = length(map_polytopes(poly).xv);
    map_polytopes(poly).obs_id = ones(1,verts)*poly; % obs_id is the same for every vertex on a single polytope
    beg_end([curpt+1,curpt+verts]) = 1; % the first and last vertices are marked with 1 and all others are 0
    curpt = curpt+verts;
end
obs_id = [map_polytopes.obs_id];
all_pts = [[map_polytopes.xv];[map_polytopes.yv];1:point_tot;obs_id;beg_end]'; % all points [x y point_id obs_id beg_end]

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

% if flag_do_plot
%     fcn_Vis_plotWalls(...
%         wall_start,wall_end,...
%         fig_num);
% end
%
% if flag_do_debug
%     fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
% end


end % Ends the function
%}

