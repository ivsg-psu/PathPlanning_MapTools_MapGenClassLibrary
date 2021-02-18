function [fig] = fcn_plot_polytopes(polytopes,fig_num,line_spec,line_width,varargin)
% FCN_PLOT_POLYTOPES plot the polytopes as specified
%
% [FIG]=FCN_PLOT_POLYTOPES(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH)
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
% [FIG]=FCN_PLOT_POLYTOPES(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,COLOR)
% allows the user to specify the input:
% COLOR: a 1-by-3 vector to specify the RGB plot colors [Red Green Blue],
% where 0 <= Red,Green,Blue <= 1
%
% [FIG]=FCN_PLOT_POLYTOPES(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,AXIS_LIMITS)
% allows the user to specify the input:
% AXIS_LIMTS: a 1-by-4 vector to specify the start and end of the x and y
% axes, [xstart xend ystart yend]
%
% [FIG]=FCN_PLOT_POLYTOPES(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,AXIS_STYLE)
% allows the user to specify the input:
% AXIS_STYLE: controls the style of the axis, such as square or equal
%
% [FIG]=FCN_PLOT_POLYTOPES(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,FILL_INFO)
% allows the user to specify the input:
% FILL_INFO: a 1-by-5 vector to specify wether or not there is fill, the
% color of fill, and the opacity of the fill [Y/N, R, G, B, alpha]
%
% [FIG]=FCN_PLOT_POLYTOPES(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,COLOR,AXIS_LIMITS,AXIS_STYLE,FILL_INFO)
% allows the user to specify any combination of all four inputs in any
% order after LINE_WIDTH
%
% Examples:
%      
%      % BASIC example
%      cur_path = pwd;
%      main_folder = '!Voronoi Tiling Obstacles - Organized';
%      parent_dir = cur_path(1:strfind(cur_path,main_folder)-2);
%      addpath([parent_dir '\' main_folder '\Map_Generation\polytope_generation'])
%      mapx = 1;
%      mapy = 1;
%      low_pt = 1;
%      high_pt = 100;
%      [polytopes] = fcn_polytope_generation_halton_voronoi_tiling(low_pt,high_pt);
%      fig1 = fcn_plot_polytopes(polytopes,[],'-',2,[0.5 0 0]);
%      fig2 = fcn_plot_polytopes(polytopes,998,'-',2,[0 0 0.5],[0 mapx 0 mapy]);
%      fig3 = fcn_plot_polytopes(polytopes,999,'-',2,[0 0.5 0],[0 mapx 0 mapy],'square');
%      fig4 = fcn_plot_polytopes(polytopes,1000,'-',2,[0 0 0],[0 mapx 0 mapy],'square',[1 0 0 0 0.5]);
%      low_pt = 101;
%      high_pt = 200;
%      [polytopes2] = fcn_polytope_generation_halton_voronoi_tiling(low_pt,high_pt);
%      fcn_plot_polytopes(polytopes2,fig1,'r-',2)
%      fcn_plot_polytopes(polytopes2,fig2,'b--',2,[0 mapx 0 mapy])
%      fcn_plot_polytopes(polytopes2,fig3,'g-',3,[0 mapx 0 mapy],'square')
%      fcn_plot_polytopes(polytopes2,fig4,'k-',3,[0 mapx 0 mapy],'square',[1 0 0 0 0.5])
%       
% 
% This function was written on 2018_12_10 by Seth Tau
% Questions or comments? sat5340@psu.edu 
%

%% chec input arguments
if nargin < 4 || nargin > 8
    error('Incorrect number of input arguments')
end

%% open figures
if isempty(fig_num)
    fig = figure; % create new figure with next default index
else
    fig = figure(fig_num); % open specific figure
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
        filler.FaceAlpha = fill_info(5);
    end
end
if plots == 1 % basic plot
    for polys = 1:size(polytopes,2) % plot each polytope
        plot(polytopes(polys).vertices(:,1),polytopes(polys).vertices(:,2),line_spec,'linewidth',line_width)
    end
else
    for polys = 1:size(polytopes,2) % plot each polytope
        plot(polytopes(polys).vertices(:,1),polytopes(polys).vertices(:,2),line_spec,'Color',color,'linewidth',line_width)
    end
end

axis(axis_limits);
axis(axis_style);


