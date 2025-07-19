function [vx,vy,h] = fcn_MapGen_generateVoronoiDiagramBetweenPolytopes(polytopes,is_nonconvex, varargin)
% fcn_MapGen_generateVoronoiDiagramBetweenPolytopes
% Wraps the matlab voronoi() function to find the voronoi diagram using
% the vertices of polytopes as the seed points, in the case of convex obstacles
% or using densely packed colinear vertices along polytope sides as seed points
% in the case of nonconvex osbtacles.  Voronoi diagram boundaries may collide
% with polytope sides but voronoi diagram boundaries bewteen obstacles are also
% calcualted.  This is similar to the process used by Masehian et al. 2004, see:
% Masehian, Ellips, and M. R. Amin‐Naseri. "A voronoi diagram‐visibility graph‐potential field compound algorithm for robot path planning." Journal of Robotic Systems 21.6 (2004): 275-300.
% To generate the medial axis between polytopes without these errant edges, you may
% wish to use the fcn_MedialAxis_* functions in PathPlanning_GridFreePathPlanners_BoundedAStar
%
%
% FORMAT:
% [vx,vy,h] = fcn_MapGen_generateVoronoiDiagramBetweenPolytopes(polytopes,is_nonconvex)
%
% INPUTS:
%     polytopes - the initial polytope field
%
%     is_nonconvex - boolean flag indicating if there are or are not non-convex polytopes
%
%     (optional inputs)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
% OUTPUTS:
%
%     Outputs are forwarded directly from Matlab's voronoi function.  See here for description:
%       https://www.mathworks.com/help/matlab/ref/voronoi.html
%
% DEPENDENCIES:
%     Matlab's voronoi function
%     fcn_MapGen_polytopesIncreaseVertexCount
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_generateVoronoiDiagramBetweenPolytopes.m
% for a full test suite.
%
% Questions or comments? contact sjh6473@psu.edu

% REVISION HISTORY:
% 2024_03_15
% -- first written by Steve Harnett
% 2025_04_16
% -- commented by Steve Harnett
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% -- fixed call to fcn_MapGen_polytopesFillFieldsFromVertices
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

if ~is_nonconvex
    [vx,vy] = voronoi([polytopes.xv],[polytopes.yv]);
    h = voronoi([polytopes.xv],[polytopes.yv]);
else
    % min of diff between all points /2 so at least every side is cut in half
    distances = diff([[polytopes.xv]',[polytopes.yv]']);
    min_distance_between_verts = min(sqrt(sum(distances.*distances,2)));
    % poly_map_stats = fcn_MapGen_statsPolytopes(polytopes);
    % % want to ensure that a side with length of 2 std dev below mean is still interpolated at least in half
    % resolution = (poly_map_stats.average_side_length - 2*poly_map_stats.std_side_length)/2;
    resolution = min_distance_between_verts/2;
    interpolated_polytopes = fcn_MapGen_polytopesIncreaseVertexCount(polytopes, resolution);
    [vx,vy] = voronoi([interpolated_polytopes.xv], [interpolated_polytopes.yv]);
    h = voronoi([interpolated_polytopes.xv], [interpolated_polytopes.yv]);
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
    figure(fig_num)
    clf;

    % LineWidth = 2;
    % fcn_MapGen_plotPolytopes(polytopes,fig_num,'r-',LineWidth);
    % fcn_MapGen_plotPolytopes(exp_polytopes,fig_num,'b-',LineWidth,'square');
    % legend('Original','Expanded')
    % box on
    % xlabel('X Position')
    % ylabel('Y Position')

end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends the function
