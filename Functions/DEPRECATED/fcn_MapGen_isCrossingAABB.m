function [isInside] = fcn_MapGen_isCrossingAABB(AABB, testPoints, varargin)

warning('on','backtrace');
warning('fcn_MapGen_isCrossingAABB is being deprecated. Use fcn_Path_findSensorHitOnWall instead.');

% fcn_MapGen_isCrossingAABB
% Checks if the line segments between every pair of test points are intersecting
% the given axis-aligned bounding box,
% AABB, returning a matrix of 1' or 0's the same number of rows and columns
% as the number of points where the value at element i,j is 1 if the line
% segment starting at point i and ending at point j intersects
% the AABB. Each line segment must be strictly within the AABB
% - e.g. this function returns "false" if a point is on the "wall" of the AABB.
%
%
%
% FORMAT:
%
%    [ ...
%    isInside ...
%    ] = ...
%    fcn_MapGen_isWithinABBB( ...
%    AABB, ...
%    testPoints, ...
%    (fig_num) ...
%    )
%
% INPUTS:
%
%     AABB: the Axis-Aligned Bounding Box, defined in form of [xmin ymin
%     xmax ymax]
%
%     testPoints: the test points to check, in form of [x y] where x and
%     y are scalar or column vectors
%
%     (optional inputs)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
%
% OUTPUTS:
%
%     isInside: a matrix of 1's or 0's, one for each pair of test points, with 1
%     meaning that the line segment from test point i to test point j is within
%     intersecting the AABB
%
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_isCrossingAABB
% for a full test suite.
%
% This function was written on 2024_03_20 by Steve Harnett
% Questions or comments? contact sjharnett@psu.edu

%
% REVISION HISTORY:
%
% 2024_03_20
% -- first written by S. Harnett
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% -- fixed call to fcn_MapGen_fillPolytopeFieldsFromVertices


% TO DO
% -- this function it's currently implemented is doing a much coarser check than is necessary.
%    It checks if a line segment COULD BE hitting an AABB.  In other words, this function indicating
%    no hit means the line semgment cannot hit the AABB but this function indicating hit means the
%    line segment may or may not hit the AABB and further collision checking is necessary.  A better
%    way to implement this is to implement line-side collision checking like this:
%    https://www.jeffreythompson.org/collision-detection/line-rect.php which requires this:
%    https://www.jeffreythompson.org/collision-detection/line-line.php

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
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
        if nargin < 2 || nargin > 3
            error('Incorrect number of input arguments')
        end

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

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  3 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
    temp = varargin{end}; % Last argument is always figure number
    if ~isempty(temp) % Make sure the user is not giving empty input
        fig_num = temp;
        flag_do_plot = 1; % Set flag to do plotting
    end
else
    if flag_do_debug % If in debug mode, do plotting but to an arbitrary figure number
        fig = figure;
        fig_for_debug = fig.Number; %#ok<NASGU>
        flag_do_plot = 1;
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

num_pts = size(testPoints,1); % number of rows in test points
xs = testPoints(:,1); % extract all x coords
ys = testPoints(:,2); % extract all y coords
% make a matrix of all x coord down the col, repeated
% also make a matrix of all x coord along the row, repeated
% then do the same for y
xcols = repmat(xs,1,num_pts);
xrows= repmat(xs',num_pts,1);
ycols = repmat(ys,1,num_pts);
yrows= repmat(ys',num_pts,1);
xmin_mat = min(xrows,xcols);
xmax_mat = max(xrows,xcols);
ymin_mat = min(yrows,ycols);
ymax_mat = max(yrows,ycols);

xmin = AABB(1);
ymin = AABB(2);
xmax = AABB(3);
ymax = AABB(4);

%% there are four ways the edge could not collide with the AABB
% the whole line segment is left of the AABB
edge_is_left = xmax_mat <= xmin;
% the whole line segment is right of the AABB
edge_is_right = xmin_mat >= xmax;
% the whole line segment is above the AABB
edge_is_above = ymin_mat >= ymax;
% the whole line segment is below the AABB
edge_is_below = ymax_mat <= ymin;
%% if none of these cases are true, the edge is either entirely within the AABB, partially within the AABB, or straddling the AABB
isInside = ~(edge_is_left | edge_is_right | edge_is_above | edge_is_below);

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



if flag_do_plot
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
