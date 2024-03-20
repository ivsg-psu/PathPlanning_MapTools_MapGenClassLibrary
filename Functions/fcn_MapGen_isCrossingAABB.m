function [isInside] = fcn_MapGen_isCrossingAABB(AABB, test_points)
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
%    test_points, ...
%    (fig_num) ...
%    )
%
% INPUTS:
%
%     AABB: the Axis-Aligned Bounding Box, defined in form of [xmin ymin
%     xmax ymax]
%
%     test_points: the test points to check, in form of [x y] where x and
%     y are scalar or column vectors
%
%     (optional inputs)
%
%     fig_num: any number that acts somewhat like a figure number output.
%     If given, this forces the variable types to be displayed as output
%     and as well makes the input check process verbose.
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
%     fcn_MapGen_checkInputsToFunctions
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_isWithinABBB
% for a full test suite.
%
% This function was written on 2021_07_11 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

%
% REVISION HISTORY:
%
% 2021_07_11 by Sean Brennan
% -- first write of function
% 2021_07_17 by Sean Brennan
% -- clarified strictness of the AABB

%
% TO DO:
%
% -- fill in to-do items here.

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

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

num_pts = size(test_points,1); % number of rows in test points
xs = test_points(:,1); % extract all x coords
ys = test_points(:,2); % extract all y coords
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
% the whole line segment i=s below the AABB
edge_is_below = ymax_mat <= ymin;
%% if none of these cases are true, the edge is either entirely within the AABB, partially within the AABB, or straddling the AABB
isInside = ~(edge_is_left | edge_is_right | edge_is_above | edge_is_below);
end
