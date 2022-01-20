function [projected_points] = ...
    fcn_MapGen_polytopeProjectVerticesOntoWalls(...,
    interior_point,...
    vertices,...
    wall_start,...
    wall_end,...
    varargin)

% From the interior point, project all_points back onto the wall. The
% results are sorted as well by polar angles around the interior point.
%
% FORMAT:
%
% [projected_points] = ...
%    fcn_MapGen_polytopeProjectVerticesOntoWalls(vertices,...
%     (fig_num))
%
% INPUTS:
%
%     interior_point: an 1x2 matrix of the [x y] point that is interior to
%     the walls, e.g. the "test point".
%
%     vertices: an Nx2 matrix of [x y] vertices
%
%     wall_start: an N x 2 vector containing the X,Y points of the
%     starting points of each "wall".
%
%     wall_end: an N x 2 vector containing the X,Y points of the
%     ending points of each "wall".
%
%    (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%     projected_points: an Mx2 matrix of [x y] points, each representing
%     the first hit of projected points onto a wall. If no wall is hit,
%     then the vertex is kept as is. Points are sorted by polar angle
%     around the interior_point.
%
% DEPENDENCIES:
%
%     fcn_MapGen_checkInputsToFunctions
%     fcn_MapGen_findIntersectionOfSegments
%
% EXAMPLES:
%
% For additional examples, see:
% script_test_fcn_MapGen_polytopeProjectVerticesOntoWalls
%
% This function was written on 2021_08_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%

% Revision History:
% 2021_08_03 - S. Brennan
% -- first write of code

% TO DO
% -- none

%% Debugging and Input checks
% set an environment variable on your machine with the getenv function...
% in the Matlab console.  Char array of '1' will be true and '0' will be false.
flag_check_inputs = getenv('ENV_FLAG_CHECK_INPUTS');  % '1' will check input args
flag_do_plot = getenv('ENV_FLAG_DO_PLOT'); % '1' will make plots
flag_do_debug = getenv('ENV_FLAG_DO_DEBUG'); % '1' will enable debugging

% if the char array has length 0, assume the env var isn't set and default to...
% dipslaying more information rather than potentially hiding an issue
if length(flag_check_inputs) = 0
    flag_check_inputs = '1';
end
if length(flag_do_plot) = 0
    flag_do_plot = '1';
end
if length(flag_do_debug) = 0
    flag_do_debug = '1';
end

% convert flag from char string to logical
flag_check_inputs = flag_check_inputs == '1';
flag_do_plot = flag_do_plot == '1';
flag_do_debug = flag_do_debug == '1';

if flag_do_debug
    fig_for_debug = 4564;
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
    if nargin < 4 || nargin > 5
        error('Incorrect number of input arguments')
    end

    % Check the interior_point input 2 columns and 1 row
    fcn_MapGen_checkInputsToFunctions(...
        interior_point, '2column_of_numbers',1);

    % Check the vertices input
    fcn_MapGen_checkInputsToFunctions(...
        vertices, '2column_of_numbers');

    % Check the wall_start input
    fcn_MapGen_checkInputsToFunctions(...
        wall_start, '2column_of_numbers');

    % Check the wall_end input
    fcn_MapGen_checkInputsToFunctions(...
        wall_end, '2column_of_numbers',length(wall_start(:,1)));
end


% Does user want to show the plots?
if  5 == nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
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


% Convert the input xy data into polar coordinates, so that we can sort by
% theta. If we don't do this, then sometimes (rarely) the walls of the
% polytope will jump around producing very weird results.
offset_points = vertices - interior_point;
[theta,~]=cart2pol(offset_points(:,1),offset_points(:,2));
[~,index_sorted] = sort(theta,'ascend');
all_points = vertices(index_sorted,:);

% Now see where the points hit via projection
projected_points = 0*all_points;
for ith_point = 1:length(all_points(:,1))

    % Define start and end of the sensor
    sensor_vector_start = interior_point;
    sensor_vector_end = all_points(ith_point,:);

    % Find hits on the walls
    [distance,location,~] = ...
        fcn_MapGen_findIntersectionOfSegments(...
        wall_start,...
        wall_end,...
        sensor_vector_start,...
        sensor_vector_end);

    % Did we hit anything?
    if ~isnan(distance)
        projected_points(ith_point,:) = location;
    else
        projected_points(ith_point,:) = sensor_vector_end;
    end

end


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
    figure(fig_num);
    grid on
    grid minor
    hold on
    axis equal

    % Plot the input vertices
    plot(vertices(:,1),vertices(:,2),'r-');

    % Plot the projected_points
    plot(projected_points(:,1),projected_points(:,2),'ko');

    % Number the points
    temp = axis;
    nudge = (temp(2)-temp(1))/100;
    for ith_point = 1:length(projected_points(:,1))
        text(projected_points(ith_point,1)+nudge,projected_points(ith_point,2),sprintf('%.0d',ith_point));
    end

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function


