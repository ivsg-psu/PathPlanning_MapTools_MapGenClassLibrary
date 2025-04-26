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
%    fcn_MapGen_polytopeProjectVerticesOntoWalls(...
%     interior_point, vertices, wall_start, wall_end, ...
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
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
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
%     fcn_DebugTools_checkInputsToFunctions
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
if (nargin==5 && isequal(varargin{end},-1))
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
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(4,5)

        % Check the interior_point input 2 columns and 1 row
        fcn_DebugTools_checkInputsToFunctions(...
            interior_point, '2column_of_numbers',1);

        % Check the vertices input
        fcn_DebugTools_checkInputsToFunctions(...
            vertices, '2column_of_numbers');

        % Check the wall_start input
        fcn_DebugTools_checkInputsToFunctions(...
            wall_start, '2column_of_numbers');

        % Check the wall_end input
        fcn_DebugTools_checkInputsToFunctions(...
            wall_end, '2column_of_numbers',length(wall_start(:,1)));
    end
end
    

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  5 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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


