function [projectedPoints] = ...
    fcn_MapGen_polytopeProjectVerticesOntoWalls(...,
    interiorPoint,...
    vertices,...
    wall_start,...
    wall_end,...
    varargin)

% From the interior point, project all_points back onto the wall. The
% results are sorted as well by polar angles around the interior point.
%
% FORMAT:
% 
% [projectedPoints] = ...
%    fcn_MapGen_polytopeProjectVerticesOntoWalls(...
%     interiorPoint, vertices, wall_start, wall_end, ...
%     (fig_num))
%
% INPUTS:
%
%     interiorPoint: an 1x2 matrix of the [x y] point that is interior to
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
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%     projectedPoints: an Nx2 matrix of [x y] points, each representing
%     the first hit of projected points onto a wall. If no wall is hit,
%     then the vertex is kept as is. Points are sorted by polar angle
%     around the interiorPoint.
%   
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_Path_findSensorHitOnWall
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
% 2025_07_10 by Sean Brennan
% -- changed fcn_MapGen_findIntersectionOfSegments to use
% fcn_Path_findSensorHitOnWall instead, as the Path function is much more
% tested/debugged and regularly updated
% -- renamed variables for clarity
% -- improved plotting
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
MAX_NARGIN = 5; % The largest Number of argument inputs to the function
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
        narginchk(4,MAX_NARGIN);

        % Check the interiorPoint input 2 columns and 1 row
        fcn_DebugTools_checkInputsToFunctions(interiorPoint, '2column_of_numbers',1);

        % Check the vertices input, 2 columns
        fcn_DebugTools_checkInputsToFunctions(vertices, '2column_of_numbers');

        % Check the wall_start input, 2 columns
        fcn_DebugTools_checkInputsToFunctions(wall_start, '2column_of_numbers');

        % Check the wall_end input, 2 columns and same length as wall_start
        fcn_DebugTools_checkInputsToFunctions(wall_end, '2column_of_numbers',length(wall_start(:,1)));

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
offset_points = vertices - interiorPoint;
[theta,~]=cart2pol(offset_points(:,1),offset_points(:,2));
[~,index_sorted] = sort(theta,'ascend');
all_points = vertices(index_sorted,:);

% Now see where the points hit via projection
projectedPoints = 0*all_points;
for ith_point = 1:length(all_points(:,1))
    
    % Define start and end of the sensor
    sensor_vector_start = interiorPoint;
    sensor_vector_end = all_points(ith_point,:);
   
    [distance, location] = ...
        fcn_Path_findSensorHitOnWall(...
        wall_start,...           % wall start
        wall_end,...             % wall end
        sensor_vector_start,...  % sensor_vector_start
        sensor_vector_end,...    % sensor_vector_end
        (0), ...                 % (flag_search_return_type) -- 0 means first hit of any results,
        (0), ...                 % (flag_search_range_type)  -- 0 means only if overlapping wall/sensor, ...
        ([]),...                 % (tolerance) -- default is eps * 1000,
        (-1));                   % (fig_num) -- -1 means to use "fast mode")
    
    % Did we hit anything?
    if ~isnan(distance)
        projectedPoints(ith_point,:) = location;
    else
        projectedPoints(ith_point,:) = sensor_vector_end;
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

if flag_do_plots
    figure(fig_num);
    grid on
    grid minor
    hold on
    axis equal

    % Plot the input interiorPoint
    plot(interiorPoint(:,1),interiorPoint(:,2),'g.','MarkerSize', 30, 'DisplayName','Input: interiorPoint');

    % Plot the input vertices
    plot(vertices(:,1),vertices(:,2),'r.-', 'DisplayName','Input: vertices');
    
    % Number the vertices
    temp = axis;
    nudge = (temp(2)-temp(1))/100;
    for ith_point = 1:length(vertices(:,1))
        NrepeatsSoFar = nnz(ismember(vertices(1:ith_point,:), vertices(ith_point,:), 'rows'))-1;
        NrepeatsAhead = nnz(ismember(vertices(ith_point:end,:), vertices(ith_point,:), 'rows'))-1;
        if NrepeatsAhead>0
            text(vertices(ith_point,1)+nudge+NrepeatsSoFar*nudge*3,vertices(ith_point,2)+3*nudge,sprintf('%.0d,',ith_point),'Color',[1 0.5 0.5]);
        else
            text(vertices(ith_point,1)+nudge+NrepeatsSoFar*nudge*3,vertices(ith_point,2)+3*nudge,sprintf('%.0d',ith_point),'Color',[1 0.5 0.5]);
        end
    end


    % Plot the projectedPoints
    plot(projectedPoints(:,1),projectedPoints(:,2),'ko', 'DisplayName','Output: projectedPoints');
    
    % Number the projectedPoints
    temp = axis;
    nudge = (temp(2)-temp(1))/100;
    for ith_point = 1:length(projectedPoints(:,1))
        NrepeatsSoFar = nnz(ismember(projectedPoints(1:ith_point,:), projectedPoints(ith_point,:), 'rows'))-1;
        NrepeatsAhead = nnz(ismember(projectedPoints(ith_point:end,:), projectedPoints(ith_point,:), 'rows'))-1;
        if NrepeatsAhead>0
            text(projectedPoints(ith_point,1)+nudge+NrepeatsSoFar*nudge*3,projectedPoints(ith_point,2),sprintf('%.0d,',ith_point));
        else
            text(projectedPoints(ith_point,1)+nudge+NrepeatsSoFar*nudge*3,projectedPoints(ith_point,2),sprintf('%.0d',ith_point));
        end
    end

    legend('Interpreter','none');
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function