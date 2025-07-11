function [verticesIncludingSelfIntersections, flag_wasIntersection] = ...
    fcn_MapGen_polytopeFindSelfIntersections(vertices,...
    varargin)
% fcn_MapGen_polytopeFindSelfIntersections finds the points where the
% vertices of a polytope create self-intersections, e.g. where the "walls"
% of the polytope cross each other. The algorithm returns all the points
% where this occurs within the vertices, in order. In other words, if
% vertex 1 to vertex 2 crosses an edge, then the output would be:
% [vertex1; crossing; vertex2; (etc)]
%
% FORMAT:
% 
%    [verticesIncludingSelfIntersections, flag_wasIntersection] = ...
%    fcn_MapGen_polytopeFindSelfIntersections(vertices, (fig_num))
%
% INPUTS:
%
%     vertices: an Nx2 matrix of [x y] vertices
%
%    (OPTIONAL INPUTS)
%
%    fig_num: a figure number to plot results. If set to -1, skips any
%    input checking or debugging, no figures will be generated, and sets
%    up code to maximize speed. As well, if given, this forces the
%    variable types to be displayed as output and as well makes the input
%    check process verbose.
%
% OUTPUTS:
%
%     verticesIncludingSelfIntersections: an Mx2 matrix of [x y] points,
%     which include both the input vertices and added crossings
%     (self-intersections) where walls intersect with other walls
%
%     flag_wasIntersection: a flag that is 0 if there were no self
%     intersections detected, 1 if there were
%   
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_Path_findSensorHitOnWall
% 
% EXAMPLES:
%      
% For additional examples, see:
% script_test_fcn_MapGen_polytopeFindSelfIntersections
%
% This function was written on 2021_08_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu 
%

% Revision History:
% 2021_08_02 - S. Brennan
% -- first write of code
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_07_10 by Sean Brennan
% -- updated variable names for clarity
% -- changed fcn_MapGen_findIntersectionOfSegments to use
% fcn_Path_findSensorHitOnWall instead, as the Path function is much more
% tested/debugged and regularly updated
% -- fixed bug with flag_wasIntersection

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==2 && isequal(varargin{end},-1))
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
        narginchk(1,2);

        % Check the vertices input
        fcn_DebugTools_checkInputsToFunctions(...
            vertices, '2column_of_numbers');

    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  2 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
    temp = varargin{end}; % Last argument is always figure number
    if ~isempty(temp) % Make sure the user is not giving empty input
        fig_num = temp;
        flag_do_plot = 1; % Set flag to do plotting
    end
else
    if flag_do_debug % If in debug mode, do plotting but to an arbitrary figure number
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

% Finds all vertices where the polytope intersects itself. Assumes the
% vertices already wrap around

% Fill in blank verticesIncludingSelfIntersections as starter
verticesIncludingSelfIntersections = [];

start_vertices = vertices(1:end-1,:);
end_vertices = vertices(2:end,:);
all_indices = (1:(length(vertices(:,1))-1))'; 

flag_wasIntersection = 0;
% Find all intersection points
for ith_point = 1:length(start_vertices(:,1))
    verticesIncludingSelfIntersections = [verticesIncludingSelfIntersections; start_vertices(ith_point,:)]; %#ok<AGROW>
    
    % Define start and end of the sensor
    sensor_vector_start = start_vertices(ith_point,:);
    sensor_vector_end = end_vertices(ith_point,:);
    
    % Define the start and end of the walls
    wall_indices = all_indices(all_indices~=ith_point);
    wall_start   = start_vertices(wall_indices,:);    
    wall_end     = end_vertices(wall_indices,:);
    midpoints = (wall_start+wall_end)/2; 
    
    if flag_do_debug
        figure(fig_for_debug);
        clf;
        hold on;
        grid on;
        axis equal;
        grid minor;
        
        % Plot the vertices
        plot(vertices(:,1),vertices(:,2),'r-');
        
        % Plot the sensor
        plot(...
            [sensor_vector_start(1,1) sensor_vector_end(1,1)],...
            [sensor_vector_start(1,2) sensor_vector_end(1,2)],...
            'g-');
        
        % Plot the walls, and number them
        Nwalls = length(wall_start(:,1));
        for ith_wall = 1:Nwalls
            handle_plot = plot(...
                [wall_start(ith_wall,1) wall_end(ith_wall,1)],...
                [wall_start(ith_wall,2) wall_end(ith_wall,2)],...
                '.-','Linewidth',5);
            line_color = get(handle_plot,'Color');
            handle_text = text(midpoints(ith_wall,1),midpoints(ith_wall,2),sprintf('%.0d',ith_wall));
            set(handle_text,'Color',line_color);
            set(handle_text,'Fontsize',20);
        end
    end


    % Call a function that determines where and if the sensor crosses the
    % walls
    [distance, location] = ...
        fcn_Path_findSensorHitOnWall(...
        wall_start,...           % wall start
        wall_end,...             % wall end
        sensor_vector_start,...  % sensor_vector_start
        sensor_vector_end,...    % sensor_vector_end
        (1), ...                 % (flag_search_return_type) -- 1 means ALL hit of any results,
        (0), ...                 % (flag_search_range_type)  -- 0 means only if overlapping wall/sensor, ...
        ([]),...                 % (tolerance) -- default is eps * 1000,
        (-1));                   % (fig_num) -- -1 means to use "fast mode")
    
    % Sort by distances
    sorting_data = [distance location];
    sorted_data = sortrows(sorting_data,1); % sort by 1st column
    distance = sorted_data(:,1);
    location = sorted_data(:,2:3);
    
    if ~isnan(distance)
        verticesIncludingSelfIntersections = [verticesIncludingSelfIntersections; location(distance~=0,:)]; %#ok<AGROW>
    end
    
end

% Get rid of duplicates. This happens when two points are both on edges,
% and at every end of connecting segments
indices_not_repeated = [~all(abs(diff(verticesIncludingSelfIntersections))<eps*10,2); 1];
verticesIncludingSelfIntersections = verticesIncludingSelfIntersections(indices_not_repeated>0,:);

if length(verticesIncludingSelfIntersections(:,1))==length(vertices(:,1))
    flag_wasIntersection = 0;
else
    flag_wasIntersection = 1;
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
    
    % Plot the vertices
    plot(vertices(:,1),vertices(:,2),'r-');
    
    % Plot the verticesIncludingSelfIntersections
    plot(verticesIncludingSelfIntersections(:,1),verticesIncludingSelfIntersections(:,2),'ko');

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function


