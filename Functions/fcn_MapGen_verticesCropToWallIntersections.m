function [croppedVertices, NwallsHit] = ...
   fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, varargin)

% fcn_MapGen_verticesCropToWallIntersections
% Crops a vertex string by the walls
%
% FORMAT:
%
%     [croppedVertices] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, (fig_num))
%
% INPUTS:
%
%     inputVertices: a list of vertex points defining, usually, a segment
%     of a polytope. These are listed as an N x 2 array, where in [x y]
%     format, where x and y are columns with at least 2 or more rows
% 
%     walls: a list of the walls of a bounding area used to crop the vertex
%     string. Only the vertices inside the walls are kept along with the
%     entry/exit points in/out of the walls. These are listed as an N x 2
%     array, where in [x y] format, where x and y are columns with at least
%     2 or more rows. Note: it is assumed that the walls are connected.
% 
%     (optional inputs)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
%
% OUTPUTS:
%
%     croppedVertices: a list of vertex points defining only the portion
%     of the vertex string within the walls, including the entry and exit
%     points.These are listed as an N x 2 array, where in [x y]
%     format, where x and y are columns
%
%     NwallsHit: the count of the different number of walls that were hit
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_Path_findSensorHitOnWall
%     inpolygon
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_verticesCropToWallIntersections
% for a full test suite.
%
% This function was written on 2021_07_15 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

%
% REVISION HISTORY:
%
% 2021_07_15 by Sean Brennan
% -- first write of function
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% -- remove dependence on test fixture
% 2025_07_11 - S. Brennan, sbrennan@psu.edu
% -- deprecated INTERNAL_fcn_geometry_findIntersectionOfSegments
%    % Using fcn_Path_findSensorHitOnWall instead (more stable)
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

        % Check the inputVertices input, 2 columns with 2 or more rows
        fcn_DebugTools_checkInputsToFunctions(...
            inputVertices, '2column_of_numbers',[2 3]);

        % Check the walls input, make sure it is '2column_of_numbers' type,
        % with 2 or more rows
        fcn_DebugTools_checkInputsToFunctions(...
            walls, '2column_of_numbers',[2 3]);

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

% STEP 1: Loop through all the vertices, looking for intersections
Nvertices = length(inputVertices(:,1));
Nwalls = length(walls(:,1))-1;
flags_in_walls = zeros(Nvertices,1);

croppedVertices = [];
walls_hit = zeros(Nwalls,1);

for ith_vertex = 1:Nvertices
    % Check if this vertex is inside a wall?
    flags_in_walls(ith_vertex) = inpolygon(inputVertices(ith_vertex,1),inputVertices(ith_vertex,2),walls(:,1),walls(:,2));
    
    if flags_in_walls(ith_vertex)
        croppedVertices = [croppedVertices; inputVertices(ith_vertex,:)]; %#ok<AGROW>
    end
    
    % Find intersection points?
    if ith_vertex~=Nvertices
        sensor_vector_start = inputVertices(ith_vertex,:);
        sensor_vector_end = inputVertices(ith_vertex+1,:);

        % Call a function that determines where and if the sensor crosses the
        % walls
        [~, wall_hit_locations, wall_that_was_hit] = ...
            fcn_Path_findSensorHitOnWall(...
            walls(1:end-1,:),...     % wall start
            walls(2:end,:),...       % wall end
            sensor_vector_start,...  % sensor_vector_start
            sensor_vector_end,...    % sensor_vector_end
            (1), ...                 % (flag_search_return_type) -- 1 means ALL hits of any results,
            (0), ...                 % (flag_search_range_type)  -- 0 means only if overlapping wall/sensor, ...
            ([]),...                 % (tolerance) -- default is eps * 1000,
            (-1));                   % (fig_num) -- -1 means to use "fast mode")

        if ~all(isnan(wall_hit_locations))
            croppedVertices = [croppedVertices; wall_hit_locations]; %#ok<AGROW>
            walls_hit(wall_that_was_hit)=1;
        end
    end
end
   
croppedVertices = unique(croppedVertices,'rows','stable');
NwallsHit = sum(walls_hit);

%§
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
    figure(fig_num);
    clf;
    hold on;
    axis equal
    grid on;
    
    % Plot the inputVertices
    plot(inputVertices(:,1),inputVertices(:,2),'ko-','LineWidth',3,'DisplayName','inputVertices');
    
    % Plot the walls
    plot(walls(:,1),walls(:,2),'b.-','LineWidth',3,'MarkerSize',20,'DisplayName','walls');
    
    if ~isempty(croppedVertices)
        % Plot the croppedVertices
        plot(croppedVertices(:,1),croppedVertices(:,2),'g.-','Linewidth',1,'MarkerSize',40,'DisplayName','croppedVertices');
    end

    legend('Interpreter','none','Location','best');
       
end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends the function

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

