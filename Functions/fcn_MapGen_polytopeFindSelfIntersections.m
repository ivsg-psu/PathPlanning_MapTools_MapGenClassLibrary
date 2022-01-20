function [all_points, flag_was_intersection] = ...
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
% [all_points, flag_was_intersection] = ...
%    fcn_MapGen_polytopeFindSelfIntersections(vertices,...
%     (fig_num))
%
% INPUTS:
%
%     vertices: an Nx2 matrix of [x y] vertices
%
%    (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%     all_points: an Mx2 matrix of [x y] points, which include the vertices
%     and point crossings
%
% DEPENDENCIES:
%
%     fcn_MapGen_checkInputsToFunctions
%     fcn_MapGen_polytopeFindVertexAngles
%     fcn_MapGen_fillPolytopeFieldsFromVertices
%
% % EXAMPLES:
%
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
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end

    % Check the vertices input
    fcn_MapGen_checkInputsToFunctions(...
        vertices, '2column_of_numbers');

end


% Does user want to show the plots?
if  2 == nargin
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

% Finds all vertices where the polytope intersects itself. Assumes the
% vertices already wrap around

% Fill in blank all_points as starter
all_points = [];

start_vertices = vertices(1:end-1,:);
end_vertices = vertices(2:end,:);
all_indices = (1:(length(vertices(:,1))-1))';

flag_was_intersection = 0;
% Find all intersection points
for ith_point = 1:length(start_vertices(:,1))
    all_points = [all_points; start_vertices(ith_point,:)]; %#ok<AGROW>

    % Define start and end of the sensor
    sensor_vector_start = start_vertices(ith_point,:);
    sensor_vector_end = end_vertices(ith_point,:);

    % Define the start and end of the walls
    wall_indices = all_indices(all_indices~=ith_point);
    wall_start   = start_vertices(wall_indices,:);
    wall_end     = end_vertices(wall_indices,:);
    midpoints = (wall_start+wall_end)/2; %#ok<NASGU>

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
    [distance,location,~] = ...
        fcn_MapGen_findIntersectionOfSegments(...
        wall_start,...
        wall_end,...
        sensor_vector_start,...
        sensor_vector_end,2);

    % Sort by distances
    sorting_data = [distance location];
    sorted_data = sortrows(sorting_data,1); % sort by 1st column
    distance = sorted_data(:,1);
    location = sorted_data(:,2:3);

    if ~isnan(distance)
        all_points = [all_points; location(distance~=0,:)]; %#ok<AGROW>
        flag_was_intersection = 1;
    end

end

% Get rid of duplicates (occurs when two points are both on edges)
indices_not_repeated = [~all(abs(diff(all_points))<eps*10,2); 1];
all_points = all_points(indices_not_repeated>0,:);


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

    % Plot the all_points
    plot(all_points(:,1),all_points(:,2),'ko');

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function


