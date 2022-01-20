function [ ...
isInside ...
] = ...
fcn_MapGen_isWithinABBB( ...
AABB, ...
test_points, ...
varargin...
)
% fcn_MapGen_isWithinABBB
% Checks if the points are within the given axis-aligned bounding box,
% AABB, returning a vector of 1' or 0's the same length as the nubmer of
% rows of points. Each point must be strictly within the AABB - e.g. this
% function returns "false" if a point is on the "wall" of the AABB.
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
%     isInside: a column of 1's or 0's, one for each test point, with 1
%     meaning that the test point is within the AABB
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
    fig_for_debug = 225;
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


if 1 == flag_check_inputs

    % Are there the right number of inputs?
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end

    % Check the AABB input, make sure it is '4column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        AABB, '4column_of_numbers',1);

    % Check the test_points input, make sure it is '2column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        test_points, '2column_of_numbers');

end

% Does user want to show the plots?
if  3== nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
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





% % See: https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection
% % for details on axis-aligned bounding boxes (AABB)

isInside = (test_points(:,1)>AABB(1,1))  & ...
    (test_points(:,2)>AABB(1,2))  & ...
    (test_points(:,1)<AABB(1,3))  & ...
    (test_points(:,2)<AABB(1,4));


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



if flag_do_plot
    figure(fig_num);
    clf;
    hold on;

    % Convert axis-aligned bounding box to wall format
    walls = [AABB(1,1) AABB(1,2); AABB(1,3) AABB(1,2); AABB(1,3) AABB(1,4); AABB(1,1) AABB(1,4); AABB(1,1) AABB(1,2)];

    % Plot the walls
    plot(walls(:,1),walls(:,2),'k-');

    % Plot the test_points

    % plot(...
    %     [test_points(:,1); test_points(1,1)],...
    %     [test_points(:,2); test_points(1,2)],...
    %     '.-');
    plot(test_points(:,1), test_points(:,2),'k.');

    % Plot the interior points with green
    plot(test_points(isInside,1),test_points(isInside,2),'go');

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




