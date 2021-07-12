function [ ...
    snap_point ...
    ] = ...
    fcn_MapGen_snapToAABB( ...
    axis_aligned_bounding_box, ...
    test_point, ...
    varargin...
    )
% fcn_MapGen_snapToAABB
% Given an axis-aligned bounding box (AABB), and a test point, returns a
% snap point representing the contact point on the closest wall to the
% test point.
%
%
%
% FORMAT:
%
%    [ ...
%    snap_point ...
%    ] = ...
%    fcn_MapGen_snapToAABB( ...
%    axis_aligned_bounding_box, ...
%    test_point, ...
%    (fig_num) ...
%    )
%
% INPUTS:
%
%     axis_aligned_bounding_box: the axis-aligned bounding box, in format
%     [xmin ymin xmax ymax]
%
%     test_point: the test point, in format [x y]
%
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a
%     figure to be drawn showing results.
%
%
% OUTPUTS:
%
%     snap_point: the resulting snap point, in format [x y]
%
%
% DEPENDENCIES:
%
%     (none)
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_snapToAABB
% for a full test suite.
%
% This function was written on 2021_07_02 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

%
% REVISION HISTORY:
%
% 2021_07_02 by Sean Brennan
% -- first write of function

%
% TO DO:
%
% -- fill in to-do items here.

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    fig_for_debug = 128;
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
    
    % Check the axis_aligned_bounding_box input, make sure it is '4column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        axis_aligned_bounding_box, '4column_of_numbers',1);
    
    % Check the test_point input, make sure it is '2column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        test_point, '2column_of_numbers',1);
    
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

if fcn_MapGen_isWithinABBB(axis_aligned_bounding_box,test_point)
    center = [mean([axis_aligned_bounding_box(1,1) axis_aligned_bounding_box(1,3)]),mean([axis_aligned_bounding_box(1,2) axis_aligned_bounding_box(1,4)])];
    vector = test_point - center;
    angle = atan2(vector(2),vector(1));
    
    snap_point = test_point;
    if angle>=-pi/4 && angle<pi/4  % This is the x-max wall
        snap_point(1,1) = axis_aligned_bounding_box(1,3);
    elseif angle>=pi/4 && angle<pi*3/4 % This is the y-max wall
        snap_point(1,2) = axis_aligned_bounding_box(1,4);
    elseif angle>=-3*pi/4 && angle<(-pi/4) % This is the y-min wall
        snap_point(1,2) = axis_aligned_bounding_box(1,2);
    else % This is the x-min wall
        snap_point(1,1) = axis_aligned_bounding_box(1,1);
    end
else % Point is outside the box - no need to snap
    snap_point = test_point;
end

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
    axis equal
    grid on;
    
    % Plot the bounding box
    box_outline = [axis_aligned_bounding_box(1,1) axis_aligned_bounding_box(1,2); axis_aligned_bounding_box(1,3) axis_aligned_bounding_box(1,2); axis_aligned_bounding_box(1,3) axis_aligned_bounding_box(1,4); axis_aligned_bounding_box(1,1) axis_aligned_bounding_box(1,4); axis_aligned_bounding_box(1,1) axis_aligned_bounding_box(1,2)];
    plot(box_outline(:,1),box_outline(:,2),'-');
    
    % Plot the test point
    plot(test_point(:,1),test_point(:,2),'o');
    
    % Plot the snap point
    plot(snap_point(:,1),snap_point(:,2),'x');
    
    % Plot the snap point
    plot([test_point(:,1) snap_point(1,1)],[test_point(:,2) snap_point(:,2)],'-');
    %
    
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


function isInside = fcn_MapGen_isWithinABBB(box,test_point)
% Checks if the point is within the AABB?
% % See: https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection
% % for details on axis-aligned bounding boxes (AABB)

isInside = (test_point(1,1)>box(1,1))  && ...
        (test_point(1,2)>box(1,2))  && ...
        (test_point(1,1)<box(1,3))  && ...
        (test_point(1,2)<box(1,4));
end


