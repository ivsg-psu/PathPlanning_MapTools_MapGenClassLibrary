function [ snapPoint, wallNumber] = fcn_MapGen_snapToAABB( ...
    axisAlignedBoundingBox, ...
    testPoint, ...
    varargin...
    )
% fcn_MapGen_snapToAABB
% Given an axis-aligned bounding box (AABB), and a test point, returns a
% snap point representing the contact point on the closest wall to the
% test point.
%
% FORMAT:
%
%    [ ...
%    snapPoint,...
%    wallNumber,...
%    ] = ...
%    fcn_MapGen_snapToAABB( ...
%    axisAlignedBoundingBox, ...
%    testPoint, ...
%    (snapType),...
%    (fig_num) ...
%    )
%
% INPUTS:
%
%     axisAlignedBoundingBox: the axis-aligned bounding box, in format
%     [xmin ymin xmax ymax]
%
%     testPoint: the test point, in format [x y]
%
%     (optional inputs)
%
%     snapType: 
%         0 - (default) snap to projection from middle of AABB wall
%         1 - snap to closest wall, 
%         2 - project the vector to closest wall. 
%             for option 2: testPoint must have 2 rows representing
%             [start; end] of vector
%
%     (optional inputs)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
% OUTPUTS:
%
%     snapPoint: the resulting snap point, in format [x y]
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_isWithinABBB
%     fcn_Path_findSensorHitOnWall
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
% 2021_07_14
% -- added snap type to allow projections from center
% 2021_07_17
% -- added vector projection snap type
% 2023_02_22
% -- fixed snapType bug, where 2 was specified but 3 was used. Made all 2.
% -- better comments
% -- error check on snap type 2 to force 2 rows
% -- switched over to fcn_DebugTools_checkInputsToFunctions
% 2025_04_24
% -- fixed comments and argument listings
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_07_14 by Sean Brennan
% -- removed internal call to geometry function, and replaced with external
%    % call to path library (this is supported and better tested)


% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
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
        narginchk(2,4);

        % Check the axisAlignedBoundingBox input, make sure it is
        % '4column_of_numbers' type
        fcn_DebugTools_checkInputsToFunctions(...
            axisAlignedBoundingBox, '4column_of_numbers',1);

        % Check the testPoint input, make sure it is '2column_of_numbers' type
        % with 1 or more rows
        fcn_DebugTools_checkInputsToFunctions(...
            testPoint, '2column_of_numbers',[1 2]);

        % Check the testPoint input, make sure it is '2column_of_numbers' type
        % with 2 or less rows
        fcn_DebugTools_checkInputsToFunctions(...
            testPoint, '2column_of_numbers',[2 1]);

    end
end

flag_snapType = 0;
if  3<= nargin
    temp = varargin{1};
    if ~isempty(temp)
        flag_snapType = temp;
        if 2 == flag_snapType && (1 == flag_check_inputs)
            % Check the flag_snapType input, make sure it is '1column_of_numbers' type
            % with 2 rows
            fcn_DebugTools_checkInputsToFunctions(flag_snapType, '1column_of_numbers',1);
        end
    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  (4 == nargin) && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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
walls = [...
    axisAlignedBoundingBox(1,1) axisAlignedBoundingBox(1,2); ...
    axisAlignedBoundingBox(1,3) axisAlignedBoundingBox(1,2); ...
    axisAlignedBoundingBox(1,3) axisAlignedBoundingBox(1,4); ...
    axisAlignedBoundingBox(1,1) axisAlignedBoundingBox(1,4); ...
    axisAlignedBoundingBox(1,1) axisAlignedBoundingBox(1,2)];

center = [mean([axisAlignedBoundingBox(1,1) axisAlignedBoundingBox(1,3)]),mean([axisAlignedBoundingBox(1,2) axisAlignedBoundingBox(1,4)])];
vector = testPoint - center;
angle = atan2(vector(2),vector(1));

% Is the point within the AABB?
if fcn_MapGen_isWithinABBB(axisAlignedBoundingBox,testPoint, -1)
    
    % Snap via projection, or nearest wall?
    if flag_snapType == 0    % Use projection from center of AABB to wall
        [rawDistance,rawSnapPoint,rawWallNumber] = ...
        fcn_Path_findSensorHitOnWall(...
        walls(1:end-1,:),...     % wall start
        walls(2:end,:),...       % wall end
        center,...               % sensor_vector_start
        testPoint,...            % sensor_vector_end
        (1), ...                 % (flag_search_return_type) -- 1 means ALL hits of any results,
        (1), ...                 % (flag_search_range_type)  -- 1 means ANY   projection of the sensor used with GIVEN projection of wall, ...
        ([]),...                 % (tolerance) -- default is eps * 1000,
        (-1));                   % (fig_num) -- -1 means to use "fast mode")
        
        [sortedDistances, sortedIndices] = sort(rawDistance);
        minIndex = find(sortedDistances>=0,1,'first');
        bestIndex = sortedIndices(minIndex);
        snapPoint = rawSnapPoint(bestIndex,:);
        wallNumber = rawWallNumber(bestIndex,:);

        
    elseif flag_snapType == 1     % Use nearest wall
        snapPoint = testPoint;
        if angle>=-pi/4 && angle<pi/4  % This is the x-max wall
            snapPoint(1,1) = axisAlignedBoundingBox(1,3);
            wallNumber = 2;
        elseif angle>=pi/4 && angle<pi*3/4 % This is the y-max wall
            snapPoint(1,2) = axisAlignedBoundingBox(1,4);
            wallNumber = 3;
        elseif angle>=-3*pi/4 && angle<(-pi/4) % This is the y-min wall
            snapPoint(1,2) = axisAlignedBoundingBox(1,2);
           wallNumber = 1;
         else % This is the x-min wall
            snapPoint(1,1) = axisAlignedBoundingBox(1,1);
           wallNumber = 4;
         end
    elseif flag_snapType == 2    % Use user-entered vector projection

        [rawDistance,rawSnapPoint,rawWallNumber] = ...
            fcn_Path_findSensorHitOnWall(...
            walls(1:end-1,:),...     % wall start
            walls(2:end,:),...       % wall end
            testPoint(1,:),...       % sensor_vector_start
            testPoint(2,:),...       % sensor_vector_end
            (1), ...                 % (flag_search_return_type) -- 1 means ALL hits of any results,
            (1), ...                 % (flag_search_range_type)  -- 1 means ANY   projection of the sensor used with GIVEN projection of wall, ...
            ([]),...                 % (tolerance) -- default is eps * 1000,
            (-1));                   % (fig_num) -- -1 means to use "fast mode")

        [sortedDistances, sortedIndices] = sort(rawDistance);
        minIndex = find(sortedDistances>=0,1,'first');
        bestIndex = sortedIndices(minIndex);
        snapPoint = rawSnapPoint(bestIndex,:);
        wallNumber = rawWallNumber(bestIndex,:);

    else
        error('Unrecognized projection type!');
    end
else % Point is outside the box - no need to snap
    [rawDistance,~,rawWallNumber] = ...
        fcn_Path_findSensorHitOnWall(...
        walls(1:end-1,:),...     % wall start
        walls(2:end,:),...       % wall end
        center,...               % sensor_vector_start
        testPoint,...            % sensor_vector_end
        (1), ...                 % (flag_search_return_type) -- 1 means ALL hits of any results,
        (1), ...                 % (flag_search_range_type)  -- 1 means ANY   projection of the sensor used with GIVEN projection of wall, ...
        ([]),...                 % (tolerance) -- default is eps * 1000,
        (-1));                   % (fig_num) -- -1 means to use "fast mode")
        
        [sortedDistances, sortedIndices] = sort(rawDistance);
        minIndex = find(sortedDistances>=0,1,'first');
        bestIndex = sortedIndices(minIndex);
        wallNumber = rawWallNumber(bestIndex,:);

    snapPoint = testPoint;
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
    box_outline = [axisAlignedBoundingBox(1,1) axisAlignedBoundingBox(1,2); axisAlignedBoundingBox(1,3) axisAlignedBoundingBox(1,2); axisAlignedBoundingBox(1,3) axisAlignedBoundingBox(1,4); axisAlignedBoundingBox(1,1) axisAlignedBoundingBox(1,4); axisAlignedBoundingBox(1,1) axisAlignedBoundingBox(1,2)];
    plot(box_outline(:,1),box_outline(:,2),'-');
    
    % Plot the test point(s)
    plot(testPoint(:,1),testPoint(:,2),'o-');
    
    % Plot the snap point
    plot(snapPoint(:,1),snapPoint(:,2),'x');
    
    % Plot the snap point to test point line
    plot([testPoint(1,1) snapPoint(1,1)],[testPoint(1,2) snapPoint(:,2)],'-');
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

