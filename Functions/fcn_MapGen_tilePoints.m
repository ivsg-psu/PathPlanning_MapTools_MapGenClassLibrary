function [tiledPoints] = ...
   fcn_MapGen_tilePoints(...
   inputPoints, tileDepth, AABB, varargin)

% fcn_MapGen_tilePoints - creates a tiling of points. 
%
% Given an input set of Nx2 vector of points that specify a points in an
% area "X", this returns a tiling of the points by repeating them but with
% coordinate displacements along the Axis-aligned Bounding Box (AABB), a
% certain tile "depth". For example, if a region "X" is specified and a
% tiling depth of 1 is used, this returns tiling points that make a 1-unit
% boundary around X, as:
%
%                 Y Y Y
%                 Y X Y
%                 Y Y Y
%
% For a tiling depth of 2, then a boundary of 2 repetitions are done around
% the region "X" as:
%
%                 Y Y Y Y Y
%                 Y Y Y Y Y
%                 Y Y X Y Y
%                 Y Y Y Y Y
%                 Y Y Y Y Y
%
% Note: a tile depth of 0 returns simply X.
%
% The points are ordered such that the resulting matrix is Kx2, with K =
% (N*(2d+1)^2), and where d is the depth. Thus a depth of 2, which produces
% a 5x5 tiling, will produce N*25x2 matrix. The original points are in the
% middle-most portion of the matrix. Specifically, the resulting tile sets,
% S1 to SK, are organized as follows, using the depth=1 case as an example:
%
%                 S3 S6 S9
%                 S2 S5 S8
%                 S1 S4 S7
%
% Thus, S1 contains the points on rows (1..N), S2 is ((N+1)...(2N)), S3 is
% ((2N+1)...(3N)), etc.
%
% FORMAT:
%
%    [tiledPoints] = ...
%    fcn_MapGen_tilePoints(...
%    inputPoints, tileDepth, AABB, (fig_num))
%
% INPUTS:
%
%     inputPoints: An Nx2 vector of points to be tiled. The points do not
%     have to be within the AABB.
%
%     tileDepth: an integer specifying how "deep" the tiling should be.
%     For example, an integer of 2 specifies that the tiling should be 2
%     deep around the inputPoints.
%
%     AABB: the axis-aligned bounding box, in format of
%     [xmin ymin xmax ymax], specifying the boundaries of the repeating
%     pattern.
%
%     (optional inputs)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
%
% OUTPUTS:
%
%     tiledPoints: the resulting vertices of the tiled inputPoints. The
%     points are ordered such that the resulting matrix is (N*(2d+1)^2)x2,
%     where d is the tileDepth.
%
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_convertAABBtoWalls
%     fcn_MapGen_isWithinABBB
%     fcn_MapGen_snapToAABB
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_tilePoints
% for a full test suite.
%
% This function was written on 2021_07_17 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu


% REVISION HISTORY:
% 2023_02_23 by Sean Brennan
% -- first write of function
% 2023_03_13 by Sean Brennan
% -- shut off debugging
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_07_10 by Sean Brennan
% -- fixed variable names for clarity

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
        narginchk(3,4);

        % Check the inputPoints input, make sure it has 2 columns
        fcn_DebugTools_checkInputsToFunctions(...
            inputPoints, '2column_of_numbers');

        % Check the tileDepth input, make sure it has 1 column, 1 row and is
        % integer
        fcn_DebugTools_checkInputsToFunctions(...
            tileDepth, 'positive_1column_of_integers',1);

        % Check the AABB input, make sure it is '4column_of_numbers' type, with
        % 1 row
        fcn_DebugTools_checkInputsToFunctions(...
            AABB, '4column_of_numbers',1);

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
% Method: the tiling is created by calculating offsets relative to the
% center of a "super matrix". These offsets are generated by using an
% ind2sub calculation of each subset matrix relative to the center. Then,
% using these offsets, the points are calculated for the entire
% supermatrix in one step.

% Calculate the number of tiles to create, e.g. the K value
Npoints = length(inputPoints);
superMatrix_height = (2*tileDepth+1);
NumTiles = superMatrix_height^2;

if NumTiles<1 || ~isnumeric(NumTiles)
    error('Tile depth producing an unusable tiling pattern.')
end

% Convert to indices and then to rows/columns in the super-array. These
% indices count down the columns, starting with column 1. They are
% converted into row/column I and J values for the superMatrix, and are
% used later to calculate the offsets.
index_numbers = (1:NumTiles)'; 
[superMatrix_I,superMatrix_J] = ind2sub([superMatrix_height superMatrix_height],index_numbers);


% Calculate origin and range from AABB
origin_point = AABB(1,1:2);
range_xy_points  = AABB(1,3:4) - origin_point;

% The center of the super matrix will be at an offset (tileDepth+1), so we
% need to subtract that offset from the superMatrix I and J indices before
% calculating the overall offsets.
offset_superMatrix_I = superMatrix_I - (tileDepth+1);
offset_superMatrix_J = superMatrix_J - (tileDepth+1);

% Make each of the offsets one "big" offset matrix
only_row_matrix_I = reshape(offset_superMatrix_I',1,NumTiles);
only_row_matrix_J = reshape(offset_superMatrix_J',1,NumTiles);

repeated_only_row_matrix_I = repmat(only_row_matrix_I,Npoints,1);
repeated_only_row_matrix_J = repmat(only_row_matrix_J,Npoints,1);

% The integer offset super matrix is the matrix that represents, for each
% point, its "integer offset" tile count from the center of the tiling. The
% geometric center is at  integer offset (0,0). The tile to the left of
% this is at (-1,0), the one "above" is at (0,1), etc. The one that is 2
% tiles to the right is at (2,0).
BIG_integer_offset_superMatrix_I = reshape(repeated_only_row_matrix_I,Npoints*NumTiles,1);
BIG_integer_offset_superMatrix_J = reshape(repeated_only_row_matrix_J,Npoints*NumTiles,1);

% NOTE: I and J indices are in matrix form, which rows are "Y" and columns
% are "X", so we reverse the order here
BIG_integer_offset_superMatrix   = [BIG_integer_offset_superMatrix_J BIG_integer_offset_superMatrix_I];

% Convert these to offsets
BIG_offsets = BIG_integer_offset_superMatrix.*range_xy_points;

% Next, shift all the points to the AABB origin, so that range operations
% are multiplicative. In other words, we want to make points that are 2
% times away as y = 2*x if x is the point. This doesn't work if there's a
% non-zero origin
shifted_points = inputPoints - origin_point;

% We are going to need points to cover the entire tiling, so we repeat the
% points across all tiles
BIG_shifted_points = repmat(shifted_points,NumTiles,1);

% Finally, calculate the tiledPoints by transforming the points
tiledPoints = (BIG_shifted_points+BIG_offsets)+origin_point;



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
    grid on;
    grid minor;
    
    % Plot the inputPoints
    plot(inputPoints(:,1),inputPoints(:,2),'bo','Markersize',5);
    
    % Label the center
    center_of_label = origin_point + range_xy_points/2;
    label_location = center_of_label;
    text(label_location(1,1), label_location(1,2),sprintf('Original\nPoints'),'Fontsize',14,'Color',[0 0 1],'HorizontalAlignment','Center');

    % Calculate the AABB, but shifted in by some eps factor so that they don't
    % lie on top of each other when tiled
    offset = 1000*eps;
    AABB_coordinates = [...
        AABB(1,1)+offset  AABB(1,2)+offset;
        AABB(1,3)-offset  AABB(1,2)+offset;
        AABB(1,3)-offset  AABB(1,4)-offset;
        AABB(1,1)+offset  AABB(1,4)-offset;
        AABB(1,1)+offset  AABB(1,2)+offset;
        ];

    % Plot the AABB
    plot(AABB_coordinates(:,1),AABB_coordinates(:,2),'b-','Linewidth',3);

    centerTile = ceil(NumTiles/2);
    % Plot the shifted points and AABBs, each in a different color
    AABB_offsets = [offset_superMatrix_J offset_superMatrix_I].*range_xy_points;
    for ith_tile = 1:NumTiles
        range = (ith_tile-1)*Npoints + (1:Npoints)';
        fig_handle = plot(tiledPoints(range,1),tiledPoints(range,2),'.','Markersize',10);
        
        if ith_tile~=centerTile
            % Plot the AABB in the same color
            offset_AABB = AABB_coordinates + AABB_offsets(ith_tile,:);
            plot(offset_AABB(:,1),offset_AABB(:,2),'-','Linewidth',2,'Color',get(fig_handle,'Color'));
            
            % Label the tile
            label_location = center_of_label + AABB_offsets(ith_tile,:);
            text(label_location(1,1), label_location(1,2),sprintf('%.0d',ith_tile),'Fontsize',14,'Color',get(fig_handle,'Color'),'HorizontalAlignment','Center');
        end
    
    end
    
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

