function [tiled_points] = ...
   fcn_MapGen_tilePoints(...
   input_points,tile_depth,AABB, varargin)

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
%    [tiled_points] = ...
%    fcn_MapGen_tilePoints(...
%    input_points,tile_depth,AABB)
%
% INPUTS:
%
%     input_points: An Nx2 vector of points to be tiled. The points do not
%     have to be within the AABB.
%
%     tile_depth: an integer specifying how "deep" the tiling should be.
%     For example, an integer of 2 specifies that the tiling should be 2
%     deep around the input_points.
%
%     AABB: the axis-aligned bounding box, in format of
%     [xmin ymin xmax ymax], specifying the boundaries of the repeating
%     pattern.
%
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a
%     figure to be drawn showing results.
%
%
% OUTPUTS:
%
%     tiled_points: the resulting vertices of the tiled input_points. The
%     points are ordered such that the resulting matrix is (N*(2d+1)^2)x2,
%     where d is the tile_depth.
%
%
% DEPENDENCIES:
%
%     fcn_MapGen_checkInputsToFunctions
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

%
% REVISION HISTORY:
%
% 2023_02_23 by Sean Brennan
% -- first write of function
% 2023_03_13 by Sean Brennan
% -- shut off debugging
%
% TO DO:
%


%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plots = 0;      % % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    fig_for_debug = 846; 
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
    narginchk(3,4);
    
    % Check the input_points input, make sure it has 2 columns
    fcn_DebugTools_checkInputsToFunctions(...
        input_points, '2column_of_numbers');
    
    % Check the tile_depth input, make sure it has 1 column, 1 row and is
    % integer
    fcn_DebugTools_checkInputsToFunctions(...
        tile_depth, 'positive_1column_of_integers',1);
       
    % Check the AABB input, make sure it is '4column_of_numbers' type, with
    % 1 row
    fcn_DebugTools_checkInputsToFunctions(...
        AABB, '4column_of_numbers',1);
    
end

% Does user want to show the plots?
if  4== nargin
    fig_num = varargin{end};
    if ~isempty(fig_num)
        flag_do_plots = 1;
    end
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
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
% Method: the tiling is created by calculating offsets relative to the
% center of a "super matrix". These offsets are generated by using an
% ind2sub calculation of each subset matrix relative to the center. Then,
% using these offsets, the points are calculated for the entire
% supermatrix in one step.

% Calculate the number of tiles to create, e.g. the K value
Npoints = length(input_points);
superMatrix_height = (2*tile_depth+1);
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

% The center of the super matrix will be at an offset (tile_depth+1), so we
% need to subtract that offset from the superMatrix I and J indices before
% calculating the overall offsets.
offset_superMatrix_I = superMatrix_I - (tile_depth+1);
offset_superMatrix_J = superMatrix_J - (tile_depth+1);

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
shifted_points = input_points - origin_point;

% We are going to need points to cover the entire tiling, so we repeat the
% points across all tiles
BIG_shifted_points = repmat(shifted_points,NumTiles,1);

% Finally, calculate the tiled_points by transforming the points
tiled_points = (BIG_shifted_points+BIG_offsets)+origin_point;



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
    grid on;
    grid minor;
    
    % Plot the input_points
    plot(input_points(:,1),input_points(:,2),'bo','Markersize',5);
    
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
        fig_handle = plot(tiled_points(range,1),tiled_points(range,2),'.','Markersize',10);
        
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

