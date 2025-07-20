function [dilatedMatrix, dilationIndicesNearby] = ...
    fcn_GridMapGen_dilateOccupancyByN(occupancyMatrix, dilationLevel, varargin)
% fcn_GridMapGen_dilateOccupancyByN  dilates a matrix by N cells
%
% given an N-by-M occupancyMatrix, returns an N-by-M dilatedMatrix which is
% occupancyMatrix dialated by n pixels.
%
% Method: finds adjacent indices to the occupancyMatrix and sets the values
% to 1 using ind2sub for speed.
%
% NOTE: this is NOT true dilation. For true dilation, see the image
% processing toolbox.
%
% NOTE: this method does NOT work for occupancyMatrix entries that are not
% binary. It only works if all values are 0 or 1.
%
% FORMAT:
%
%     [dilatedMatrix, dilationIndexShift] = ...
%        fcn_GridMapGen_dilateOccupancyByN(occupancyMatrix, dilationLevel, ...
%        (dilationIndexShift), (fig_num));
%
% INPUTS:
%
%     occupancyMatrix: N-by-M matrix
%
%     dilationLevel: number of cells (or pixels) to dilate by
%
%     (optional inputs)
%
%     dilationIndexShift: a precomputed index shifting matrix. Entering
%     this, if dilations are repeated, can significantly increase
%     processing speeds. See outputs listed below. If dilationIndexShift is
%     entered, the dilationLevel entry is not used.
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%     dilatedMatrix: a matrix where all elements of occupancyMatrix>=1 have
%     been dilated outward by n pixels
%
%     dilationIndexShift: the px(NxM) matrix that is used to calculate
%     which p points nearby each of the NxM entries of the matrix need to
%     be set to 1, for a given dilation level
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     occupancyMatrix = rand(100,100)<0.1; % about 10 percent occupied
%     percentOccupied = fcn_GridMapGen_dilateOccupancyStats(occupancyMatrix*1.0, (28282));
%     dilationLevel = 1;
%
%     % Call the function
%     dilatedMatrix = fcn_GridMapGen_dilateOccupancyByN(occupancyMatrix, dilationLevel, (fig_num));
%
% See the script: script_test_fcn_GridMapGen_dilateOccupancyByN
% for a full test suite.
%
% This function was written on 2008_10_18 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

%
% REVISION HISTORY:
%
% 2008_10_18 by Sean Brennan
% -- first write of function
% 2025_07_17 by Sean Brennan
% -- imported into MapGen library with updates to formatting

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 4; % The largest Number of argument inputs to the function
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

        % Check the occupancyMatrix input, must be [2+ 2+] in size
        fcn_DebugTools_checkInputsToFunctions(occupancyMatrix*1.0, 'positive_2orMorecolumn_of_numbers',[2 3]);

        % Check the dilationLevel input, must be [1 1] integer
        fcn_DebugTools_checkInputsToFunctions(dilationLevel, 'positive_1column_of_integers',[1 1]);

    end
end

% check variable argument mColumns
flag_calculateIndexShift = 1;
if 3 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        dilationIndicesNearby = temp;
        flag_calculateIndexShift = 0;
        if flag_check_inputs
            % Check the inputs
            nRows = size(occupancyMatrix,1);
            mColumns = size(occupancyMatrix,2);
            assert(isequal(size(dilationIndicesNearby,2),nRows*mColumns));
            fcn_DebugTools_checkInputsToFunctions(dilationIndicesNearby, 'positive_');
        end
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

% Do the multipliers need to be calculated?
if 1==flag_calculateIndexShift
    [n,m] = size(occupancyMatrix);
    dilationIndicesNearby = fcn_INTERNAL_calcNearbyIndices(n,m, dilationLevel);
end

% Perform the dilation
indOnes=find(occupancyMatrix);
whichIndiciesToFillWithOnes = dilationIndicesNearby(:,indOnes);

% Two methods to remove outliers: (no longer needed because we get rid of
% that in the initialize code of the dialation matrix).
% # 1
% indDialatedPoints=unique(indDialatedPoints);
% if indDialatedPoints(1,1)==0 %fiter the 0's
%    indDialatedPoints=indDialatedPoints(2:end);
% end

% # 2
% indDialatedPoints=indDialatedPoints(find(indDialatedPoints));
dilatedMatrix = occupancyMatrix;
dilatedMatrix(whichIndiciesToFillWithOnes(whichIndiciesToFillWithOnes>0))=ones;


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

    % Display dialated on top of original by making testOccupancy
    % pixels have value of 3, dilated have values of 2, 1 otherwise
    % Then use a colormap to map 3 values to white, 2 to grey, 1 to
    % white
    image(2*(occupancyMatrix>0) + 1.0*(dilatedMatrix>0) + 1);
    colormap([1 1 1;0.5 0.5 0.5;0 0 0])

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

%% fcn_INTERNAL_calcNearbyIndices
function dilationIndicesNearby = fcn_INTERNAL_calcNearbyIndices(n,m, dilationLevel)

% Create a matrix of indices where each entry is the actual index value
% in the subscript to ind form. In other words, each entry is its row
% value added with (columnNumber - 1 ) * rowDepth
% NOTE: this works in MATLAB. In Python, the ind mode is by columns
% first, then rows. So the logic would need to change here.
matrixOfInd=zeros(n,m);
for ith_col=1:n
    for jth_row=1:m
        matrixOfInd(ith_col,jth_row)=ith_col+(jth_row-1)*n;
    end
end

% Next, for each ith/jth element, find which nearby elements are within
% the dilation distance. Note that dilationLevel of X will only affect
% entries X distance from the center, hence we only need to search from
% -dilationLevel to +dilationLevel. The maximum possible elements is
% (2*dilationLevel + 1)^2. We create a mask that represents the
% difference between the center index and nearby indices that are
% within the radius. The indexShift matrix thus contains the offsets
% from the current index, in each direction, that can be used to
% calculate which OTHER indices need to be changed.
maxNearby = (2*dilationLevel + 1);
indexShift = zeros(maxNearby,maxNearby);
dilationLevels = -dilationLevel:dilationLevel;
centerIndex = find(dilationLevels==0);
for ith_col = 1:length(dilationLevels)
    ith_dist = dilationLevels(ith_col);
    for jth_row = 1:length(dilationLevels)
        jth_dist = dilationLevels(jth_row);

        if (ith_dist^2 + jth_dist^2)<=dilationLevel^2
            indexShift(jth_row, ith_col) = (ith_col-centerIndex)*n + (jth_row-centerIndex);
        end
    end
end

% allShiftedIndices is a kx1 matrix of index offsets, relative to
% current index
allShiftedIndices = indexShift(indexShift~=0);
KshiftedIndices = length(allShiftedIndices);

% For each index in the resulting array, list which OTHER pixels this
% index will change, using the allShiftedIndices listing
allIndicesCountingUpVector = (1:n*m);
% Copy allIndicesCountingUpVector down by rows to make same length as
% allShiftedIndices. This produces a kx(n*m) matrix where each column
% counts up 1, 2, 3, etc.
allIndicesCountingUpMatrix = repmat(allIndicesCountingUpVector,KshiftedIndices,1);

% Fill in offsets for each index value. This produces a k x (n*m)
% matrix.
offsetsToAdd = repmat(allShiftedIndices,1,n*m);

% Loop through all the indices near edges, and find offsets for each.
% Fix boundary offsets to avoid wrap-around.

% Fix rows
edgeValuesRows = [1:centerIndex, n-centerIndex+1:n];
edgeValuesRows = unique(edgeValuesRows);

for ith_row_index = 1:length(edgeValuesRows)
    ith_row = edgeValuesRows(ith_row_index);
    for jth_col = 1:m
        thisIndex = sub2ind([n m],ith_row,jth_col);
        thisOffsets = offsetsToAdd(:,thisIndex);

        lowRowsValid = 1;
        highRowsValid = maxNearby;

        proximityToRowStart = ith_row;
        if proximityToRowStart<centerIndex
            % The following crops rows that are low
            lowRowsValid = centerIndex-proximityToRowStart+1;
        end
        proximityToRowEnd = n-ith_row+1;
        if proximityToRowEnd<centerIndex
            % The following crops rows that are high
            highRowsValid = centerIndex+proximityToRowEnd-1;
        end

        lowColsValid = 1;
        highColsValid = maxNearby;

        proximityToColStart = jth_col;
        if proximityToColStart<centerIndex
            % The following crops columns that are low
            lowColsValid = centerIndex-proximityToColStart+1;
        end
        proximityToColEnd = m-jth_col+1;
        if proximityToColEnd<centerIndex
            % The following crops columns that are high
            highColsValid = centerIndex+proximityToColEnd-1;
        end

        maskPortionToUse = indexShift(lowRowsValid:highRowsValid,lowColsValid:highColsValid);
        validValuesInMask = maskPortionToUse(maskPortionToUse~=0);

        % Returns values in thisOffsets that are NOT in
        % validValuesInMask
        [~,badIndices] = setdiff(thisOffsets,validValuesInMask);
        goodVector = thisOffsets;
        goodVector(badIndices,1) = -inf;
        offsetsToAdd(:,thisIndex) = goodVector;

    end
end


% Fix columns
edgeValuesCols = [1:centerIndex, m-centerIndex+1:m];
edgeValuesCols = unique(edgeValuesCols);

for jth_col_index = 1:length(edgeValuesCols)
    jth_col = edgeValuesCols(jth_col_index);
    for ith_row = 1:n
        thisIndex = sub2ind([n m],ith_row,jth_col);
        thisOffsets = offsetsToAdd(:,thisIndex);

        lowRowsValid = 1;
        highRowsValid = maxNearby;

        proximityToRowStart = ith_row;
        if proximityToRowStart<centerIndex
            % The following crops rows that are low
            lowRowsValid = centerIndex-proximityToRowStart+1;
        end
        proximityToRowEnd = n-ith_row+1;
        if proximityToRowEnd<centerIndex
            % The following crops rows that are high
            highRowsValid = centerIndex+proximityToRowEnd-1;
        end

        lowColsValid = 1;
        highColsValid = maxNearby;

        proximityToColStart = jth_col;
        if proximityToColStart<centerIndex
            % The following crops columns that are low
            lowColsValid = centerIndex-proximityToColStart+1;
        end
        proximityToColEnd = m-jth_col+1;
        if proximityToColEnd<centerIndex
            % The following crops columns that are high
            highColsValid = centerIndex+proximityToColEnd-1;
        end

        maskPortionToUse = indexShift(lowRowsValid:highRowsValid,lowColsValid:highColsValid);
        validValuesInMask = maskPortionToUse(maskPortionToUse~=0);

        % Returns values in thisOffsets that are NOT in
        % validValuesInMask
        [~,badIndices] = setdiff(thisOffsets,validValuesInMask);
        goodVector = thisOffsets;
        goodVector(badIndices,1) = -inf;
        offsetsToAdd(:,thisIndex) = goodVector;

    end
end

dilationIndicesNearby= allIndicesCountingUpMatrix + offsetsToAdd;

% Allowable indices cannot be less than 1 or greater than n*m
dilationIndicesNearby(dilationIndicesNearby<1) = 0;
dilationIndicesNearby(dilationIndicesNearby>n*m) = 0;

end % Ends fcn_INTERNAL_calcNearbyIndices


