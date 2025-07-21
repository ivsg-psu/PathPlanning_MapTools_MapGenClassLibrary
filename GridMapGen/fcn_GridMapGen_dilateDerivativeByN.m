function [rowDerivative, colDerivative, leftDilationMultiplier, rightDilationMultiplier] = ...
    fcn_GridMapGen_dilateDerivativeByN(inputMatrix, derivativeOrder, varargin)
% fcn_GridMapGen_dilateDerivativeByN  dilates a matrix by N cells
% 
% given an N-by-M inputMatrix, returns the m-th derivatives by dilating
% a 2nd order derivative appxoiamtion by n pixels.
%
% Method: performs diagonal plus off-diagonal left- and right- matrix
% multiplication of the derivative operator raised to powers of the derivativeOrder,
% thereby giving the derivative approximations at increasing powers
% 
% FORMAT:
% 
%     [rowDerivative, colDerivative, fullDerivative, leftDilationMultiplier, rightDilationMultiplier] = ...
%        fcn_GridMapGen_dilateDerivativeByN(inputMatrix, derivativeOrder, ...
%        (leftDilationMultiplier), (rightDilationMultiplier), (fig_num));
% 
% INPUTS:
% 
%     inputMatrix: N-by-M matrix
%   
%     derivativeOrder: number of cells (or pixels) to dilate by
% 
%     (optional inputs)
%
%     leftDilationMultiplier: an N-by-N matrix that is precalculated from
%     prior calls to the function, passed in to speed up the code. This
%     input supercedes any dilation level values. See outputs listed below.
%
%     rightDilationMultiplier: an M-by-M matrix that is precalculated from
%     prior calls to the function, passed in to speed up the code. This
%     input supercedes any dilation level values. See outputs listed below.
%     NOTE: both left and right multipliers must be entered for them to be
%     used. Otherwise, the matrix is recalculated.
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
% 
% OUTPUTS:
% 
%     rowDerivative: a matrix of the derivatives in the row direction
%
%     colDerivative: a matrix of the derivatives in the col direction
%
%     leftDilationMultiplier: the NxN left multiplier of the
%     inputMatrix to produce the fullDerivative. This is output so that,
%     if future dilations are required for the map at the same dilation
%     level, this calculation does not need to be repeated, as this
%     calculation is a significant portion of the time required for the
%     function.
%
%     rightDilationMultiplier: the MxM right multiplier of the
%     inputMatrix to produce the fullDerivative. This is output so that,
%     if future dilations are required for the map at the same dilation
%     level, this calculation does not need to be repeated, as this
%     calculation is a significant portion of the time required for the
%     function.
% 
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
% 
% EXAMPLES:
%
%     inputMatrix = rand(100,100)<0.1; % about 10 percent occupied
%     percentOccupied = fcn_GridMapGen_dilateOccupancyStats(inputMatrix*1.0, (28282));
%     derivativeOrder = 1;
%     
%     % Call the function
%     fullDerivative = fcn_GridMapGen_dilateDerivativeByN(inputMatrix, derivativeOrder, (fig_num));
%
% See the script: script_test_fcn_GridMapGen_dilateDerivativeByN
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
        narginchk(2,MAX_NARGIN);

        % Check the inputMatrix input, must be [2+ 2+] in size
        fcn_DebugTools_checkInputsToFunctions(inputMatrix*1.0, 'positive_2orMorecolumn_of_numbers',[2 3]);

        % Check the derivativeOrder input, must be [1 1] integer
        fcn_DebugTools_checkInputsToFunctions(derivativeOrder, 'positive_1column_of_integers',[1 1]);
        
    end
end

% check variable argument mColumns
flag_calculateMultipliers = 1;
if 4 <= nargin
    temp1 = varargin{1};
    temp2 = varargin{2};
    if flag_check_inputs
        if ~isempty(temp1) && isempty(temp2)
            error('The leftDilationMultiplier is specified but rightDilationMultiplier was not. Code will not continue with just one entry.')
        end
        if ~isempty(temp2) && isempty(temp1)
            error('The rightDilationMultiplier is specified but leftDilationMultiplier was not. Code will not continue with just one entry.')
        end
    end
    if ~isempty(temp1) && ~isempty(temp2)
        leftDilationMultiplier = temp1;
        rightDilationMultiplier = temp2;
        flag_calculateMultipliers = 0;
        if flag_check_inputs
            % Check the inputs
            nRows = size(inputMatrix,1);
            mColumns = size(inputMatrix,2);
            assert(isequal(size(leftDilationMultiplier),[nRows nRows]));
            fcn_DebugTools_checkInputsToFunctions(leftDilationMultiplier, 'positive_');
            assert(isequal(size(rightDilationMultiplier),[mColumns mColumns]));
            fcn_DebugTools_checkInputsToFunctions(rightDilationMultiplier, 'positive_');

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
if 1==flag_calculateMultipliers
    [n,m] = size(inputMatrix);
    if 1==0 % Change to n==m if need speed for square matrices
    
        derivativeOperator = fcn_INTERNAL_fillOpeator(m);
        
        leftDilationMultiplier  = derivativeOperator^derivativeOrder;
        rightDilationMultiplier = leftDilationMultiplier;
    else
        leftOperator  = fcn_INTERNAL_fillOpeator(n);
        rightOperator = fcn_INTERNAL_fillOpeator(m);

        leftDilationMultiplier    = leftOperator^derivativeOrder;
        rightDilationMultiplier   = rightOperator^derivativeOrder;
    end
end

% fullDerivative = leftDilationMultiplier*inputMatrix*rightDilationMultiplier';
rowDerivative  = leftDilationMultiplier*inputMatrix;
colDerivative  = inputMatrix*rightDilationMultiplier';

% Last column will not be correct because not possible to find derivative
% there
colDerivative(:,end)  = colDerivative(:,end-1);
% Last row will not be correct because not possible to find derivative
% there
rowDerivative(end,:)  = rowDerivative(end-1,:);



colorMin = min([rowDerivative colDerivative],[],"all");
colorMax = max([rowDerivative colDerivative],[],"all");


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

    numColors = 256;
    cmap = parula(numColors);

    subplot(2,2,1);
    colorizedImage = floor(rescale(inputMatrix,1,numColors));
    image(colorizedImage);
    colormap(gca,cmap);
    title('inputMatrix','FontSize',10);

    
    subplot(2,2,2);
    colorizedImage = fcn_INTERNAL_rescaleToColors(rowDerivative, colorMin, colorMax, numColors);
    image(colorizedImage);
    colormap(gca,cmap);
    title('rowDerivative','FontSize',10);

    subplot(2,2,3);
    colorizedImage = fcn_INTERNAL_rescaleToColors(colDerivative, colorMin, colorMax, numColors);
    image(colorizedImage);
    colormap(gca,cmap);
    title('colDerivative','FontSize',10);

    if 1==0
        % Create Sobel kernels
        sobel_horizontal = fspecial('sobel');  % Emphasizes horizontal edges
        sobel_vertical = sobel_horizontal';    % Transpose to emphasize vertical edges

        % Apply the filters using convolution
        Gx = imfilter(inputMatrix, sobel_horizontal);
        Gy = imfilter(inputMatrix, sobel_vertical);

        subplot(2,2,4);
        colorizedImage = fcn_INTERNAL_rescaleToColors(-Gy, colorMin*7, colorMax*7, numColors);
        image(colorizedImage);
        colormap(gca,cmap);
        title('Gx','FontSize',10);
    end

    % subplot(2,2,4);
    % colorizedImage = fcn_INTERNAL_rescaleToColors(fullDerivative, min(fullDerivative,[],'all'), max(fullDerivative,[],'all'), numColors);
    % image(colorizedImage);
    % colormap(gca,cmap);
    % title('fullDerivative','FontSize',10);


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

%% fcn_INTERNAL_fillOpeator
function derivativeOperator = fcn_INTERNAL_fillOpeator(m)
derivativeOperator = -diag(ones(m,1),0) +diag(ones(m-1,1),1);

% derivativeOperator = diag(0.5*ones(m-1,1),1) -diag(0.5*ones(m-1,1),-1);

%derivativeOperator(1,1)   = -0.5;
% derivativeOperator(1,2)   =  1;
% derivativeOperator(m,m-1) = -1;
derivativeOperator(m,m)   =  0;
end % Ends fcn_INTERNAL_fillOpeator


%% fcn_INTERNAL_rescaleToColors
function rescaledToColorIntegers = fcn_INTERNAL_rescaleToColors(randomMatrixDilated, colorMin, colorMax, numColors)

colorInterval = (colorMax-colorMin)/(numColors-1); % The interval represents the "jump" between different colors
offsetRemovedMatrix = randomMatrixDilated - colorMin; % Remove the offset
countMatrix = floor(offsetRemovedMatrix./colorInterval)+1;
rescaledToColorIntegers = min(max(countMatrix,1),numColors);

end % Ends fcn_INTERNAL_rescaleToColors
