function [centroid, area] = fcn_MapGen_polytopeCentroidAndArea( vertices, varargin)
% fcn_MapGen_polytopeCentroidAndArea
% calculates the centroid and area of a closed polytope.
% 
% FORMAT:
% 
%    [centroid, area] = fcn_MapGen_polytopeCentroidAndArea( vertices, (fig_num))
% 
% INPUTS:
% 
%     vertices: the list of verticies used to perform calculation, in 
%     format [x y] where x and y are column vectors. X: x coordinates of 
%     the polytope (with the same first and last point)  Y: y coordinates 
%     of the polytope (with the same first and last point)
% 
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a 
%     figure to be drawn showing results.
% 
% 
% OUTPUTS:
% 
%     centroid: the calculated centroid of the polytope, given as 
%     [x-coordinate y_coordinate]
% 
%     area: the unsigned area enclosed by the polytope
% 
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
% 
% 
% EXAMPLES:
% 
% See the script: script_test_fcn_MapGen_polytopeCentroidAndArea
% for a full test suite.
% 
% This function was written on 2021_07_02 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

% 
% REVISION HISTORY:
% 
% 2021_02_23 by Seth Tau
% -- Added comments
% 2021_03_02 by Seth Tau
% -- Removed old add path stuff on 
% 2021_07_02
% -- Cleaned up arguments a bit to compactify x,y coordinate convention
% -- rebased code to MapGen format
% 2022_02_17 - S. Brennan
% -- Fixed bug when someone passes a line segment or repeated point
% sequence, causing area to be zero. Get a divide-by-zero problem. This is
% fixed now via an if-statement check which uses the mean of points if A=0.
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_07_16 by Sean Brennan
% -- cleaned up header
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
MAX_NARGIN = 2; % The largest Number of argument inputs to the function
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
        narginchk(1,MAX_NARGIN);

        % Check the vertices input, make sure it has 2 columns
        fcn_DebugTools_checkInputsToFunctions(vertices, '2column_of_numbers');

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

% current points
xi = vertices(1:end-1,1); 
yi = vertices(1:end-1,2);

% next points
xip1 = vertices(2:end,1); 
yip1 = vertices(2:end,2);

% signed area
A = sum(xi.*yip1 - xip1.*yi)/2; 

% centroid calculation
if A>0
    Cx = sum((xi+xip1).*(xi.*yip1 - xip1.*yi))/(6*A); % centroid x coordinate
    Cy = sum((yi+yip1).*(xi.*yip1 - xip1.*yi))/(6*A); % centroid x coordinate
    centroid = [Cx, Cy];
else
    centroid = [mean(xi) mean(yi)];
end

area = abs(A); % unsigned area


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
    figure(fig_num)
    clf;
    hold on
    
    plot(vertices(:,1),vertices(:,2),'b-','linewidth',2);
    
    plot(centroid(:,1),centroid(:,2),'go','Markersize',10);
    
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

    
