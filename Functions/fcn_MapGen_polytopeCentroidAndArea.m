function [ ...
Centroid, ...
Area ...
] = ...
fcn_MapGen_polytopeCentroidAndArea( ...
vertices, ...
varargin...
)
% fcn_MapGen_polytopeCentroidAndArea
% calculates the centroid and area of a closed polytope.
%
%
%
% FORMAT:
%
%    [ ...
%    Centroid, ...
%    Area ...
%    ] = ...
%    fcn_MapGen_polytopeCentroidAndArea( ...
%    vertices, ...
%    (fig_num) ...
%    )
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
%     Centroid: the calculated centroid of the polytope, given as
%     [x-coordinate y_coordinate]
%
%     Area: the unsigned area enclosed by the polytope
%
%
% DEPENDENCIES:
%
%     fcn_MapGen_checkInputsToFunctions
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
% sequence, causing Area to be zero. Get a divide-by-zero problem. This is
% fixed now via an if-statement check which uses the mean of points if A=0.

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
if length(flag_check_inputs) == 0
    flag_check_inputs = '1';
end
if length(flag_do_plot) == 0
    flag_do_plot = '1';
end
if length(flag_do_debug) == 0
    flag_do_debug = '1';
end

% convert flag from char string to logical
flag_check_inputs = flag_check_inputs == '1';
flag_do_plot = flag_do_plot == '1';
flag_do_debug = flag_do_debug == '1';

if flag_do_debug
    fig_for_debug = 60;
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
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end

    % Check the vertices input, make sure it is '2column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        vertices, '2column_of_numbers');

end

% Does user want to show the plots?
if  2== nargin
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

% current points
xi = vertices(1:end-1,1);
yi = vertices(1:end-1,2);

% next points
xip1 = vertices(2:end,1);
yip1 = vertices(2:end,2);

% signed area
A = sum(xi.*yip1 - xip1.*yi)/2;

% Centroid calculation
if A>0
    Cx = sum((xi+xip1).*(xi.*yip1 - xip1.*yi))/(6*A); % centroid x coordinate
    Cy = sum((yi+yip1).*(xi.*yip1 - xip1.*yi))/(6*A); % centroid x coordinate
    Centroid = [Cx, Cy];
else
    Centroid = [mean(xi) mean(yi)];
end

Area = abs(A); % unsigned area


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
    figure(fig_num)
    clf;
    hold on

    plot(vertices(:,1),vertices(:,2),'b-','linewidth',2);

    plot(Cx,Cy,'go','Markersize',10);

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


