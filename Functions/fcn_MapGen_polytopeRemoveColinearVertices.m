function cleanedVertices = fcn_MapGen_polytopeRemoveColinearVertices(inputVertices, varargin)

% Given a set of vertices that may have multiple points in a line on the
% same edge, this function returns an equivalent polytope with no repeats
% on any edge.
%
% FORMAT:
% 
%     cleanedVertices = fcn_MapGen_polytopeRemoveColinearVertices(inputVertices, (fig_num))
%
% INPUTS:
%
%     inputVertices: an Nx2 matrix of [x y] vertices
%
%     (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%     cleanedVertices: an Mx2 matrix of [x y] points, where no vertices
%     are repeated
%   
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
% 
% EXAMPLES:
%      
% For additional examples, see:
% script_test_fcn_MapGen_polytopeRemoveColinearVertices
%
% This function was written on 2021_08_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu 
%

% Revision History:
% 2021_08_03 - S. Brennan
% -- first write of code
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
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

        % Check the vertices input, 2 columns
        fcn_DebugTools_checkInputsToFunctions(inputVertices, '2column_of_numbers');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Remove repeats
[uniqueVertices,~,~] = unique(round(inputVertices,7),'rows','stable');

% Use the cross-product to eliminate co-linear points
Npoints = length(uniqueVertices(:,1));
crossResults = nan(Npoints,1);
for ith_point = 1:Npoints
    
    if ith_point>1
        previous_vector = [uniqueVertices(ith_point,:)-uniqueVertices(ith_point-1,:), 0];
    else % must wrap around
        previous_vector = [uniqueVertices(end,:)-uniqueVertices(1,:), 0];
    end
    
    if ith_point<Npoints
        subsequent_vector = [uniqueVertices(ith_point+1,:)-uniqueVertices(ith_point,:), 0];
    else % must wrap around
        subsequent_vector = [uniqueVertices(1,:)-uniqueVertices(Npoints,:), 0];
    end
    
    cross_result = cross(previous_vector,subsequent_vector);
    crossResults(ith_point,1) = cross_result(1,3);
    
end

% Keep only points that "bend", e.g. have a non-zero cross-product
tolerance = 1E-7; % Roughly floating point error
keepIndices = abs(crossResults)>tolerance;

% Final polytope
cleanedVertices = uniqueVertices(keepIndices,:);



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

if flag_do_plots
    figure(fig_num);

    hold on
    axis equal
    
    % Plot the input vertices
    plot(inputVertices(:,1),inputVertices(:,2),'r.-','MarkerSize',20,'DisplayName','inputVerticies');
    
    % Plot the cropped_vertices
    plot(cleanedVertices(:,1),cleanedVertices(:,2),'ko','MarkerSize',10,'DisplayName','cleanedVertices');

    legend('Interpreter','none','Location','best');

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function


