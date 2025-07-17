function shrunkPolytope = fcn_MapGen_polytopeShrinkToRadius(shrinker, newRadius, varargin)
% fcn_MapGen_polytopesShrinkToRadius shrinks the polytopes to achieve the
% specified maximum radius. The vertices are all porportionally pulled
% toward the centroid location such that the new maximum radius matches the
% specified maximum radius.
%
% FORMAT:
%     shrunkPolytope = fcn_MapGen_polytopeShrinkToRadius(shrinker, newRadius, (fig_num))
%
% INPUTS:
%
%     shrinker: original polytope with same fields as shrunkPolytopes
%     below
%
%     newRadius: desired polytope radius
%     
%    (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%     shrunkPolytope: a polytope structure containing the shrunken input
%     polytope
%   
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_polytopesFillFieldsFromVertices
%     fcn_MapGen_plotPolytopes
% 
% EXAMPLES:
%      
% For additional examples, see: script_test_fcn_MapGen_polytopeShrinkToRadius
%
% This function was written on 2019_08_29 by Seth Tau
% Questions or comments? sat5340@psu.edu 
%

% Revision History:
% 2021-06-08 - S. Brennan
% -- revised function to prep for MapGen class 
% -- added plotting option
% -- added comments, added debugging option
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% -- fixed call to fcn_MapGen_polytopesFillFieldsFromVertices
% 2025_07_16 by Sean Brennan
% -- cleaned up documentation and typos
% -- removed tolerance input as it wasn't being used
% -- turned on fast mode for all internal calls to external functions
% 2025_07_17 by Sean Brennan
% -- standardized Debugging and Input checks area, Inputs area
% -- made codes use MAX_NARGIN definition at top of code, narginchk
% -- made plotting flag_do_plots and code consistent across all functions

% TO DO
% -- Vectorize the for loop if possible
% -- check inputs are positive numbers where appropriate (e.g. make a
%    % "positive number" chec

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 3; % The largest Number of argument inputs to the function
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

        % Check the shrinker input
        fcn_DebugTools_checkInputsToFunctions(shrinker, 'polytopes');

        % Check the newRadius input
        fcn_DebugTools_checkInputsToFunctions(newRadius, 'positive_1column_of_numbers',1);

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

% Initialize the result:
shrunkPolytope = shrinker;

% pull values
vertices = shrinker.vertices;
centroid = shrinker.mean;
rad = shrinker.max_radius;

% determine scale factor
scale = newRadius/rad;

% calculation error can sometimes make the scale greater than 1, so if we
% are doing shrinking, check that actual shrinking is happening!
if scale < 1 
    % find new vertices
    new_vert = centroid + scale*(vertices-centroid);
    shrunkPolytope.vertices = new_vert;
end

% fill in other fields from the vertices field
shrunkPolytope = fcn_MapGen_polytopesFillFieldsFromVertices(shrunkPolytope, [], -1);


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
    
    % Plot the cetroid in black
    plot(centroid(:,1),centroid(:,2),'ko','Markersize',10);
    
    % Plot the input shrinker in red
    % fcn_MapGen_plotPolytopes(shrinker,fig_num,'r',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [1 0 0];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(shrinker, (plotFormat), (fillFormat), (fig_num)); %#ok<NASGU>
    
    % plot the output polytope in blue
    % fcn_MapGen_OLD_plotPolytopes(shrunkPolytope,fig_num,'b',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [0 0 1];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(shrunkPolytope, (plotFormat), (fillFormat), (fig_num)); %#ok<NASGU>

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function




