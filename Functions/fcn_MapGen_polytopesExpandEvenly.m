function expandedPolytopes  = fcn_MapGen_polytopesExpandEvenly( polytopes, expansionDistance, varargin)
% fcn_MapGen_polytopesExpandEvenly
% Expands an obstacle out by expansionDistance on all sides.  This function works as intended
% with convex polytopes.  
% 
% NOTE: There can be counterintuitive behavior with non-convex polytopes
% that can be avoided by using fcn_MapGen_polytopesExpandEvenlyForConcave
% which implements MATLAB's polyshape object and polybuffer method (object
% function).
%
% FORMAT:
%
%     expandedPolytopes  = fcn_MapGen_polytopesExpandEvenly( polytopes, expansionDistance,(fig_num))
%
% INPUTS:
%
%     polytopes: the structure of 'polytopes' type that stores the
%     polytopes to be expanded
%
%     expansionDistance: distance to expand the obstacle
%
%     (optional inputs)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%     expandedPolytopes: structure of expanded polytopes
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_fillPolytopeFieldsFromVertices
%     fcn_MapGen_plotPolytopes
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_polytopesExpandEvenly
% for a full test suite.
%
% This function was written on 2018_11_17, Adjusted example code on 2021_04_28 by Seth Tau, Rebased on 2021_06_26 by S. Brennan by Seth Tau
% Questions or comments? contact sbrennan@psu.edu and sat5340@psu.edu

%
% REVISION HISTORY:
%
% 2018_11_17, Seth Tau
% -- first write of script
% 2021_04_28, Seth Tau
% -- Adjusted example code ,
% 2021_06_26 S. Brennan
% -- Rebased code
% -- Rewrote for clarity
% 2021_07_06 S. Brennan
% -- Vectorized plotting into array structure, to better support legends
% (rather than plotting all polytopes individually)
% 2024_02_14 S.J. Harnett
% -- Updated docstring comment to point to fcn_MapGen_polytopesExpandEvenlyForConcave
%    for non-convex obstacles
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% -- fixed call to fcn_MapGen_fillPolytopeFieldsFromVertices
% 2025_07_16 by Sean Brennan
% -- cleaned up header comments
% -- replaced exp_dist with expansionDistance
% -- changed plotPolytopes to new format
% -- turned on fast mode on sub-function calls

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
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
        narginchk(2,3);

        % Check the polytopes input, make sure it is 'polytopes' type
        fcn_DebugTools_checkInputsToFunctions(...
            polytopes, 'polytopes');

        % Check the expansionDistance input, make sure it is 'positive_column_of_numbers' type
        fcn_DebugTools_checkInputsToFunctions(...
            expansionDistance, 'positive_1column_of_numbers',1);

    end
end


% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  3 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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

expandedPolytopes = polytopes; % both structures will be the same size

for ith_poly = 1:size(polytopes,2) % check each obstacle

    % pull values
    vertices = polytopes(ith_poly).vertices;
    centroid = polytopes(ith_poly).mean;
    rad = polytopes(ith_poly).max_radius;

    % Calculate scale
    scale = (rad+expansionDistance)/rad;

    % Calculate new vertices
    expandedPolytopes(ith_poly).vertices = centroid + scale*(vertices-centroid);

    % fill in other fields from the vertices field
    expandedPolytopes(ith_poly) = fcn_MapGen_fillPolytopeFieldsFromVertices(expandedPolytopes(ith_poly),[], -1);

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
    figure(fig_num)
    clf;

    % LineWidth = 2;
    % fcn_MapGen_OLD_plotPolytopes(polytopes,fig_num,'r-',LineWidth);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [1 0 0];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat), (fillFormat), (fig_num)); 
    set(h_plot,'DisplayName','polytopes')

    % fcn_MapGen_OLD_plotPolytopes(expandedPolytopes,fig_num,'b-',LineWidth,'square');
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [0 0 1];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(expandedPolytopes, (plotFormat), (fillFormat), (fig_num)); 
    set(h_plot,'DisplayName','expandedPolytopes')

    legend('Original','Expanded')
    box on
    xlabel('X Position')
    ylabel('Y Position')

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

