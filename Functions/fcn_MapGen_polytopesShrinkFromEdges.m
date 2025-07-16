function [shrunkPolytopes] = ...
    fcn_MapGen_polytopesShrinkFromEdges(...
    polytopes,...
    des_gap_size,...
    varargin)
% fcn_MapGen_polytopesShrinkFromEdges shrinks the polytopes to achieve the
% desired gap size between polytopes
%
% FORMAT:
%
% [shrunkPolytopes] = ...
%     fcn_MapGen_polytopesShrinkFromEdges(...
%     polytopes,...
%     des_gap_size,...
%     (fig_num))
%
% INPUTS:
%
%     POLYTOPES: original polytopes with same fields as shrunkPolytopes
%
%     DES_GAP_SIZE: desired normal gap size between polytopes, note this is
%     relative, not absolute e.g. if the input polytopes have no gap and a
%     des_gap_size of 0.005 km is input, the resulting gap size will be
%     0.005 km.  However if the input polytopes already have a gap size of
%     0.002 km the resulting gap size will be 0.007 km for a des_gap_size
%     of 0.005 km.
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
%     shrunkPolytopes: a 1-by-n seven field structure of shrunken polytopes,
%     where n <= number of polytopes with fields:
%       vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%         the number of the individual polytope vertices
%       xv: a 1-by-m vector of vertice x-coordinates
%       yv: a 1-by-m vector of vertice y-coordinates
%       distances: a 1-by-m vector of perimeter distances from one point to the
%         next point, distances(i) = distance from vertices(i) to vertices(i+1)
%       mean: centroid xy coordinate of the polytope
%       area: area of the polytope
%       max_radius: distance from the mean to the farthest vertex
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_plotPolytopes
%     fcn_MapGen_polytopeShrinkFromEdges
%
% EXAMPLES:
%
% For additional examples, see: script_test_fcn_MapGen_polytopesShrinkFromEdges
%
% This function was written on 2022_01_17 by Steve Harentt using fcn_MapGen_polytopesShrinkToRadius
% as a template
% Questions or comments? sjh6473@psu.edu
%

% Revision History:
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% -- fixed call to fcn_MapGen_fillPolytopeFieldsFromVertices

% TO DO
% copied from TODOs in fcn_MapGen_polytopesShrinkToRadius
% -- Vectorize the for loop if possible
% -- check inputs are positive numbers where appropriate (e.g. make a
% "positive number" check
% -- add non uniform shrinking and allow specification of a gap size variance,
%    similar to how a radius variance is allowed in fcn_MapGen_polytopesShrinkToRadius.m

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
    if flag_check_inputs
        % Are there the right number of inputs?
        if nargin < 2 || nargin > 3
            error('Incorrect number of input arguments')
        end

        % Check the polytopes input
        fcn_DebugTools_checkInputsToFunctions(...
            polytopes, 'polytopes');

        % Check the des_radius input
        fcn_DebugTools_checkInputsToFunctions(...
            des_gap_size, 'positive_1column_of_numbers',1);
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
        fig_for_debug = fig.Number; 
        flag_do_plot = 1;
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

initial_stats = fcn_MapGen_polytopesStatistics(polytopes);
initial_average_max_rad = initial_stats.average_max_radius;

if flag_do_debug
    fprintf(1,'Target statistics:\n');
    fprintf(1,'\tGap size: %.4f\n',des_gap_size);

    fprintf(1,'Input field statistics:\n');
    fprintf(1,'\tAvg Max Rad: %.4f\n',initial_average_max_rad);

    % fcn_MapGen_OLD_plotPolytopes(polytopes,fig_for_debug,'b',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [0 0 1];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat), (fillFormat), (fig_for_debug)); %#ok<NASGU>

end

%% shrink polytopes to achieve the distribution
% Check that the old polytopes are large enough to shrink and
% achieve the new radius distribution. Want all the changes to be smaller
% than -2 times the minimum radius, to ensure we do not get singular
% polytopes.

% Initialize the shrunk polytopes structure array, and tolerance for
% distance between vertices, below which vertices are merged into one.
clear shrunkPolytopes
shrinker = polytopes(1); % obstacle to be shrunk
temp = ...
    fcn_MapGen_polytopeShrinkFromEdges(...
    shrinker,des_gap_size/2);
shrunkPolytopes(length(polytopes)) = temp;

% Loop through each polytope, shrinking it to the reference size
for ith_poly = 1:length(polytopes)
    shrinker = polytopes(ith_poly); % obstacle to be shrunk

    if isnan(shrinker.vertices(1,1)) % Degenerate
        shrunkPolytopes(ith_poly)=shrinker;
    else
        % assign to shrunkPolytopes
        % gap_size over 2 is the normal distance to pull edges in
        shrunkPolytopes(ith_poly) = ...
            fcn_MapGen_polytopeShrinkFromEdges(...
            shrinker,des_gap_size/2);
    end
end



if flag_do_debug
    final_stats = fcn_MapGen_polytopesStatistics(shrunkPolytopes);
    final_average_max_rad = final_stats.average_max_radius;
    fprintf(1,'Final distrubution statistics:\n');
    fprintf(1,'\tAvg mag rad: %.4f\n',final_average_max_rad);
end

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

if flag_do_plot
    figure(fig_num);
    hold on

    % Plot the input polytopes in red
    % fcn_MapGen_OLD_plotPolytopes(polytopes,fig_num,'r',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [1 0 0];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat), (fillFormat), (fig_num)); %#ok<NASGU>

    % plot the shrunk in blue
    fcn_MapGen_plotPolytopes(shrunkPolytopes,fig_num,'b',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [0 0 1];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(shrunkPolytopes, (plotFormat), (fillFormat), (fig_num)); %#ok<NASGU>


end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function
