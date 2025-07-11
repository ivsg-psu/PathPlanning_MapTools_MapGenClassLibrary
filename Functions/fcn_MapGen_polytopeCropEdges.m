function [trim_polytopes] = ...
    fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,varargin)
% fcn_MapGen_polytopeCropEdges removes polytopes that extend
% beyond the boundaries specified
%
% FORMAT:
%
% [trim_polytopes] = ...
%     fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,(fig_num))
%
% INPUTS:
%
%     polytopes: the original polytopes with the same fields as trim_polytopes
%
%     bounding_box: a 2 x 2 matrix of [xlow ylow; xhigh yhigh] in which all
%     the polytopes must exist, e.g. the corner coordinates of the
%     axis-aligned bounding box.
%
%    (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
%
% OUTPUTS:
%
%     TRIM_POLYTOPES: a 1-by-n seven field structure of polytopes within the
%     boundaries, where n <= number of polytopes with fields:
%     vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%     the number of the individual polytope vertices
%     xv: a 1-by-m vector of vertice x-coordinates
%     yv: a 1-by-m vector of vertice y-coordinates
%     distances: a 1-by-m vector of perimeter distances from one point to the
%     next point, distances(i) = distance from vertices(i) to vertices(i+1)
%     mean: centroid xy coordinate of the polytope
%     area: area of the polytope
%
%
% EXAMPLES:
%
%   polytopes = fcn_MapGen_haltonVoronoiTiling([1 1000]);
%   fig = fcn_plot_polytopes(polytopes,[],'b',2,[0 1 0 1]);
%
%   bounding_box = [0,1; 0,1];
%   trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box);
%   fcn_plot_polytopes(trim_polytopes,fig,'g',2,[0 1 0 1]);
%
% For additional examples, see: script_test_fcn_MapGen_polytopeCropEdges
%
% This function was written on 2019_06_13 by Seth Tau
% Questions or comments? sat5340@psu.edu
%

% REVISIONS
% 2021-06-06
% -- rewrote function from fc_polytope_editing_remove_edge_polytopes
% -- added a test script
% 2021-06-10
% -- updated comments
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions


% TO DO
% -- Vectorize the for loop if possible

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

        % Check the bounding_box input
        fcn_DebugTools_checkInputsToFunctions(...
            bounding_box, '2column_of_numbers',2);

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlow = bounding_box(1,1);
ylow = bounding_box(1,2);
xhigh = bounding_box(2,1);
yhigh = bounding_box(2,2);

Npolys = length(polytopes);

% Preallocate the polytopes
trim_polytopes(Npolys) = fcn_MapGen_polytopeFillEmptyPoly((-1));

keep = 0; % number of polytopes to keep
for poly = 1:Npolys % check each polytope within polytopes
    xv = polytopes(poly).xv;
    yv = polytopes(poly).yv;
    if sum((xv<xlow)+(xv>xhigh)+(yv<ylow)+(yv>yhigh))==0 % if the x or y vertices are inside of the bounds
        keep = keep + 1;
        trim_polytopes(keep) = polytopes(poly);
    end
end
trim_polytopes = trim_polytopes(1:keep); % Save only the ones that were filled

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
    fcn_MapGen_plotPolytopes(polytopes,fig_num,'r',2,[xlow xhigh ylow yhigh]);

    % plot the tiled_polytopes in blue
    fcn_MapGen_plotPolytopes(trim_polytopes,fig_num,'b',2,[xlow xhigh ylow yhigh]);

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends the function



