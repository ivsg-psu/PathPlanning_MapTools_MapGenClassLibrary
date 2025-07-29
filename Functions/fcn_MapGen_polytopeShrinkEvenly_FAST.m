function [shrunkPolytope, newVertices, newProjectionVectors, cutDistance] = ...
    fcn_MapGen_polytopeShrinkEvenly_FAST(...
    unshrunkPolytope,...
    edgeCut,...
    varargin)
% fcn_MapGen_polytopeShrinkEvenly_FAST cuts edges off the polytopes
% Each edge is cut so that the entire polytope is trimmed exactly the same
% amount from each edge. Uses the vertex skeleton method.
%
% This is implemented in three steps:
% 1. Calculate the polytope skeleton, or use the skeleton if the user
% provides this as inputs.
%
% 2. Using the cut distance, find the template that is less than or equal
% to the cut distance. If less, then calculate the additional cut and
% project the verticies to their new points based on residual cut.
%
% 3. Convert the resulting verticies into the standard polytope form.
%
% FORMAT:
%
% [shrunkPolytope, (newVertices, newProjectionVectors, cutDistance)]= ...
%     fcn_MapGen_polytopeShrinkEvenly_FAST(...
%     unshrunkPolytope,...
%     edgeCut,...
%     (fig_num))
%
% OR:
%
% [shrunkPolytope, (newVertices, newProjectionVectors, cutDistance)]= ...
%     fcn_MapGen_polytopeShrinkEvenly_FAST(...
%     unshrunkPolytope,...
%     edgeCut,...
%     (precalcVertices, precalcProjectionVectors, precalcCutDistance),...
%     (fig_num))
%
% INPUTS:
%
%     unshrunkPolytope: original polytope with same fields as shrunkPolytopes
%     below
%
%     edgeCut: desired cut distance from each edge
%
%    (OPTIONAL INPUTS)
%
%    [precalcVertices, precalcProjectionVectors, precalcCutDistance] :
%    outputs from the function:
%    fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num) 
%    or outputs from previous calls, used to speed up code since this
%    skeleton calculation is by far the slowest code and only needs to be
%    calculated once per polytope.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
% OUTPUTS:
%
%     shrunkPolytopeS: a 1-by-n seven field structure of shrunken polytopes,
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
%    [newVertices, newProjectionVectors, cutDistance] : outputs from
%    the function: fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num)
%    or outputs from previous calls, used to speed up code since this
%    skeleton calculation is by far the slowest part.
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_polytopeFindVertexAngles
%     fcn_MapGen_fillPolytopeFieldsFromVertices
%
% % EXAMPLES:
%
% For additional examples, see: script_test_fcn_MapGen_polytopeShrinkEvenly
%
% This function was written on 2021_08_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%

% Revision History:
% 2021_08_02 - S. Brennan
% -- first write of code
% 2022_02_13 - S.Brennan
% -- supress MATLAB's warning about flags
% 2023_01_15 - S.Brennan
% -- clean up comments
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_07_29 by Sean Brennan
% -- added test script

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 6; % The largest Number of argument inputs to the function
flag_max_speed = 0;
if (nargin==MAX_NARGIN && isequal(varargin{end},-1)) || (nargin==3 && isequal(varargin{end},-1))
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

        % Check the unshrunkPolytope input
        fcn_DebugTools_checkInputsToFunctions(...
            unshrunkPolytope, 'polytopes');

        % Check the edgeCut input
        fcn_DebugTools_checkInputsToFunctions(...
            edgeCut, 'positive_1column_of_numbers',1);

    end
end

% Does user want to input skeleton values?
flag_use_user_skeleton = 0;
if  3 < nargin % Only way this happens is if user specifies skeleton
    if nargin<5
        error('Incorrect number of input arguments');
    end
    temp = varargin{1};
    if ~isempty(temp)
        newVertices = varargin{1};
        newProjectionVectors = varargin{2};
        cutDistance = varargin{3};

        % Check the cutDistance input
        % fcn_DebugTools_checkInputsToFunctions(...
        %     newVertices, '1column_of_numbers');

        flag_use_user_skeleton = 1;
    end

end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  (0==flag_max_speed) && ((MAX_NARGIN == nargin) || (3==nargin)) % Only create a figure if NOT maximizing speed
    temp = varargin{end}; % Last argument is always figure number
    if ~isempty(temp) % Make sure the user is not giving empty input
        fig_num = temp;
        flag_do_plot = 1; % Set flag to do plotting
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

% Initialize variables we may need
vertices = unshrunkPolytope.vertices;

%% STEP 1. Calculate the polytope skeleton,
% or use the skeleton if the user provides this as inputs.

% Do we need to calculate skeleton values?
if 0 == flag_use_user_skeleton
    [newVertices, newProjectionVectors, cutDistance] = ...
        fcn_MapGen_polytopeFindVertexSkeleton(vertices, -1);
end

%% STEP 2. Using the cut distance, find the template
% We want to use the one that is less than or equal to the cut distance. If
% less, then calculate the additional cut and project the verticies to
% their new points based on residual cut.


% Find the shape that is less than or equal to the cut
shape_index = find(cutDistance<=edgeCut,1,'last');

% Grab vertices to start from, cut to start from
template_vertices = newVertices{shape_index};
template_start_cut = cutDistance(shape_index);

% Calculate projection distance
additional_cutDistance = edgeCut - template_start_cut;

% Determine final vertices
final_vertices = template_vertices + newProjectionVectors{shape_index}*additional_cutDistance;

%% STEP 3. Convert the resulting verticies into the standard polytope form
% Fill in the results
shrunkPolytope.vertices = final_vertices;

% SKIPPING FOR SPEED: fill in other fields from the vertices field
% shrunkPolytope = fcn_MapGen_polytopesFillFieldsFromVertices(shrunkPolytope, [], -1);


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
    grid on
    hold on
    axis equal

    % Plot the input unshrunkPolytope in red
    % fcn_MapGen_plotPolytopes(unshrunkPolytope,fig_num,'r',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [1 0 0];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(unshrunkPolytope, (plotFormat), (fillFormat), (fig_num)); 
    set(h_plot,'DisplayName','Input: unshrunkPolytope');

    % plot the output polytope in blue
    % fcn_MapGen_OLD_plotPolytopes(shrunkPolytope,fig_num,'b',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [0 0 1];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(shrunkPolytope, (plotFormat), (fillFormat), (fig_num)); 
    set(h_plot,'DisplayName','Output: shrunkPolytope');

    legend('Interpreter','none','Location','best');

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function

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

% (none)

