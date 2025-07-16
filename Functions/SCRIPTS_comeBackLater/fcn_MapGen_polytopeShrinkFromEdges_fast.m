function [shrunkPolytope, new_vertices, new_projection_vectors, cut_distance] = ...
    fcn_MapGen_polytopeShrinkFromEdges_fast(...
    shrinker,...
    edge_cut,...
    varargin)
% fcn_MapGen_polytopeShrinkFromEdges_fast cuts edges off the polytopes
% Each edge is cut so that the entire polytope is trimmed exactly the same
% amount from each edge. This fast variant does NOT return the fully
% populated polytope, nor does it check inputs. It simply calculates the
% verticies field without updating areas, etc.
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

% [shrunkPolytope, (new_vertices, new_projection_vectors, cut_distance)]= ...
%     fcn_MapGen_polytopeShrinkFromEdges(...
%     shrinker,...
%     edge_cut,...
%     (fig_num))
%
% OR:
%
% [shrunkPolytope, (new_vertices, new_projection_vectors, cut_distance)]= ...
%     fcn_MapGen_polytopeShrinkFromEdges(...
%     shrinker,...
%     edge_cut,...
%     (new_vertices, new_projection_vectors, cut_distance),...
%     (fig_num))
%
% INPUTS:
%
%     shrinker: original polytope with same fields as shrunkPolytopes
%     below
%
%     edge_cut: desired cut distance from each edge
%
%    (OPTIONAL INPUTS)
%
%    [new_vertices, new_projection_vectors, cut_distance] : outputs from
%    the function: fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num)
%    or outputs from previous calls, used to speed up code since this
%    skeleton calculation is by far the slowest part.
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
%    [new_vertices, new_projection_vectors, cut_distance] : outputs from
%    the function: fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num)
%    or outputs from previous calls, used to speed up code since this
%    skeleton calculation is by far the slowest part.
%
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_polytopeFindVertexAngles
%     fcn_MapGen_fillPolytopeFieldsFromVertices
%
% % EXAMPLES:
%
%
% For additional examples, see: script_test_fcn_MapGen_polytopeShrinkFromEdges
%
% This function was written on 2021_08_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%

% Revision History:
% 2021_08_02 - S. Brennan
% -- first write of code
% 2022_02_13 - S.Brennan
% -- supress MATLAB's warning about flags
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions


% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
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
        narginchk(2,6);

        % Check the shrinker input
        fcn_DebugTools_checkInputsToFunctions(...
            shrinker, 'polytopes');

        % Check the edge_cut input
        fcn_DebugTools_checkInputsToFunctions(...
            edge_cut, 'positive_1column_of_numbers',1);

    end
end

% Does user want to input skeleton values?
flag_use_user_skeleton = 0;
if  3 < nargin % Only way this happens is if user specifies skeleton
    if nargin<5
        error('Incorrect number of input arguments');
    end

    new_vertices = varargin{1};
    new_projection_vectors = varargin{2};
    cut_distance = varargin{3};

    % if flag_check_inputs
    %     % Check the cut_distance input
    %     fcn_DebugTools_checkInputsToFunctions(...
    %         new_vertices, '1column_of_numbers');
    % end

    flag_use_user_skeleton = 1;

end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  ((6 == nargin) || (3==nargin)) && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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

% Initialize variables we may need
vertices = shrinker.vertices;

%% STEP 1. Calculate the polytope skeleton,
% or use the skeleton if the user provides this as inputs.

% Do we need to calculate skeleton values?
if 0 == flag_use_user_skeleton
    if flag_do_plot
        [new_vertices, new_projection_vectors, cut_distance] = ...
            fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num);
    else
        [new_vertices, new_projection_vectors, cut_distance] = ...
            fcn_MapGen_polytopeFindVertexSkeleton(vertices);
    end
end

%% STEP 2. Using the cut distance, find the template
% We want to use the one that is less than or equal to the cut distance. If
% less, then calculate the additional cut and project the verticies to
% their new points based on residual cut.


% Find the shape that is less than or equal to the cut
shape_index = find(cut_distance<=edge_cut,1,'last');

% Grab vertices to start from, cut to start from
template_vertices = new_vertices{shape_index};
template_start_cut = cut_distance(shape_index);

% Calculate projection distance
additional_cut_distance = edge_cut - template_start_cut;

% Determine final vertices
final_vertices = template_vertices + new_projection_vectors{shape_index}*additional_cut_distance;



%% STEP 3. Convert the resulting verticies into the standard polytope form
% Fill in the results
shrunkPolytope.vertices = final_vertices;

% SKIPPING FOR SPEED: fill in other fields from the vertices field
% shrunkPolytope = fcn_MapGen_fillPolytopeFieldsFromVertices(shrunkPolytope);



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
    grid minor
    hold on
    axis equal

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

