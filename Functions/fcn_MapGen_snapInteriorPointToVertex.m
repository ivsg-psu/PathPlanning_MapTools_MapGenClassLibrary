function output_pts = fcn_MapGen_snapInteriorPointToVertex(polytopes, pts_to_test, varargin)
% fcn_MapGen_snapInteriorPointToVertex
% if a point is in a polytope, it places the point on the nearest vertex of polytope
%
%
%
% FORMAT:
%
% output_pts = fcn_MapGen_snapInteriorPointToVertex(polytopes, pts_to_test, (fig_num))
%
% INPUTS:
%
%     polytopes: the structure of 'polytopes' type that stores the
%     polytopes to be expanded
%
%     pts_to_test: Nx2 matrix of (X,Y) points for N points to snap to vertices, if interior to
%           a polytope in 'polytopes'
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
%     output_pts: Nx2 matrix of N (X,Y) positions.  For each of the N points in pts_to_test,
%       the row in 'output_pts' is either the nearest polytope verticex if the input point was internal
%       to a polytope in 'polytopes' or the output point is the same as the input point if the input
%       is not within any polytope
%
%
% DEPENDENCIES:
%
%     MATLAB's polyshape object and isInterior polyshape function (method)
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_snapInteriorPointToVertex
% for a full test suite.
%
% This function was written 5 Feb. 2024 by Steve Harnett
% Questions or comments? contact sjharnett@psu.edu

%
% REVISION HISTORY:
%
% 2024_02_05, Steve Harnett
% -- first write of script
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

        % % Check the polytopes input, make sure it is 'polytopes' type
        % fcn_DebugTools_checkInputsToFunctions(...
        %     polytopes, 'polytopes');
        % 
        % 
        % % Check the exp_dist input, make sure it is 'positive_column_of_numbers' type
        % fcn_DebugTools_checkInputsToFunctions(...
        %     exp_dist, 'positive_1column_of_numbers',1);

    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  (3 == nargin) && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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


output_pts = pts_to_test;
for p = 1:length(polytopes)
    these_verts = polytopes(p).vertices;
    this_polyshape = polyshape(these_verts);
    % is point in but not on polyshape?
    [is_in,is_on] = isinterior(this_polyshape,pts_to_test); %#ok<ASGLU>
    pts_in = find(is_in);
    if isempty(pts_in)
        continue
    end
    for pt_in_idx = 1:length(pts_in)
        pt_to_test = pts_to_test(pts_in(pt_in_idx),:);
        % if it is, get distance to all vertices
        vert_to_start_deltas = these_verts - pt_to_test;
        vert_to_start_distances = vert_to_start_deltas(:,1).^2 + vert_to_start_deltas(:,2).^2;
        [min_value, idx_of_min] = min(vert_to_start_distances); %#ok<ASGLU>
        % set this point to the nearest vertex
        output_pts(pts_in(pt_in_idx),:) = these_verts(idx_of_min,:);
    end
end
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

    LineWidth = 2;
    fcn_MapGen_plotPolytopes(polytopes,fig_num,'r-',LineWidth);
    fcn_MapGen_plotPolytopes(exp_polytopes,fig_num,'b-',LineWidth,'square');
    legend('Original','Expanded')
    box on
    xlabel('X Position')
    ylabel('Y Position')

end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends main function

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


