function [ ...
convex_hull_overlap_ratio,...
A_overlap,...
A_occupied...
] = ...
fcn_MapGen_calculateConvexHullOverlapRatio( ...
polytopes, ...
varargin...
)
% fcn_MapGen_calculateConvexHullOverlapRatio
% calculates the convex hull of every obstacle.  The area of the overlap between
% these hulls relative to the total occupied area
%
% FORMAT:
%
% function [ ...
% convex_hull_overlap_ratio...
% ] = ...
% fcn_MapGen_calculateConvexHullOverlapRatio( ...
% polytopes, ...
% (fig_num)...
% )
%
% INPUTS:
%
%     polytopes: the structure of 'polytopes' type that stores the
%     polytopes to be expanded
%
%     (optional inputs)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
%
% OUTPUTS:
%
%     covex_hull_overlap_ratio: portion of overlapping convex hull area to total obstacle area
%
%
% DEPENDENCIES:
%
%     MATLAB's polyshape object and union object function (method)
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_calculateConvexHullOverlapRatio
% for a full test suite.
%
% This function was written 23 Feb. 2024 by Steve Harnett
% Questions or comments? contact sjharnett@psu.edu

%
% REVISION HISTORY:
%
% 2024_02_23, Steve Harnett
% -- first write of function
% 2025_04_16, Steve Harnett
% -- add legend to plotting
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_04_29 by Sean Brennan
% -- Fixed script mis-labeling in the docstrings

% TO DO
% -- rewrite to move plotting to debug area only

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==2 && isequal(varargin{end},-1))
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
        if nargin < 1 || nargin > 2
            error('Incorrect number of input arguments')
        end

        % Check the polytopes input, make sure it is 'polytopes' type
        fcn_DebugTools_checkInputsToFunctions(...
            polytopes, 'polytopes');


    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  2 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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

clear exp_polytopes;
A_overlap = 0;
A_occupied = 0;
conv_hull_polyshapes = [];
for p = 1:length(polytopes)
    this_polytope = polytopes(p); % look at one polytope
    if flag_do_plot
        figure(fig_num); hold on; box on; fill(this_polytope.vertices(:,1)',this_polytope.vertices(:,2),[0 0 1],'FaceAlpha',1);
        if p == 1
            leg_str = {'obstacles'};
        else
            leg_str{end+1} = ''; %#ok<AGROW>
        end
    end
    these_vertices = this_polytope.vertices(1:(end-1),:); % grab only non-repeating vertices
    k = convhull(these_vertices); % find convex hull of vertices
    convex_hull_vertices = [these_vertices(k,1),these_vertices(k,2)];
    this_conv_hull_polyshape = polyshape(convex_hull_vertices); % convert it to matlab polyshape
    if flag_do_plot
        figure(fig_num); hold on; box on; plot(this_conv_hull_polyshape,'FaceColor','green','FaceAlpha',0.2);
        if p == 1
            leg_str{end+1} = 'obs. convex hull'; %#ok<AGROW>
        else
            leg_str{end+1} = ''; %#ok<AGROW>
        end
    end
    conv_hull_polyshapes = [conv_hull_polyshapes; this_conv_hull_polyshape]; %#ok<AGROW>
    A_occupied = A_occupied+this_polytope.area;
end
flag_havent_plotted = 1;
for i = 1:length(polytopes)
    j = 1;
    while j <= length(polytopes)
        if j<=i
            % if we checked 1,2 we don't need to check 2,1 so ignore j<i
            % also don't want to check j=i because there is meaningless overlap
            j = j+1;
            continue;
        end
        overlap_polyshape = intersect(conv_hull_polyshapes(i),conv_hull_polyshapes(j));
        if flag_do_plot
            figure(fig_num); hold on; box on; plot(overlap_polyshape,'FaceColor','red');
            if flag_havent_plotted
                leg_str{end+1} = 'overlap'; %#ok<AGROW>
                flag_havent_plotted = 0;
            end
        end
        A_overlap = A_overlap + area(overlap_polyshape);
        j = j+1;
    end
end

convex_hull_overlap_ratio = A_overlap/A_occupied;

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
    title_str = sprintf('total obstacle area: %.3f\noverlapping convex hull area: %.3f\nconvex hull overlap ratio: %.3f',A_occupied,A_overlap,convex_hull_overlap_ratio);
    figure(fig_num); title(title_str);
    legend(leg_str)
end

end % Ends the function
