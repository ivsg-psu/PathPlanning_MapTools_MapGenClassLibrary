function [convexHullOverlapRatio, areaOverlap, areaOccupied] = ...
    fcn_MapGen_calculateConvexHullOverlapRatio(polytopes, varargin)

% fcn_MapGen_calculateConvexHullOverlapRatio
% calculates the convex hull of every obstacle.  The area of the overlap between
% these hulls relative to the total occupied area
%
% FORMAT:
%
%     [convexHullOverlapRatio, areaOverlap, areaOccupied] = ...
%     fcn_MapGen_calculateConvexHullOverlapRatio( polytopes, (fig_num))
%
% INPUTS:
%
%     polytopes: the structure of 'polytopes' type that stores the
%     polytopes to be expanded. 
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
%     covex_hull_overlap_ratio: portion of overlapping convex hull area to total obstacle area
%
%     areaOverlap: area of overlap
%
%     areaOccupied: area occupied
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     MATLAB's polyshape object and union object function (method)
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_calculateConvexHullOverlapRatio
% for a full test suite.
%
% This function was written 2024_02_23 by Steve Harnett
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
% 2025_07_15 by Sean Brennan
% -- Fixed missing output descriptions in docstrings
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

        % Check the polytopes input, make sure it is 'polytopes' type
        fcn_DebugTools_checkInputsToFunctions(polytopes, 'polytopes');
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

Npolytopes = length(polytopes);

% Initialize variables
areaOverlap = 0;
areaOccupied = 0;
conv_hull_polyshapes(Npolytopes) = polyshape; 

% Loop through the polytopes, adding up area occupied
for ith_polytope = 1:Npolytopes
    this_polytope = polytopes(ith_polytope); % look at one polytope

    these_vertices = this_polytope.vertices(1:(end-1),:); % grab only non-repeating vertices

    % find convex hull of vertices
    k = convhull(these_vertices); 
    convex_hull_vertices = [these_vertices(k,1),these_vertices(k,2)];


    this_conv_hull_polyshape = polyshape(convex_hull_vertices); % convert it to matlab polyshape
    conv_hull_polyshapes(ith_polytope) = this_conv_hull_polyshape; 
    areaOccupied = areaOccupied + this_polytope.area;
end

% Look at overlap by intersecting the convex hulls with each other
NumOverlaps = ((Npolytopes-1)+1)/2 * (Npolytopes-1);
overlap_polyshapes(NumOverlaps) = polyshape;
Noverlaps = 0;
for ith_polytope = 1:Npolytopes-1
    % if we checked 1,2 we don't need to check 2,1 so ignore j<i
    % also don't want to check j=i because there is meaningless overlap
    for jth_overlap = (ith_polytope+1):Npolytopes
        overlap_polyshape = intersect(conv_hull_polyshapes(ith_polytope),conv_hull_polyshapes(jth_overlap));  

        Noverlaps = Noverlaps + 1;
        overlap_polyshapes(Noverlaps) = overlap_polyshape; 

        areaOverlap = areaOverlap + area(overlap_polyshape);
    end
end

% Calculate overlap ratio
convexHullOverlapRatio = areaOverlap/areaOccupied;

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
if flag_do_plots

    % Prep the figure
    figure(fig_num); 
    hold on; 
    box on;
 
    % Loop through the polytopes, shading them
    for ith_polytope = 1:Npolytopes
        this_polytope = polytopes(ith_polytope); % look at one polytope
        this_conv_hull_polyshape = conv_hull_polyshapes(ith_polytope);

        h_fill = fill(this_polytope.vertices(:,1)',this_polytope.vertices(:,2),[0 0 1],'FaceAlpha',1);
        if ith_polytope == 1
            set(h_fill,'DisplayName',sprintf('Obstacles, Area: %.3f',areaOccupied));
        else
            set(h_fill,'HandleVisibility','off');
        end

        h_plot = plot(this_conv_hull_polyshape,'FaceColor','green','FaceAlpha',0.2);
        if ith_polytope == 1
            set(h_plot,'DisplayName','Convex hull of obstacles'); 
        else
            set(h_plot,'HandleVisibility','off');
        end

    end

    for ith_overlap = 1:length(overlap_polyshapes)
        h_overlap = plot(overlap_polyshape(ith_overlap),'FaceColor','red');
        if 1==ith_overlap
            set(h_overlap,'DisplayName',sprintf('Overlap, Area: %.3f, ratio: %.3f',areaOverlap,convexHullOverlapRatio));
        else
            set(h_overlap,'HandleVisibility','off');
        end
    end

    legend('Interpreter','none');

end

end % Ends the function
