function r_lc_straight_through = fcn_MapGen_polytopesPredictLengthCostRatioStraightPath(...
    pre_shrink_polytopes,polytopes,des_gap_size,start_x,start_y,finish_x,finish_y, varargin)
% fcn_MapGen_polytopesPredictLengthCostRatioStraightPath
% Given an polytope field, predict the length cost ratio of driving through
% the entire field in a straight line
%
%
%
% FORMAT:
% r_lc_straight_through = fcn_MapGen_polytopesPredictLengthCostRatioStraightPath(...
%   pre_shrink_polytopes,polytopes,des_gap_size,start_x,start_y,finish_x,finish_y, (fig_num))
%
% INPUTS:
%
% pre_shrink_polytopes - the fully tiled polytope array prior to shrinking
% polytopes - the polytopes with the desired shrinking applied
% des_gap_size - the gap size applied to go from pre_shrink_polytopes to polytopes
% start_x, start_y, finish_x, finish_y - the (x,y) coordinates of the start and finish points
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
%     r_lc_straight_through - the estimated length cost ratio of driving straight through the field
%
% DEPENDENCIES:
%
%     fcn_MapGen_polytopesPredictUnoccupancyRatio
%
% EXAMPLES:
%
% See the script: script_planning_performed_at_multiple_costs.m
% in the repo PathPlanning_GridFreePathPlanners_BoundedAStar
% for a comparison of predicted costs from this function to costs from a straight path planner
%
% Questions or comments? contact sjh6473@psu.edu

% REVISION HISTORY:
% 2022_05_20
% -- first written by Steve Harnett
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_07_17 by Sean Brennan
% -- standardized Debugging and Input checks area, Inputs area
% -- made codes use MAX_NARGIN definition at top of code, narginchk
% -- made plotting flag_do_plots and code consistent across all functions

% TO DO
% -- add better input checking

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 8; % The largest Number of argument inputs to the function
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
        narginchk(7,MAX_NARGIN);

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

% find area unoccupancy ratio
unocc_ests = fcn_MapGen_polytopesPredictUnoccupancyRatio(pre_shrink_polytopes,polytopes,des_gap_size);
map_stats = fcn_MapGen_polytopesStatistics(polytopes);
r_D = map_stats.avg_r_D; %#ok<NASGU>
L_unocc = unocc_ests.L_unocc_est_avg_circle_min_rad_est_1;
% convert to occupancy ratio
L_occ = 1-L_unocc;
% find average cost for polytope field
cost_avg = mean(extractfield(polytopes,'cost'));
% find distance from start to finish
a_b = ((finish_x - start_x)^2 + (finish_y - start_y)^2)^0.5;
% scale percent of occupied distance from start to finish by 1+cost
cost_in_polys = a_b * L_occ *(1+cost_avg);
% add on distance outside of polytopes as unscaled length cost
total_cost = cost_in_polys + (1-L_occ)*a_b;
% find cost ratio as cost/distance travelled
r_lc_straight_through = total_cost/a_b;

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
    figure(fig_num)
    clf;


end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends the function

