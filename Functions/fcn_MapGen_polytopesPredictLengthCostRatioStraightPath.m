function r_lc_straight_through = fcn_MapGen_polytopesPredictLengthCostRatioStraightPath(...
    pre_shrink_polytopes,polytopes,des_gap_size,start_x,start_y,finish_x,finish_y)
    % fcn_MapGen_polytopesPredictLengthCostRatioStraightPath
    % Given an polytope field, predict the length cost ratio of driving through
    % the entire field in a straight line
    %
    %
    %
    % FORMAT:
    % r_lc_straight_through = fcn_MapGen_polytopesPredictLengthCostRatioStraightPath(...
    %   pre_shrink_polytopes,polytopes,des_gap_size,start_x,start_y,finish_x,finish_y)
    %
    % INPUTS:
    %
    % pre_shrink_polytopes - the fully tiled polytope array prior to shrinking
    % polytopes - the polytopes with the desired shrinking applied
    % des_gap_size - the gap size applied to go from pre_shrink_polytopes to polytopes
    % start_x, start_y, finish_x, finish_y - the (x,y) coordinates of the start and finish points
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


    %% Debugging and Input checks
    flag_check_inputs = 1; % Set equal to 1 to check the input arguments
    flag_do_plot = 0;      % Set equal to 1 for plotting
    flag_do_debug = 0;     % Set equal to 1 for debugging
    if flag_do_debug
        st = dbstack; %#ok<*UNRCH>
        fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
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


    if 1 == flag_check_inputs

        % Are there the right number of inputs?
        if nargin < 7 || nargin > 7
            error('Incorrect number of input arguments')
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
    r_D = map_stats.avg_r_D;
    if r_D < 0.15 || r_D >= 0.53
        L_unocc = unocc_ests.L_unocc_est_poly_fit;
    else
        L_unocc = unocc_ests.L_unocc_est_slant_AABB_width;
    end
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
end
