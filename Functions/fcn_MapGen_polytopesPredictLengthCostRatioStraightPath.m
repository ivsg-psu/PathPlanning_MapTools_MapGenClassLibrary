function r_lc_straight_through = fcn_MapGen_polytopesPredictLengthCostRatioStraightPath(polytopes,des_gap_size,start,finish)
    % find area unoccupancy ratio
    unocc_ests = fcn_MapGen_polytopesPredictUnoccupancyRatio(polytopes,des_gap_size);
    A_unocc = unocc_ests.A_unocc_est_parallelogram;
    % convert to occupancy ratio
    A_occ = 1-A_unocc;
    % convert from area to linear
    L_occ = sqrt(A_occ);
    % find average cost for polytope field
    cost_avg = mean(polytopes.cost);
    % find distance from start to finish
    a_b = ((finish(1) - start(1))^2 + (finish(2) - start(2))^2)^0.5;
    % scale percent of occupied distance from start to finish by 1+cost
    cost_in_polys = a_b * L_occ *(1+cost_avg);
    % add on distance outside of polytopes as unscaled length cost
    total_cost = cost_in_polys + (1-L_occ)*a_b;
    % find cost ratio as cost/distance travelled
    r_lc_straight_through = total_cost/a_b;
end
