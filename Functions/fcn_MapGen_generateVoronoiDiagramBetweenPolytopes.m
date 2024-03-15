function [vx,vy,h] = fcn_MapGen_generateVoronoiDiagramBetweenPolytopes(polytopes,is_nonconvex)
    % fcn_MapGen_increasePolytopeVertexCount
    % Given polytope field and a desired resolution distance, n, returns an equivalent
    % polytope field with colinear vertices added to each polytope side such that
    % there is a vertex every n units
    % The utility of this is that if path planning is restricted to using polytope
    % vertices as waypoints, this increases the number of options the planner has
    % while keeping the obstacle field the same
    %
    %
    %
    % FORMAT:
    % interpolated_polytopes = fcn_MapGen_increasePolytopeVertexCount(polytopes,resolution)
    %
    % INPUTS:
    %     polytopes - the initial polytope field
    %     resolution - the desired linear spacing between vertices along each polytope side
    %
    % OUTPUTS:
    %
    %
    %     interpolated_polytopes - a polytope field equivalent to the input but with vertices added
    %       along the polytopes sides every RESOLUTION units such that each polytope now has more vertices
    %
    % DEPENDENCIES:
    %
    % EXAMPLES:
    %
    % See the script: script_fcn_MapGen_increasePolytopeVertexCount.m
    % for a full test suite.
    %
    % Questions or comments? contact sjh6473@psu.edu

    % REVISION HISTORY:
    % 2021_10_13
    % -- first written by Steve Harnett
    if ~is_nonconvex
        [vx,vy] = voronoi([polytopes.xv],[polytopes.yv]);
        h = voronoi([polytopes.xv],[polytopes.yv]);
    else
        % TODO @sjharnett min of diff between all points /2 so at least every side is cut in half
        distances = diff([[polytopes.xv]',[polytopes.yv]']);
        min_distance_between_verts = min(sqrt(sum(distances.*distances,2)));
        % poly_map_stats = fcn_MapGen_polytopesStatistics(polytopes);
        % % want to ensure that a side with length of 2 std dev below mean is still interpolated at least in half
        % resolution = (poly_map_stats.average_side_length - 2*poly_map_stats.std_side_length)/2;
        resolution = min_distance_between_verts/2;
        interpolated_polytopes = fcn_MapGen_increasePolytopeVertexCount(polytopes, resolution);
        [vx,vy] = voronoi([interpolated_polytopes.xv], [interpolated_polytopes.yv]);
        h = voronoi([interpolated_polytopes.xv], [interpolated_polytopes.yv]);
    end
end
