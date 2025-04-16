function [vx,vy,h] = fcn_MapGen_generateVoronoiDiagramBetweenPolytopes(polytopes,is_nonconvex)
    % fcn_MapGen_increasePolytopeVertexCount
    % Wraps the matlab voronoi() function to find the voronoi diagram using
    % the vertices of polytopes as the seed points, in the case of convex obstacles
    % or using densely packed colinear vertices along polytope sides as seed points
    % in the case of nonconvex osbtacles.  Voronoi diagram boundaries may collide
    % with polytope sides but voronoi diagram boundaries bewteen obstacles are also
    % calcualted.  This is similar to the process used by Masehian et al. 2004, see:
    % Masehian, Ellips, and M. R. Amin‐Naseri. "A voronoi diagram‐visibility graph‐potential field compound algorithm for robot path planning." Journal of Robotic Systems 21.6 (2004): 275-300.
    % To generate the medial axis between polytopes without these errant edges, you may
    % wish to use the fcn_MedialAxis_* functions in PathPlanning_GridFreePathPlanners_BoundedAStar
    %
    %
    % FORMAT:
    % [vx,vy,h] = fcn_MapGen_generateVoronoiDiagramBetweenPolytopes(polytopes,is_nonconvex)
    %
    % INPUTS:
    %     polytopes - the initial polytope field
    %     is_nonconvex - boolean flag indicating if there are or are not non-convex polytopes
    %
    % OUTPUTS:
    %
    %     Outputs are forwarded directly from Matlab's voronoi function.  See here for description:
    %       https://www.mathworks.com/help/matlab/ref/voronoi.html
    %
    % DEPENDENCIES:
    %     Matlab's voronoi function
    %     fcn_MapGen_increasePolytopeVertexCount
    %
    % EXAMPLES:
    %
    % See the script: script_fcn_MapGen_generateVoronoiDiagramBetweenPolytopes.m
    % for a full test suite.
    %
    % Questions or comments? contact sjh6473@psu.edu

    % REVISION HISTORY:
    % 2024_03_15
    % -- first written by Steve Harnett
    % 2025_04_16
    % -- commented by Steve Harnett
    if ~is_nonconvex
        [vx,vy] = voronoi([polytopes.xv],[polytopes.yv]);
        h = voronoi([polytopes.xv],[polytopes.yv]);
    else
        % min of diff between all points /2 so at least every side is cut in half
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
