function output_pts = fcn_MapGen_snapInteriorPointToVertex(polytopes, pts_to_test)
% fcn_MapGen_snapInteriorPointToVertex
% if a point is in a polytope, it places the point on the nearest vertex of polytope
%
%
%
% FORMAT:
%
%    [ ...
%    exp_polytopes ...
%    ] = ...
%    fcn_MapGen_polytopesExpandEvenlyForConcave( ...
%    polytopes, ...
%    delta, ...
%    exp_dist, ...
%    (fig_num) ...
%    )
%
% INPUTS:
%
%     polytopes: the structure of 'polytopes' type that stores the
%     polytopes to be expanded
%
%     exp_dist: distance to expand the obstacle
%
%     (optional inputs)
%
%     fig_num: any number that acts somewhat like a figure number output.
%     If given, this forces the variable types to be displayed as output
%     and as well makes the input check process verbose.
%
%
% OUTPUTS:
%
%     exp_polytopes: structure of expanded polytopes
%
%
% DEPENDENCIES:
%
%     fcn_MapGen_checkInputsToFunctions
%     fcn_MapGen_fillPolytopeFieldsFromVertices
%     fcn_MapGen_plotPolytopes
%     MATLAB's polyshape object and polybuffer object function (method)
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_polytopesExpandEvenlyForConcave
% for a full test suite.
%
% This function was written 5 Feb. 2024 by Steve Harnett
% Questions or comments? contact sjharnett@psu.edu

%
% REVISION HISTORY:
%
% 2024_02_05, Steve Harnett
% -- first write of script

%
% TO DO:
%
% -- fill in to-do items here.
    % TODO this needs to go in MapGen
    % TODO @sjharnett make function for snapping point to nearest vertex of polytope
    % enlarging polytopes may have put the midway start inside a polytope
    % for each polytope, check if this point is inside the polytope
    output_pts = pts_to_test;
    for p = 1:length(polytopes)
        these_verts = polytopes(p).vertices;
        this_polyshape = polyshape(these_verts);
        % is point in but not on polyshape?
        [is_in,is_on] = isinterior(this_polyshape,pts_to_test);
        pts_in = find(is_in);
        if isempty(pts_in)
            continue
        end
        for pt_in_idx = 1:length(pts_in)
            pt_to_test = pts_to_test(pts_in(pt_in_idx),:);
            % if it is, get distance to all vertices
            vert_to_start_deltas = these_verts - pt_to_test;
            vert_to_start_distances = vert_to_start_deltas(:,1).^2 + vert_to_start_deltas(:,2).^2;
            [min_value, idx_of_min] = min(vert_to_start_distances);
            % set this point to the nearest vertex
            output_pts(pts_in(pt_in_idx),:) = these_verts(idx_of_min,:);
        end
    end
end
