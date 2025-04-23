function output_pts = fcn_MapGen_snapInteriorPointToVertex(polytopes, pts_to_test)
% fcn_MapGen_snapInteriorPointToVertex
% if a point is in a polytope, it places the point on the nearest vertex of polytope
%
%
%
% FORMAT:
%
% output_pts = fcn_MapGen_snapInteriorPointToVertex(polytopes, pts_to_test)
%
% INPUTS:
%
%     polytopes: the structure of 'polytopes' type that stores the
%     polytopes to be expanded
%
%     pts_to_test: Nx2 matrix of (X,Y) points for N points to snap to vertices, if interior to
%           a polytope in 'polytopes'
%
%
%
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

%
% TO DO:
%
% -- fill in to-do items here.
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
