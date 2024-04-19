function output_pts = fcn_MapGen_snapInteriorPointToVertex(polytopes, pts_to_test)
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
        for pt_in = 1:length(pts_in)
            pt_to_test = pts_to_test(pts_in(pt_in));
            % if it is, get distance to all vertices
            vert_to_start_deltas = these_verts - pt_to_test;
            vert_to_start_distances = vert_to_start_deltas(:,1).^2 + vert_to_start_deltas(:,2).^2;
            [min_value, idx_of_min] = min(vert_to_start_distances);
            % set this point to the nearest vertex
            output_pts(pts_in(pt_in),:) = these_verts(idx_of_min,:);
        end
    end
end
