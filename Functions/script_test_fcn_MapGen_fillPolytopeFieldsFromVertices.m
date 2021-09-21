% script_test_fcn_MapGen_fillPolytopeFieldsFromVerticies
% Tests: fcn_MapGen_fillPolytopeFieldsFromVerticies

%
% REVISION HISTORY:
%
% 2021_07_02 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§




clear polytopes
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes);

assert(isequal(polytopes.vertices,[0,0;4,2;2,4;0,0]));
assert(isequal(polytopes.xv,[0,4,2]));
assert(isequal(polytopes.yv,[0,2,4]));
assert(isequal(round(polytopes.distances,4),[4.4721;2.8284;4.4721]));
assert(isequal(polytopes.mean,[2,2]));
assert(isequal(polytopes.area,6));
assert(isequal(round(polytopes.max_radius,4),2.8284));
