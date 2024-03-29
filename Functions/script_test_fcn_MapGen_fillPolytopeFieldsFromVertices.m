% script_test_fcn_MapGen_fillPolytopeFieldsFromVerticies
% Tests: fcn_MapGen_fillPolytopeFieldsFromVerticies

%
% REVISION HISTORY:
%
% 2021_07_02 by Sean Brennan
% -- first write of script
% 2023_01_15 by Sean Brennan
% -- added figure number
%%%%%%%%%%%%%%§



fig_num = 2;
line_width = 3;
clear polytopes
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes);
fcn_MapGen_plotPolytopes(polytopes,fig_num,'r-',line_width);

assert(isequal(polytopes.vertices,[0,0;4,2;2,4;0,0]));
assert(isequal(polytopes.xv,[0,4,2]));
assert(isequal(polytopes.yv,[0,2,4]));
assert(isequal(round(polytopes.distances,4),[4.4721;2.8284;4.4721]));
assert(isequal(polytopes.mean,[2,2]));
assert(isequal(polytopes.area,6));
assert(isequal(round(polytopes.max_radius,4),2.8284));

polytopes(2).vertices = [10 10; 14 21; 12 41; 10 10];
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes);
fcn_MapGen_plotPolytopes(polytopes,fig_num,'r-',line_width);
