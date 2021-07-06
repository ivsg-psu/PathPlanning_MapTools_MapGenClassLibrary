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

