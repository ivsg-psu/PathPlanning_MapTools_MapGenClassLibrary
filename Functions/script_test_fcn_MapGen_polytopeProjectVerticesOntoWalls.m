% script_test_fcn_MapGen_polytopeFindSelfIntersections
% Tests function: fcn_MapGen_polytopeFindSelfIntersections

% REVISION HISTORY:
% 2021_08_03
% -- first written by S. Brennan

close all;


%% Basic example of self-intersection
fig_num = 1;
vertices = [0 0; 1 0; 0.5 1.5; 1 1; 0 1; 0 0];
vertices_with_self_intersects = fcn_MapGen_polytopeFindSelfIntersections(...
    vertices,fig_num);

fig_num = 11;
interior_point = [0.5 0.5];
[projected_points] = ...
    fcn_MapGen_polytopeProjectVerticesOntoWalls(...,
    interior_point,...
    vertices_with_self_intersects,...
    vertices_with_self_intersects(1:end-1,:),...
    vertices_with_self_intersects(2:end,:),...
    fig_num);

    assert(isequal(round(projected_points,4),[0,0; 0,0; 1,0; 0.75,0.75; 0.6667,1; 0.6667,1; 0.5,1; 0,1]));
