% script_test_fcn_MapGen_calculateConvexHullOverlapRatio
% Tests: fcn_MapGen_calculateConvexHullOverlapRatio

%
% REVISION HISTORY:
%
% 2024_02_28 by Steve Harnett
% -- first write of script
% 2025_04_16 by Steve Harnett
% -- remove dependence on test fixture
%%%%%%%%%%%%%%§

close all; 

%% basic polytope case
fig_num = 1;
figure(fig_num);
clf;

polytopes(1).vertices = [0 0; 10 0; 10 1; 0 1; 0 0];
polytopes(2).vertices = polytopes(1).vertices+[8,0];
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes);

[ ...
    convex_hull_overlap_ratio,...
    A_overlap,...
    A_occupied...
    ] = ...
    fcn_MapGen_calculateConvexHullOverlapRatio( ...
    polytopes, ...
    fig_num...
    );

assert(isequal(0.1000            ,round(convex_hull_overlap_ratio,4)))
assert(isequal(2                 ,round(A_overlap,4)))
assert(isequal(20                ,round(A_occupied,4)))


%% 
fig_num = 2;
figure(fig_num);
clf;

polytopes(1).vertices = [0 0; 5 0; 7, 0.5; 5 1; 0 1; 0 0];
polytopes(2).vertices = [6 0; 10 0; 10 1; 6 1; 8 0.5; 6 0];
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes,1,1009);
fig_num = fig_num + 1;
[ ...
    convex_hull_overlap_ratio,...
    A_overlap,...
    A_occupied...
    ] = ...
    fcn_MapGen_calculateConvexHullOverlapRatio( ...
    polytopes, ...
    fig_num...
    );

assert(isequal(0.0278            ,round(convex_hull_overlap_ratio,4)))
assert(isequal(0.25              ,round(A_overlap,4)))
assert(isequal(9                 ,round(A_occupied,4)))

