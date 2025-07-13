% script_test_polytope_editing_set_all_costs
% Tests function: polytope_editing_set_all_costs

% REVISION HISTORY:
% 2025_04_28
% -- first written by S. Harnett
fig_num = 1;
bounding_box = [0,0; 1,1];

seedGeneratorNames = 'haltonset';
seedGeneratorRanges = [1 20];
AABBs = [0 0 1 1];
mapStretchs = [1 1];

[polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (-1));

URHERE

polytopes = fcn_MapGen_haltonVoronoiTiling([1 20]);
polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box);

des_cost = 0.1;

new_cost_polytopes = ...
    fcn_polytope_editing_set_all_costs(polytopes, des_cost);

new_costs = extractfield(new_cost_polytopes, 'cost');
assert(isequal(new_costs , des_cost*ones(1,length(polytopes))));
