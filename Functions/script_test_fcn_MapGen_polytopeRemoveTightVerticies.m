% script_test_fcn_MapGen_polytopeRemoveTightVerticies
% Tests: fcn_MapGen_polytopeRemoveTightVerticies

%
% REVISION HISTORY:
%
% 2021_07_02 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§




%% Set up variables
fig_num = 1;
polytopes = fcn_MapGen_haltonVoronoiTiling([1 100],[1 1], fig_num);


bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);

%% Pick a random polytope
Npolys = length(trim_polytopes);
rand_poly = 1+floor(rand*Npolys);
shrinker = trim_polytopes(rand_poly);

%% Basic example of uniform shrinking
fig_num = 11;
orig_radius = shrinker.max_radius;
ratio = 0.5;
des_rad = orig_radius*ratio;
tolerance = 1e-5;
shrunk_polytope =...
    fcn_MapGen_polytopeShrinkToRadius(...
    shrinker,des_rad,tolerance,fig_num);

%% Iterative example of uniform shrinking
fig_num = 2;
orig_radius = shrinker.max_radius;
ratios = (0.99:-0.05:0);

for ith_ratio = 1:length(ratios)
    des_rad = orig_radius*ratios(ith_ratio);
    tolerance = 1e-5;
    shrunk_polytope =...
        fcn_MapGen_polytopeShrinkToRadius(...
        shrinker,des_rad,tolerance,fig_num);
    pause(0.01);
end

%% Show results of increasing the tolerance distance to merge points.
% So that points merge earlier than they would, thus allowing one to see
% the effects of merging points around the polytope.
fig_num = 3;
figure(fig_num);
clf;
orig_radius = shrinker.max_radius;
ratios = (0.99:-0.05:0);

for ith_ratio = 1:length(ratios)
    des_rad = orig_radius*ratios(ith_ratio);
    tolerance = 0.02;
    shrunk_polytope =...
        fcn_MapGen_polytopeShrinkToRadius(...
        shrinker,des_rad,tolerance);
    cleaned_polytope = fcn_MapGen_polytopeRemoveTightVerticies(...
        shrunk_polytope, tolerance,fig_num);
    pause(0.01);
end

% make some assertion tests based on expectations of the last cleaned polytope
assert(isequal(round(cleaned_polytope.distances,4),[0;0;0]));
assert(isnan(cleaned_polytope.mean(1)));
assert(isequal(cleaned_polytope.area,0));
assert(isnan(cleaned_polytope.max_radius));

% plot the last output polytope in black
fcn_MapGen_plotPolytopes(shrunk_polytope,fig_num,'k.',2);


