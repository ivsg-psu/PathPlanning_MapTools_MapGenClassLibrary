% script_test_fcn_MapGen_polytopesShrinkToRadius
% Tests function: fcn_MapGen_polytopesShrinkToRadius

% REVISION HISTORY:
% 2021_06_09
% -- first written by S. Brennan using
% script_test_fcn_MapGen_polytopeCropEdges as a template
%% main code §

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
        shrinker,des_rad,tolerance,fig_num);
    pause(0.01);
end

% plot the last output polytope in black
fcn_MapGen_plotPolytopes(shrunk_polytope,fig_num,'k.',2);
  