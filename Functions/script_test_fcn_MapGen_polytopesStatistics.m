% script_test_fcn_MapGen_polytopesStatistics
% Tests: fcn_MapGen_polytopesStatistics

% 
% REVISION HISTORY:
% 
% 2021_07_12 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§

close all;

% Generate a set of polytopes from the Halton set
fig_num = 12;
Halton_range = [200 1455]; % range of Halton points to use to generate the tiling
polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1]);

fcn_MapGen_polytopesStatistics(...
    polytopes,...
    fig_num)



% Shrink all polytopes by a gap of 0.001
shrinkage = 0.015;
shrunk_polytopes = polytopes;
for ith_poly = 1:length(polytopes)
    orig_radius = polytopes(ith_poly).max_radius;
    des_rad = orig_radius - shrinkage;

    tolerance = 1e-5; % This is the edge distance below which vertices are merged together in the polytope
    shrunk_polytopes(ith_poly) =...
        fcn_MapGen_polytopeShrinkToRadius(...
        polytopes(ith_poly),des_rad,tolerance);
end



% des_rad = 0.03; 
% sigma_radius = 0; 
% min_rad = 0.001;
% shrunk_polytopes2=fcn_MapGen_polytopesShrinkToRadius(...
%     polytopes,des_rad,sigma_radius,min_rad);

fig_num = 13;
fcn_MapGen_polytopesStatistics(...
    shrunk_polytopes,...
    fig_num)