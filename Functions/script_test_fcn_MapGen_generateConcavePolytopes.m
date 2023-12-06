clear all; close all; clc;
%% generate polytopes
line_width = 3;
axis_box = [0 1 0 1];
halton_range = [1 100]; % Sets the indices to use from halton set
[polytopes] = fcn_MapGen_haltonVoronoiTiling(halton_range);

%% concavify
concave_polytopes = polytopes;
for i = 1:length(concave_polytopes)
    centroid = concave_polytopes(i).mean;
    xv = concave_polytopes(i).xv(2);
    yv = concave_polytopes(i).yv(2);
    new_xv = xv + (centroid(1) - xv)*(1);
    new_yv = yv + (centroid(2) - yv)*(1);
    concave_polytopes(i).xv(2) = new_xv;
    concave_polytopes(i).yv(2) = new_yv;
    concave_polytopes(i).vertices(2,:)= [new_xv new_yv];
end

%% shrink to radius
des_rad = 0.05; sigma_radius = 0; min_rad = 0.001;
fig_num = 1;
shrunk_polytopes=...
    fcn_MapGen_polytopesShrinkToRadius(...
    polytopes,des_rad,sigma_radius,min_rad,fig_num);
fig_num = fig_num + 1;
shrunk_concave_polytopes=...
    fcn_MapGen_polytopesShrinkToRadius(...
    concave_polytopes,des_rad,sigma_radius,min_rad,fig_num);

%% shrink from edges
% fig_num = 1; % Will default to the next figure, usually 1
% gap_size = 0.05;
% shrunk_polytopes = fcn_MapGen_polytopesShrinkFromEdges(...
%     polytopes,gap_size);
% fig_num = fig_num + 1;
% shrunk_concave_polytopes = fcn_MapGen_polytopesShrinkFromEdges(...
%     concave_polytopes,gap_size);

%% plot

fig_num = fig_num + 1;
fig1 = fcn_MapGen_plotPolytopes(polytopes,fig_num,'r-',line_width);
hold on; box on;

fig_num = fig_num + 1;
fig2 = fcn_MapGen_plotPolytopes(concave_polytopes,fig_num,'m-',line_width);
hold on; box on;

fig_num = fig_num + 1;
fig3 = fcn_MapGen_plotPolytopes(shrunk_polytopes,fig_num,'b-',line_width);
hold on; box on;

fig_num = fig_num + 1;
fig4 = fcn_MapGen_plotPolytopes(shrunk_concave_polytopes,fig_num,'g-',line_width);
hold on; box on;
