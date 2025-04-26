%% script_test_fcn_MapGen_generateVoronoiDiagramBetweenPolytopes
% test script of making a Voronoi diagram from points along side sof polytopes

clear; close all;

%% test convex
fig_num = 1;
figure(fig_num);
clf;

Halton_range = [200 301]; % range of Halton points to use to generate the tiling
polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1]);
des_gap_size = 0.025;
[shrunk_polytopes] = fcn_MapGen_polytopesShrinkFromEdges(polytopes, des_gap_size);
is_nonconvex = 0;
[vx,vy,h] = fcn_MapGen_generateVoronoiDiagramBetweenPolytopes(shrunk_polytopes,is_nonconvex);
hold on; box on;
for j = 1:length(shrunk_polytopes)
     fill(shrunk_polytopes(j).vertices(:,1)',shrunk_polytopes(j).vertices(:,2),[0 0 1],'FaceAlpha',0.5)
end
xlabel('x [km]')
ylabel('y [km]')

axis([-0.5 1.5 -0.5 1.5])


%% test nonconvex
fig_num = 2;
figure(fig_num);
clf;

load('flood_plain_4.mat'); % Data is in testFixtures subfolder
shrunk_polytopes = flood_plain_4;
is_nonconvex = 1;
[vx,vy,h] = fcn_MapGen_generateVoronoiDiagramBetweenPolytopes(shrunk_polytopes,is_nonconvex);
hold on; box on;
for j = 1:length(shrunk_polytopes)
     fill(shrunk_polytopes(j).vertices(:,1)',shrunk_polytopes(j).vertices(:,2),[0 0 1],'FaceAlpha',0.5);
end
xlabel('x [km]')
ylabel('y [km]')
