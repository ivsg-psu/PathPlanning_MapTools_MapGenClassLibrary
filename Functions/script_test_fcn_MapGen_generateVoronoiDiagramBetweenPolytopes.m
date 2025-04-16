clear; close all; clc
% script_test_fcn_MapGen_generateVoronoiDiagramBetweenPolytopes
% test script of making a Voronoi diagram from points along side sof polytopes

addpath(strcat(pwd,'\..\..\PathPlanning_PathTools_PathClassLibrary\Functions'));
addpath(strcat(pwd,'\..\..\PathPlanning_MapTools_MapGenClassLibrary\Functions'));
addpath(strcat(pwd,'\..\..\Errata_Tutorials_DebugTools\Functions'));

%% test convex
Halton_range = [200 301]; % range of Halton points to use to generate the tiling
polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1]);
des_gap_size = 0.025;
[shrunk_polytopes] = fcn_MapGen_polytopesShrinkFromEdges(polytopes, des_gap_size);
is_nonconvex = 0;
[vx,vy,h] = fcn_MapGen_generateVoronoiDiagramBetweenPolytopes(shrunk_polytopes,is_nonconvex)
hold on; box on;
for j = 1:length(shrunk_polytopes)
     fill(shrunk_polytopes(j).vertices(:,1)',shrunk_polytopes(j).vertices(:,2),[0 0 1],'FaceAlpha',0.5)
end
xlabel('x [km]')
ylabel('y [km]')


%% test nonconvex
load(strcat(pwd,'\..\testFixtures\flood_plain_4.mat'))
shrunk_polytopes = flood_plain_4;
is_nonconvex = 1;
[vx,vy,h] = fcn_MapGen_generateVoronoiDiagramBetweenPolytopes(shrunk_polytopes,is_nonconvex)
hold on; box on;
for j = 1:length(shrunk_polytopes)
     fill(shrunk_polytopes(j).vertices(:,1)',shrunk_polytopes(j).vertices(:,2),[0 0 1],'FaceAlpha',0.5)
end
xlabel('x [km]')
ylabel('y [km]')
