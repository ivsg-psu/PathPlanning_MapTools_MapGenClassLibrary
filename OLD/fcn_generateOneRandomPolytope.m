function [exp_polytopes] = fcn_generateOneRandomPolytope
% FCN_POLYTOPE_EDITING_EXPAND_POLYTOPES_EVENLY  expands an obstacle out by exp_dist on all sides
%
% [EXP_POLYTOPES]=FCN_POLYTOPE_EDITING_EXPAND_POLYTOPES_EVENLY(POLYTOPES,DELTA,EXP_DIST)
% returns:
% EXP_POLYTOPES: a 1-by-n seven field structure of expanded polytopes, where 
%   n <= NUM_POLY, with fields:
% vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%   the number of the individual polytope vertices
% xv: a 1-by-m vector of vertice x-coordinates
% yv: a 1-by-m vector of vertice y-coordinates
% distances: a 1-by-m vector of perimeter distances from one point to the
%   next point, distances(i) = distance from vertices(i) to vertices(i+1)
% mean: average xy coordinate of the polytope
% area: area of the polytope
% max_radius: distance of the farthest vertex from the mean
%
% with inputs:
% POLYTOPES: the original 1-by-n seven field structure with the same fields
% DELTA: a small number relative to vehicle size to determine the inside of
%   an obstacle
% EXP_DIST: distance to expand the obstacle
%
% Examples:
%      
%      xv = [-2 -1 1 2 2 1 -1 -2];
%      yv = [-1 -2 -2 -1 1 2 2 1];
%      polytopes.vertices = [[xv xv(1)]' [yv yv(1)]'];
%      polytopes.xv = xv;
%      polytopes.yv = yv;
%      polytopes.distances = fcn_general_calculation_euclidean_point_to_point_distance(polytopes(1).vertices(1:end-1,:),polytopes(1).vertices(2:end,:));
%      [Cx,Cy,polytopes.area] = fcn_polytope_calculation_centroid_and_area([xv xv(1)],[yv yv(1)]);
%      polytopes.mean = [Cx, Cy];
%      polytopes.max_radius = max(fcn_general_calculation_euclidean_point_to_point_distance(polytopes.vertices(1:end-1,:),ones(length(xv),1)*polytopes.mean));
%      delta = 0.01;
%      exp_dist = 1;
%      exp_polytopes=fcn_polytope_editing_expand_polytopes_evenly(polytopes,delta,exp_dist);
%      fcn_plot_polytopes(polytopes,99,'r-',2);
%      fcn_plot_polytopes(exp_polytopes,99,'b-',2,[-4 4 -4 4],'square');
%      legend('Original','Expanded')
%      box on
%      xlabel('X Position')
%      ylabel('Y Position')
%
% This function was written on 2018_11_17 by Seth Tau
% Adjusted example code on 2021_04_28 by Seth Tau
% Questions or comments? sat5340@psu.edu 
%


%% main code §
%% Set up variables
polytopes = fcn_MapGen_haltonVoronoiTiling([1 100]);

fig_num = 1;
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);

%% Pick a random polytope
Npolys = length(trim_polytopes);
rand_poly = 1+floor(rand*Npolys);
shrinker = trim_polytopes(rand_poly);

% §
% Debug
%
% Functions §