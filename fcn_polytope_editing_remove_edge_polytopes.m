function [trim_polytopes] = fcn_polytope_editing_remove_edge_polytopes(polytopes,xlow,xhigh,ylow,yhigh)
% FCN_POLYTOPE_EDITING_REMOVE_EDGE_POLYTOPES removes polytopes that extend 
% beyond the boundaries specified
%
% [TRIM_POLYTOPES]=FCN_POLYTOPE_EDITING_REMOVE_EDGE_POLYTOPES(POLYTOPES,XLOW,XHIGH,YLOW,YHIGH)
% returns:
% TRIM_POLYTOPES: a 1-by-n seven field structure of polytopes within the 
% boundaries, where n <= number of polytopes with fields:
%   vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%     the number of the individual polytope vertices
%   xv: a 1-by-m vector of vertice x-coordinates
%   yv: a 1-by-m vector of vertice y-coordinates
%   distances: a 1-by-m vector of perimeter distances from one point to the
%     next point, distances(i) = distance from vertices(i) to vertices(i+1)
%   mean: centroid xy coordinate of the polytope
%   area: area of the polytope
%
% with inputs:
% POLYTOPES: the original polytopes with the same fields as trim_polytopes
% XLOW: lower x boundary
% XHIGH: upper x boundary
% YLOW: lower y boundary
% YHIGH: upper y boundary
%
% Examples:
%   cur_path = pwd;
%   main_folder = '!Voronoi Tiling Obstacles - Organized';
%   parent_dir = cur_path(1:strfind(cur_path,main_folder)-2);
%   addpath([parent_dir '\' main_folder '\Plotting'])
%   addpath([parent_dir '\' main_folder '\Map_Generation\polytope_generation'])
%   addpath([parent_dir '\' main_folder '\Map_Generation\polytope_editing'])
%   polytopes = fcn_polytope_generation_halton_voronoi_tiling_obstacles(1,1000);
%   fig = fcn_plot_polytopes(polytopes,[],'b',2,[0 1 0 1]);
%   trim_polytopes = fcn_polytope_editing_remove_edge_polytopes(polytopes,0,1,0,1);
%   fcn_plot_polytopes(trim_polytopes,fig,'g',2,[0 1 0 1]);
%
% This function was written on 2019_06_13 by Seth Tau
% Questions or comments? sat5340@psu.edu 
%

keep = 0; % number of polytopes to keep
trim_polytopes(1) = struct('vertices',[],'xv',[],'yv',[],'distances',[],'mean',[],'area',[],'max_radius',[]);
for poly = 1:size(polytopes,2) % check each polytope within polytopes
    xv = polytopes(poly).xv;
    yv = polytopes(poly).yv;
    if sum((xv<xlow)+(xv>xhigh)+(yv<ylow)+(yv>yhigh))==0 % if the x or y vertices are inside of the bounds
        keep = keep + 1;
        trim_polytopes(keep) = polytopes(poly);
    end
end

% %%%%%%% troubleshooting
% fcn_plot_polytopes(trim_polytopes,[],'b',2,[0 1 0 1]);
% %%%%%%%%%%%%%%%%%%%%%%%