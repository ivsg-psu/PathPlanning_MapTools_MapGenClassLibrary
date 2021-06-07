function [shrunk_polytopes,mu_final,sigma_final] = fcn_MapGen_polytopeShrinkToRadius(polytopes,des_radius,sigma_radius,min_rad)
% FCN_POLYTOPE_EDITING_SHRINK_TO_AVERAGE_MAX_RADIUS_WITH_VARIANCE shrinks the 
% polytopes to achieve the desired lookout limit with specified variance
%
% [SHRUNK_POLYTOPES,MU_FINAL,SIGMA_FINAL]=FCN_POLYTOPE_EDITING_SHRINK_TO_AVERAGE_MAX_RADIUS_WITH_VARIANCE(POLYTOPES,DES_RADIUS,SIGMA_RADIUS,MIN_RAD)
% returns:
% SHRUNK_POLYTOPES: a 1-by-n seven field structure of shrunken polytopes, 
% where n <= number of polytopes with fields:
%   vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%     the number of the individual polytope vertices
%   xv: a 1-by-m vector of vertice x-coordinates
%   yv: a 1-by-m vector of vertice y-coordinates
%   distances: a 1-by-m vector of perimeter distances from one point to the
%     next point, distances(i) = distance from vertices(i) to vertices(i+1)
%   mean: centroid xy coordinate of the polytope
%   area: area of the polytope
%   max_radius: distance from the mean to the farthest vertex
% MU_FINAL: final average maximum radius achieved
% SIGMA_FINAL: final variance achieved
%
% with inputs:
% POLYTOPES: original polytopes with same fields as shrunk_polytopes
% DES_RAD: desired average max radius   
% SIGMA_RADIUS: desired variance in the radii 
% MIN_RAD: minimum acceptable radius
%
% Examples:
%
% cur_path = pwd;
% main_folder = '!Voronoi Tiling Obstacles - Organized';
% parent_dir = cur_path(1:strfind(cur_path,main_folder)-2);
% addpath([parent_dir '\' main_folder '\Plotting'])
% addpath([parent_dir '\' main_folder '\General_Calculation'])
% addpath([parent_dir '\' main_folder '\Map_Generation\polytope_generation'])
% addpath([parent_dir '\' main_folder '\Map_Generation\polytope_editing'])
% addpath([parent_dir '\' main_folder '\Map_Generation\polytope_calculation'])
% low_pt = 1; high_pt = 100;
% polytopes = fcn_polytope_generation_halton_voronoi_tiling(low_pt,high_pt);
% xlow = 0; xhigh = 1; ylow = 0; yhigh = 1;
% trim_polytopes = fcn_polytope_editing_remove_edge_polytopes(polytopes,xlow,xhigh,ylow,yhigh);
% fcn_plot_polytopes(trim_polytopes,[],'b',2,[0 1 0 1]);
% % uniform shrinking
% des_rad = 0.05; sigma_radius = 0; min_rad = 0.001;
% shrunk_polytopes1=fcn_polytope_editing_shrink_to_average_max_radius_with_variance(trim_polytopes,des_rad,sigma_radius,min_rad);
% fcn_plot_polytopes(shrunk_polytopes1,[],'b',2,[0 1 0 1]);
% % non-uniform shrinking
% scale = ones(1,size(trim_polytopes,2));
% des_rad = 0.05; sigma_radius = 0.01; min_rad = 0.001;
% shrunk_polytopes2=fcn_polytope_editing_shrink_to_average_max_radius_with_variance(trim_polytopes,des_rad,sigma_radius,min_rad);
% fcn_plot_polytopes(shrunk_polytopes2,[],'b',2,[0 1 0 1]);
%   
% This function was written on 2019_08_29 by Seth Tau
% Questions or comments? sat5340@psu.edu 
%


%% find current distribution
radii = [polytopes.max_radius];

%%% trouble shooting
% fcn_plot_polytopes(polytopes,[],'b',2);
% 
% figure
% histogram(radii,20)
% r_sigma = std(radii);
%%%

r_mu = mean(radii);
r_size = length(radii);

if r_mu < des_radius
    error('cannot achieve the desired radius')
end

%% determine desired distribution
r_dist = normrnd(des_radius,sigma_radius,[r_size,1]);

%%% troubleshoooting
% figure
% histogram(r_dist,20)
%%%

r_dist = r_dist + (des_radius-mean(r_dist)); % adjust to ensure the mean value is mu
max_r_dist = max(radii); % largest possible radius
min_r_dist = min_rad; % smallest possible radius
if sum((r_dist>max_r_dist)+(r_dist<min_r_dist)) > 0
    warning('standard deviation skewed due to truncated distribution')
end
r_dist(r_dist>max_r_dist) = max_r_dist; % truncate any values that are too large
r_dist(r_dist<min_r_dist) = min_r_dist; % truncate any values that are too small
while abs(mean(r_dist) - des_radius) > 1e-10
    r_dist = r_dist + (des_radius-mean(r_dist));
    r_dist(r_dist>max_r_dist) = max_r_dist;
    r_dist(r_dist<min_r_dist) = min_r_dist;
end

mu_final = mean(r_dist);
sigma_final = std(r_dist);

%%% troubleshoooting
% figure
% histogram(r_dist,20)
%%%

%% shrink polytopes to achieve the distribution
[new_rads,ob_ind] = sort(r_dist);
if sum((sort(radii)'-sort(r_dist))>=-2*min_rad) < r_size
    error('distribution is unachievable with generated map')
end

shrunk_polytopes = polytopes;
for idx = 1:length(new_rads)
    shrinker = polytopes(ob_ind(idx)); % obstacle to be shrunk
    
    % pull values
    vertices = shrinker.vertices;
    centroid = shrinker.mean;
    rad = shrinker.max_radius;
    
    % determine scale factor
    scale = new_rads(idx)/rad;
    
    if scale < 1 % calculation error can sometimes make this greater than 1
        % find new vertices
        new_vert = centroid + scale*(vertices-centroid);

        % adjust polytopes
        shrinker.vertices = new_vert;
        shrinker.xv = new_vert(1:end-1,1)';
        shrinker.yv = new_vert(1:end-1,2)';
        shrinker.distances = fcn_geometry_euclideanPointsToPointsDistance(new_vert(1:end-1,:),new_vert(2:end,:));
        shrinker.area = shrinker.area*scale^2;
        shrinker.max_radius = shrinker.max_radius*scale;
    end
    
    % assign to shrunk_polytopes
    shrunk_polytopes(ob_ind(idx)) = shrinker;
end


