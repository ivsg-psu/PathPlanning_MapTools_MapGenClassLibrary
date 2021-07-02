function [exp_polytopes] = fcn_polytope_editing_expand_polytopes_evenly(polytopes,delta,exp_dist)
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

%% Check input arguments
if nargin ~= 3
    error('Incorrect number of arguments');
end

%% main code §

exp_polytopes = polytopes; % both structures will be the same size
for polys = 1:size(polytopes,2) % check each obstacle
    %% pull values
    xvert = polytopes(polys).xv; % pull vertice values
    yvert = polytopes(polys).yv;
    numverts = length(polytopes(polys).xv); % find number of vertices
    xv = zeros(1,numverts); % pre-allocate xv and yv 
    yv = xv;
    for vert1 = 1:numverts % repeat for each vert
        %% assign 3 vertex indices
        if vert1 < numverts-1 
            vert2 = vert1 + 1;
            vert3 = vert2 + 1;
        elseif vert1 < numverts
            vert2 = vert1 + 1;
            vert3 = 1;
        else
            vert2 = 1;
            vert3 = vert2 + 1;
        end
        %% find which side of the line is in the obstacle
        ang1 = atan2(yvert(vert2)-yvert(vert1),xvert(vert2)-xvert(vert1)); % angle of first segment
        if vert1 == 1 % only execute once per obstacle
            mid_point = [(xvert(vert2)+xvert(vert1))/2; (yvert(vert2)+yvert(vert1))/2; 0]; % segment midpoint
            % check if turning ccw puts the point in the obstacle
            if inpolygon(mid_point(1)+delta*cos(ang1+pi/2),mid_point(2)+delta*sin(ang1+pi/2),xvert,yvert)
                turn = -pi/2; % turn the other way
            else % if puts point out of the obstacle
                turn = pi/2; % keep turning that way
            end
        end
        ang2 = atan2(yvert(vert3)-yvert(vert2),xvert(vert3)-xvert(vert2)); % angle of second segment
        %% force both angles to be positive
        if ang1 < 0
            ang1 = ang1 + 2*pi;
        end
        if ang2 < 0
            ang2 = ang2 + 2*pi;
        end
        %% move vertices
        exp_ang = (ang1 + ang2 + 2*turn)/2; % direction to move vertex
        inner_ang = exp_ang - ang2 - turn; % angle for determining distance
        dist = exp_dist/cos(inner_ang); % how far to move to ensure both sides move out by exp_dist
        xv(vert2) = polytopes(polys).xv(vert2) + dist*cos(exp_ang); % change x of vert2
        yv(vert2) = polytopes(polys).yv(vert2) + dist*sin(exp_ang); % change y of vert2
    end
    %% assign new exp_polytopes values
    exp_polytopes(polys).vertices = [[xv xv(1)]' [yv yv(1)]'];
    exp_polytopes(polys).xv = xv;
    exp_polytopes(polys).yv = yv;
    exp_polytopes(polys).distances = sum((exp_polytopes(polys).vertices(1:end-1,:)-exp_polytopes(polys).vertices(2:end,:)).^2,2).^0.5;
    [Cx,Cy,exp_polytopes(polys).area] = fcn_MapGen_polytopeCentroidAndArea([xv xv(1)]', [yv yv(1)]');
    exp_polytopes(polys).mean = [Cx, Cy];        
    exp_polytopes(polys).max_radius = max(sum((exp_polytopes(polys).vertices(1:end-1,:)-ones(length(xv),1)*exp_polytopes(polys).mean).^2,2).^0.5);
end

% §
% Debug
%
% Functions §