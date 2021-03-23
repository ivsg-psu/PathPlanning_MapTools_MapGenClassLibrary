function [sensing_polytopes] = fcn_polytope_editing_expand_polytopes_for_sensing(polytopes,delta,exp_dist)
% expands polytopes while cutting off corners where radii would be to get
% a conservative estimate of when a polytope would be sensed

% %% size sensing polytopes
% sensing_polytopes = polytopes;
% 
% %% expand poltyopes evenly as normal
% exp_polytopes = fcn_polytope_editing_expand_polytopes_evenly(polytopes,delta,exp_dist);
% 
% %% modify polytopes to add psudo radii
% for polys = 1:size(exp_polytopes,2)
%     % pull values
%     xv = exp_polytopes(polys).xv;
%     yv = exp_polytopes(polys).yv;
%     
%     % reorginize for easier use
%     x = ([xv(end) xv xv(1)]-centroid(1))*scale; % x values to modify
%     y = [yv(end) yv yv(1)]; % y values to modify
%     
%     
%     % modify corners
%     for vert = 1:length(xv)
        
        
%% Check input arguments
if nargin ~= 3
    error('Incorrect number of arguments');
end
if exp_dist < 0
    error('exp_dist cannot be less than 0')
end

%% main code
sensing_polytopes = polytopes; % both structures will be the same size
if exp_dist > 0 % there is an expansion
    for polys = 1:size(polytopes,2) % check each obstacle
        %% pull values
        xvert = polytopes(polys).xv; % pull vertice values
        yvert = polytopes(polys).yv;
        numverts = length(polytopes(polys).xv); % find number of vertices
        xv = zeros(1,numverts*3); % pre-allocate xv and yv 
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
            if dist < 0
                exp_ang = exp_ang+pi;
            end
            xv(vert2*3-2) = polytopes(polys).xv(vert2) + exp_dist*cos(ang1+turn); % change x of vert2
            xv(vert2*3-1) = polytopes(polys).xv(vert2) + exp_dist*cos(exp_ang); % change x of vert2
            xv(vert2*3) = polytopes(polys).xv(vert2) + exp_dist*cos(ang2+turn); % change x of vert2
            yv(vert2*3-2) = polytopes(polys).yv(vert2) + exp_dist*sin(ang1+turn); % change y of vert2
            yv(vert2*3-1) = polytopes(polys).yv(vert2) + exp_dist*sin(exp_ang); % change y of vert2
            yv(vert2*3) = polytopes(polys).yv(vert2) + exp_dist*sin(ang2+turn); % change y of vert2
        end
        %% assign new exp_polytopes values
        sensing_polytopes(polys).vertices = [[xv xv(1)]' [yv yv(1)]'];
        sensing_polytopes(polys).xv = xv;
        sensing_polytopes(polys).yv = yv;
        sensing_polytopes(polys).distances = fcn_general_calculation_euclidean_point_to_point_distance(sensing_polytopes(polys).vertices(1:end-1,:),sensing_polytopes(polys).vertices(2:end,:));
        [Cx,Cy,sensing_polytopes(polys).area] = fcn_polytope_calculation_centroid_and_area([xv xv(1)],[yv yv(1)]);
        sensing_polytopes(polys).mean = [Cx Cy];
        sensing_polytopes(polys).max_radius = max(fcn_general_calculation_euclidean_point_to_point_distance(sensing_polytopes(polys).vertices(1:end-1,:),ones(length(xv),1)*sensing_polytopes(polys).mean));
    end
% else % exp_dist == 0 % no need to do anything
end