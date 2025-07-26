function [path,distance]=path_CalculatePathMethodDStar(m,n,map,startPoint,quadrantPos,goalPoint)
%[path]=path_CalculatePathMethodDStar(m, n, map, startPoint, quadrantPos,
% goalPoint) generates a path from the startPoint to goalPoint using the
% DStar method. The input and output variables are defined as follows:
%
% PATH  a nx2 vector of i,j indeces in the matrix to visit that represents
% the best path.
%
% N,M   the numbers of rows, columns in the map array, respectively.
%
% MAP   the N-by-M matrix containing 1's and 0's, with 1 being occupied and
% 0 clear. 
%
% STARTPOINT  the starting position of the robot.
%
% QUADRANT_POS   ?.
% 
% GOALPOINT  the goal, e.g. where the robot is trying to get to.
% 
% Examples:
%        Set size of array... change to something small (20 by 30) if debug mode!
%        M = 300;
%        N = 300;
% 
%        % Fill in initial random map
%        occupancy_map = fcn_generate_random_occupancy(M,N);
%        occupancy_map = +occupancy_map;
%        image(occupancy_map+1);
%        colormap([1 1 1;0 0 0]);
% 
%        % Fill in the start position
%        goodpoint = false;
%        while goodpoint==false
%           startPointrow = ceil(rand*M);
%            startPointcol = ceil(rand*N);
%            if occupancy_map(startPointrow,startPointcol)==0
%                goodpoint = true;
%            end
%        end
%        startPoint = [startPointrow startPointcol];
%        text(startPointcol,startPointrow,sprintf('S'),'Color',[0 1 0]); hold on;
% 
%        % Fill in the goal position
%        goodpoint = false;
%        while goodpoint==false
%            goalPointrow = ceil(rand*M);
%            goalPointcol = ceil(rand*N);
%            if occupancy_map(goalPointrow,goalPointcol)==0
%                goodpoint = true;
%            end
%        end
%        goalPoint = [goalPointrow goalPointcol];
%        text(goalPointcol,goalPointrow,sprintf('G'),'Color',[1 0 0]); hold on;
% 
%        quadrantPos = 0;
% 
%        % Call function to calculate path
%        [path]=path_CalculatePathMethodDstar(M,N,occupancy_map,startPoint,quadrantPos,goalPoint);
% 
%        % Plot map results
%        tempmap = occupancy_map;
%        for i=1:length(path(:,1))
%            tempmap(path(i,1),path(i,2)) = 2;
%        end
%        image(tempmap+1);
%        colormap([1 1 1;0 0 0;1 0 0]); hold on;
%        text(startPoint(1,2),startPoint(1,1),sprintf('S'),'Color',[0 0 0]); 
%        text(goalPoint(1,2),goalPoint(1,1),sprintf('G'),'Color',[1 0 0]);
%
% This function was written on 2009_03_26 by Dr. Sean Brennan
% Questions or comments? sbrennan@psu.edu 

do_debug = 0;

% Fill in the first point
path(1,:)=startPoint;

% Not sure what the following code does?
% searchForGoalPoint=false;
% goalPointCheck=goalPoint-(quadrantPos-1)*100;
% if goalPointCheck>0
%     if goalPointCheck<300
%         goalPoint=goalPointCheck;
%         searchForGoalPoint=true;
%     end
% end

%% Initialize distance map
% FOR DEBUGGING: Initialize distance matrix to -1 (will use -1 for unexplored)
distance_map = map-1*ones(m,n);
distance_map(goalPoint(1),goalPoint(2))=1;

if do_debug == 1
    plotit(map,distance_map,startPoint,goalPoint,path);
end
%% Initialize neighbors
distance = 1;
[neighbors,distance_map] = find_neighbors(goalPoint,m,n,distance,distance_map);
neighbors = sortrows(neighbors,1);

if do_debug == 1
    plotit(map,distance_map,startPoint,goalPoint,path);
end
    
found_start = false;
while ~isempty(neighbors)
    old_neighbors = neighbors;
    neighbors = [];
    for i=1:length(old_neighbors(:,1))
        current_point = [old_neighbors(i,2) old_neighbors(i,3)];
        [new_neighbors,distance_map] = find_neighbors(current_point,m,n,old_neighbors(i,1),distance_map);
        neighbors = [neighbors; new_neighbors]; %#ok<AGROW>
    end

    if do_debug == 1
        % FOR DEBUGGING:
        plotit(map,distance_map,startPoint,goalPoint,path);
    end
    
    neighbors = sortrows(neighbors,1);
    
    i = 0;
    while i<length(neighbors)
        i=i+1;
        if all(neighbors(i,2:3)==startPoint)
            neighbors = [];
            found_start = true;
        end
    end
end    

% When the D-star is done, find best path back to goal and save
if found_start == true
    % Fill in the first point
    path(1,:)=startPoint;
    while ~all(path(end,:)==goalPoint)
        best_distance = distance_map(path(end,1),path(end,2));
        best_pos_i = path(end,1);
        best_pos_j = path(end,2);
        for i=max(path(end,1)-1,1):min(path(end,1)+1,m)  % Loop through all rows, from min to max, avoiding edges
            for j=max(path(end,2)-1,1):min(path(end,2)+1,n)  % Loop through all columns, from min to max, avoiding edges
                if (distance_map(i,j)~=-1) && (distance_map(i,j)~=0)
                    if distance_map(i,j)<best_distance
                        best_distance = distance_map(i,j);
                        best_pos_i = i;
                        best_pos_j = j;
                    end
                end
            end
        end
        path = [path; best_pos_i best_pos_j]; %#ok<AGROW>
    end
    distance = distance_map(startPoint(1),startPoint(2));
end

if do_debug == 1
    plotit(map,distance_map,startPoint,goalPoint,path);
end

% path(:,1)=path(:,1)+ones(length(path),1).*(quadrantPos(1)-1)*100;
% path(:,2)=path(:,2)+ones(length(path),1).*(quadrantPos(2)-1)*100;

end

function [neighbors,distance_map] = find_neighbors(Point,m,n,distance,distance_map)
neighbors = zeros(8,3);
count = 0;
for i=max(Point(1)-1,1):min(Point(1)+1,m)  % Loop through all rows, from min to max, avoiding edges
    for j=max(Point(2)-1,1):min(Point(2)+1,n)  % Loop through all columns, from min to max, avoiding edges
        if distance_map(i,j)==-1  % This is an undiscovered point
            % Fill in the distance
            if i~=Point(1)&&j~=Point(2)  % This is in a cell on a diagonal
                cell_dist = distance+1.4142;
            else
                cell_dist = distance+1;
            end
            distance_map(i,j) = cell_dist;

            %Update the neighbors list
            count = count+1;
            neighbors(count,:) = [cell_dist i j];
        end
    end
end

if count>0
    neighbors = neighbors(1:count,:);
else
    neighbors = [];
end


% % The following works, but takes 71% of time according to profiler
% neighbors = [];
% for i=max(Point(1)-1,1):min(Point(1)+1,m)  % Loop through all rows, from min to max, avoiding edges
%     for j=max(Point(2)-1,1):min(Point(2)+1,n)  % Loop through all columns, from min to max, avoiding edges
%         if distance_map(i,j)==-1  % This is an undiscovered point
%             % Fill in the distance
%             if i~=Point(1)&&j~=Point(2)  % This is in a cell on a diagonal
%                 cell_dist = distance+1.4142;
%             else
%                 cell_dist = distance+1;
%             end
%             distance_map(i,j) = cell_dist;
%             
%             %Update the neighbors list
%             neighbors = [neighbors; [cell_dist i j]]; %#ok<AGROW>
%         end
%     end
% end

end


function plotit(map, distance_map,startPoint,goalPoint,path)
tempmap = map;
for i=1:length(path(:,1))
    tempmap(path(i,1),path(i,2)) = 2;
end
image(tempmap+1);
colormap([1 1 1;0 0 0;1 0 0]); hold on;
for i=1:length(distance_map(:,1))
    for j=1:length(distance_map(1,:))
        if map(i,j)==0
            if isequal([i j],startPoint)
                text(j,i,sprintf('S'),'Color',[0 0 0]); 
            elseif isequal([i j],goalPoint)
                text(j,i,sprintf('G'),'Color',[1 0 0]);
            elseif distance_map(i,j)==-1
              %  text(j,i,sprintf('%d',floor(distance_map(i,j))),'Color',[1 0 1]); hold on;
            else
                text(j,i,sprintf('%d',floor(distance_map(i,j))),'Color',[0 0 1]);
            end
        end
    end
end
end