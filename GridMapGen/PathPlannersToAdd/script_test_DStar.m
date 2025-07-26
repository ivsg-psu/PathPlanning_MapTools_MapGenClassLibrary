%script_test_DStar

% Set size of array... change to something small (20 by 30) if debug mode!
M = 300;
N = 300;

% Fill in initial random map
occupancy_map = fcn_generate_random_occupancy(M,N);
occupancy_map = +occupancy_map;
   image(occupancy_map+1);
colormap([1 1 1;0 0 0]);

% Fill in the start position
goodpoint = false;
while goodpoint==false
    startPointrow = ceil(rand*M);
    startPointcol = ceil(rand*N);
    if occupancy_map(startPointrow,startPointcol)==0
        goodpoint = true;
    end
end
startPoint = [startPointrow startPointcol];
text(startPointcol,startPointrow,sprintf('S'),'Color',[0 1 0]); hold on;

% Fill in the goal position
goodpoint = false;
while goodpoint==false
    goalPointrow = ceil(rand*M);
    goalPointcol = ceil(rand*N);
    if occupancy_map(goalPointrow,goalPointcol)==0
        goodpoint = true;
    end
end
goalPoint = [goalPointrow goalPointcol];
text(goalPointcol,goalPointrow,sprintf('G'),'Color',[1 0 0]); hold on;

quadrantPos = 0;

% Call function to calculate path
[path]=path_CalculatePathMethodDstar(M,N,occupancy_map,startPoint,quadrantPos,goalPoint);

% Plot map results
tempmap = occupancy_map;
for i=1:length(path(:,1))
    tempmap(path(i,1),path(i,2)) = 2;
end
image(tempmap+1);
colormap([1 1 1;0 0 0;1 0 0]); hold on;
text(startPoint(1,2),startPoint(1,1),sprintf('S'),'Color',[0 0 0]); 
text(goalPoint(1,2),goalPoint(1,1),sprintf('G'),'Color',[1 0 0]); 