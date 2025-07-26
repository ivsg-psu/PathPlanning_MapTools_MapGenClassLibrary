% Example script demonstrating how to get timing of a fucntion
clc;
close all;
fig_num = 1;
figure(fig_num); clf;

% Initialize vector containing recorded times.
times = [];
path=[];
showSteps=true;
mapWait=load('mapWait.txt');

% Input variables
startPoint=[10,10];
goalPoint=[80,80];

colormap([1 1 1; 0 0 0;1 0 0]);

for i=1:5
    fileName=strcat('Map',num2str(i),'.txt');
    map=load(fileName);

    if showSteps==true
        %% SHOW map and the diferent paths calculated
        image(map+1);
        title(sprintf('CASE MAP # %d/5   - press any key to continue -',i));
        figure(fig_num);
        pause

        % POTENTIAL FIELD
        path=path_CalculatePathMethodPfield(map,startPoint,goalPoint);
        mapPlot=map;
        for j=1:length(path)
            mapPlot(path(j,1),path(j,2))=2;
        end
        figure(fig_num);
        image(mapPlot+1)
        title(sprintf('P.FIELD METHOD - press any key to continue - map # %d/5',i));
        pause

        % A STAR
        path=path_CalculatePathMethodAstar(map,startPoint,goalPoint);
        mapPlot=map;
        for j=1:length(path)
            mapPlot(path(j,1),path(j,2))=2;
        end
        figure(fig_num);
        image(mapPlot+1);
        title(sprintf('A* METHOD - press any key to continue - map # %d/5',i));
        pause
        
        % D STAR
        m = length(map(:,1));
        n = length(map(1,:));
        quadrantPos = 0;
        [path,distance]=path_CalculatePathMethodDStar(m,n,map,startPoint,quadrantPos,goalPoint);
        mapPlot=map;
        for j=1:length(path)
            mapPlot(path(j,1),path(j,2))=2;
        end
        figure(fig_num);
        image(mapPlot+1);
        title(sprintf('D* METHOD - press any key to continue - map # %d/5',i));
        pause;
        
    end

    helpdlg('Timing all functions now');
    
    %% TIMING THE FUNCTIONS
    f1 = @() path_CalculatePathMethodPfield(map,startPoint,goalPoint);
    f2 = @() path_CalculatePathMethodAstar(map,startPoint,goalPoint);
    f3 = @() path_CalculatePathMethodDStar(m,n,map,startPoint,quadrantPos,goalPoint);

    times(i,1) = timeit(f1); %#ok<AGROW>
    times(i,2) = timeit(f2); %#ok<AGROW>
    times(i,3) = timeit(f3); %#ok<AGROW>

end

figure(fig_num);
bar(times)
title 'COMPARISON PATH PLANNERS'
legend('Potential Field','A star','D star')


