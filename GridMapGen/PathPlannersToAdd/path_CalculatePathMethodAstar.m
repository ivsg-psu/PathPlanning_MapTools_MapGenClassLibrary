function path=path_CalculatePathMethodAstar(map,start,goal)
%function [path,open,closed]=astar_intro(map,start,goal)


%% OPEN and CLOSED sets creation

    % open is a 5 cols array of priority queue (xpos, ypos, g, rank, parentx, parenty)

    open(1,:)=[start(1),start(2),0,1.41421356237*min(abs(goal(1)-start(1)),abs(goal(2)-start(2)))+abs((abs(goal(1)-start(1))-abs(goal(2)-start(2)))),0,0];
    closed=zeros(0,6);
    lowest=1;
    sizemap=size(map);
    isinopen=false;
    isinclosed=false;
    index=0;
    
    coef=0.01;
    
    while open(lowest,1) ~= goal(1) || open(lowest,2) ~= goal(2)      % while lowest rank in OPEN is not the goal...

            current=open(lowest,:);                                 % current = remove lowest rank item from OPEN
            closed(end+1,:)=current(:);                             % add current to closed
            %open=removerows(open,lowest);
            open(lowest,:)=[];

        for i=current(1)-1:1:current(1)+1                       % for neighbors of current...
            for j=current(2)-1:1:current(2)+1
                if i~=current(1)||j~=current(2)

                    if map(i,j)==1 || i>sizemap(1)-1 || j>sizemap(2)-1 || i<2 || j< 2      %if neighbor is on an obstacle or outside the map, cost = infinite
                        cost=Inf;
                    else
                        cost=current(3)+sqrt((current(1)-i)^2+(current(2)-j)^2);
                    end

                    [a b]=size(open);
                    to_remove=[];
                    for k1=a:-1:1                        % if neighbor is in open
                        if i==open(k1,1)&&j==open(k1,2)
                            isinopen=true;
                            if cost<open(k1,3)
                                isinopen=false;
                                open(k1,:)=[];
                                %open=removerows(open,k1);
                                break
                            end
                        end
                    end


                    [a b]=size(closed);
                    for k1=a:-1:1                    % if neighbor is in closed
                        if i==closed(k1,1)&&j==closed(k1,2)
                            isinclosed=true;
                            if cost<closed(k1,3)
                                closed(k1,:)=[];
                                %closed=removerows(closed,k1);
                                isinclosed=false;
                                break
                            end
                        end
                    end


                    if ~isinopen && ~isinclosed                 % if neighbor is not in open neither in closed
                        [a b]=size(open);
                        k=a;
                        open(k+1,1)=i;
                        open(k+1,2)=j;
                        open(k+1,3)=cost;
                        open(k+1,4)=coef*cost+sqrt((goal(1)-i)^2+(goal(2)-j)^2);
                        open(k+1,5)=current(1);
                        open(k+1,6)=current(2);
                    end
                    isinopen=false;
                    isinclosed=false;
                end
            end
        end

        [XXX,lowest]=min(open(:,4));

    end     %ends the while loop
    
    closed(1,:)=[];
    %removerows(closed,1);


%% path vector creation

    path(1,:)=goal;
    for k=length(open):-1:1
        if open(k,1:2)==goal
            path(end+1,:)=open(k,5:6);
            break
        end
    end

    for k=length(closed):-1:1
        if closed(k,1:2)==path(end,:)
            path(end+1,:)=closed(k,5:6);
        end
        if path(end,:)==start
            break
        end
    end
    
    path=flipud(path);

    
%% visualization
% path_map=zeros(size(map));
% open_map=path_map;
% closed_map=path_map;
% 
% for i=1:length(path) 
%     open_map(path(i,1),path(i,2))=1; 
% end
% for i=1:length(open) 
%     path_map(open(i,1),open(i,2))=1; 
% end
% for i=1:length(closed) 
%     closed_map(closed(i,1),closed(i,2))=1; 
% end
% figure
% image(map+2*path_map+3*open_map+4*closed_map+1);pause(0.1);
% colormap([1 1 1; 0 0 0; 1 0 0;0 0.8 0;0 0.5 0]);
% end

