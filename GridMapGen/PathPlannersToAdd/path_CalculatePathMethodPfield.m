function [path]=path_CalculatePathMethodPfield(map,actualPos,goalPoint)
% by Joan Singla - jzs207@psu.edu

[m,n]=size(map);
path(1,:)=actualPos;
p=1;
endPathPlaner=false;
removeLake=false;

%% PREPARE PFIELD MAP

% create the smoothing filters (left and right) generator (just first time)
left_filter=diag(0.5*ones(m,1))+diag(0.25*ones(m-1,1),1)+diag(0.25*ones(m-1,1),-1);
left_filter(1,1)=0.75;
left_filter(m,m)=0.75;
right_filter=diag(0.5*ones(n,1))+diag(0.25*ones(n-1,1),1)+diag(0.25*ones(n-1,1),-1);
right_filter(1,1)=0.75;
right_filter(n,n)=0.75;
% create the slope from a start point to the end point
[x,y]=meshgrid(1:m,1:n);
U=(x-goalPoint(2)).^2+(y-goalPoint(1)).^2;
U=U/max(max(U));           % values between 0 and 1


order=14;           % times the filter to be applied
map=left_filter^order*map*right_filter^order;
map=map+U;


%% CALCULATE THE PATH

while endPathPlaner==false
    searchZone=[path(p,1)-1 path(p,1)+1;path(p,2)-1 path(p,2)+1];

    if searchZone(1,1)<1
        searchZone(1,1)=1;
    elseif searchZone(1,2)>n
        searchZone(1,2)=n;
    end
    if searchZone(2,1)<1
        searchZone(2,1)=1;
    elseif searchZone(2,2)>m
        searchZone(2,2)=m;
    end

    [minXCols,posRowsMin]=min(map(searchZone(1,1):searchZone(1,2),searchZone(2,1):searchZone(2,2)));
    [valueMin,posColMin]=min(minXCols);

    if posRowsMin(posColMin)==2 && posColMin==2
        if path(p,1)==goalPoint(1) && path(p,2)==goalPoint(2)
            endPathPlaner=true;
        else %minimum local reached - time to fill it
            map(path(p,1),path(p,2))=1;
            removeLake=true;
        end
    elseif path(p,1)==2 || path(p,1)==m-1 || path(p,2)==2 || path(p,2)==n-1
        endPathPlaner=true;    
    else
        p=p+1;
        path(p,2)=path(p-1,2)+posColMin-2;
        path(p,1)=path(p-1,1)+posRowsMin(posColMin)-2;
    end

end

%% REMOVE THE LAKES CREATED BY THE WATER FILLING EFFECT

if removeLake==true
    k=length(path);
    mapWithPath=zeros(m,n);
    for i=1:k
        mapWithPath(path(i,1),path(i,2))=i;
    end
    q=1;
    newPath(q,:)=path(1,:);
    
    while ~isequal(newPath(q,:),path(k,:))
        searchZone=[newPath(q,1)-1 newPath(q,1)+1;newPath(q,2)-1 newPath(q,2)+1];

        [maxXCols,posRowsMax]=max(mapWithPath(searchZone(1,1):searchZone(1,2),searchZone(2,1):searchZone(2,2)));
        [valueMax,posColMax]=max(maxXCols);

        q=q+1;
        newPath(q,2)=newPath(q-1,2)+posColMax-2;
        newPath(q,1)=newPath(q-1,1)+posRowsMax(posColMax)-2;
    end

    path=newPath;
end


%% PLOT RESULTS

% for i=1:length(path)
%     map(path(i,1),path(i,2))=1.5;    
% end
% figure
% pcolor(map);pause(0.1);

%    %Checks if its tring to search out of the matrix (just for DEBUGGING
%    with maps that don't have barriers)
%     [m,n]=size(map);
%     if searchZone(1,1)<1
%         searchZone(1,1)=1;
%     elseif searchZone(1,2)>n
%         searchZone(1,2)=n;
%     end
%     if searchZone(2,1)<1
%         searchZone(2,1)=1;
%     elseif searchZone(2,2)>m
%         searchZone(2,2)=m;
%     end