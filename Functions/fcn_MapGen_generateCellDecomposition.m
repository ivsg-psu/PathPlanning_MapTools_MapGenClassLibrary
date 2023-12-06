% cells = fcn_MapGen_generateCellDecomposition(polytopes)
    % sweep vertical bar from left to right
    % generate all points
    % find left most point
    % does this polytope have points both left and right of this point?
    % if not, draw line all the way up and down
    % if so, are the other points above or below this one?
    % if above, draw line below
    % if below, draw line above

    % if we have a list of sides, we can draw until the nearest side

clear;

n = input('Input the number of vertices of the 2D configuration space:');
config_sp = zeros(n,2);
axis([0,10,0,10]);
daspect([1 1 1]);
for i=1:n
[x,y] = ginput(1);

if ismember(round(x,2),round(config_sp(:,1),2))

x=x +double(0.5);
end
% for j=1:size(config_sp,1)
%
% if x == config_sp(j,1)
% x=x +double(0.5);
%
%
% end
%
% end
config_sp(i,:) = [x,y];
hold on;
plot(x,y,'rx');
hold off;

end

hold on;
plot(polyshape(config_sp(:,1),config_sp(:,2)));
hold off;
Space=config_sp;
vertices = config_sp;

n_o = input('Input the number of obstacles: ');
%obstacles={};
O={};
ov=[];
ovs=[];
ove=[];
for i=1:n_o
n_oi = input(strcat('Input the number of vertices in', ' obstacle', num2str(i)));
obstacle = zeros(n_oi,2);
for j=1:n_oi
[x,y] = ginput(1);
if j==1
ovs=[ovs;x,y]
end

if j==n_oi
ove=[ove;x,y];

end
%obstacle(j,:) = [x,y];

%vertices = [vertices;x,y];
if ismember(round(x,2),round(vertices(:,1),2))

x=x +double(0.5);
end
% for m=1:size(vertices,1)
%
% if x == vertices(m,1)
% x=x +0.5;
%
% end
%
% end

ov=[ov;x,y];
vertices = [vertices;x,y];
obstacle(j,:) = [x,y];
hold on;
plot(x,y,'bx');
hold off;


end
obs(i,:)={obstacle};
hold on;
pg=plot(polyshape(obstacle(:,1),obstacle(:,2)));
pg.FaceColor = 'c';
hold off;
%obstacles ={obstacles;obstacle};
end
tic
svertices =sortrows(vertices);
incell={};
save1=[];
for i=1:size(svertices,1)
vi = svertices(i,:);
% for j=1:size(svertices,1)
% vj = svertices(j,:);
% if i~=j
%
% end
% end


poly1 = polyshape(Space(:,1),Space(:,2));
lineseg = [vi(1) -10;vi(1) 10];
edge=[vi(1) -10 vi(1) 10];
in=intersectEdgePolygon(edge,Space);
in=round(in,2);
size123=size(in);
if size123(1)==1 |(in(1,2)==in(2,2) & size123(1)==2)
in=[];
end
in=unique(in,'rows')
%[in,out] = intersect(poly1,lineseg);
incell(i,:)={in};
%plot(poly1)
hold on
%plot(in(:,1),in(:,2),'b',out(:,1),out(:,2),'r')
%plot(in(:,1),in(:,2),'b')
hold off
poly2=polyshape();
for m=1:n_o
y=obs{m,:};
polyx = polyshape(y(:,1),y(:,2));
poly2=union(poly2,polyx)
end
for m=1:n_o
y1=obs{m,:};
lineseg = [vi(1) -10;vi(1) 10];
edge1=[vi(1) -10 vi(1) 10];
in1=intersectEdgePolygon(round(edge1,2),round(y1,2));
if isempty(in1)
in1=save1;
continue;
else
in1=round(in1,2);
size1234=size(in1);
if size1234(1)==1 |(in1(1,2)==in1(2,2) & size1234(1)==2)


in1=[];
end
in1=unique(in1,'rows');
%[in,out] = intersect(poly1,lineseg);
save1=in1;
incell1(i,:)={in1};
end
end


%[in1,out2] = intersect(poly2,lineseg);
%incell1(i,:)={in1};
hold on
%plot(in1(:,1),in1(:,2),'c')
hold off
if isempty(in1)
inner=1;
else
inner=0;
end

if isempty(in)
outer=1;
end
if isempty(in1)
hold on;
if ~isempty(in)
plot(in(:,1),in(:,2),'b');
line(i,:)={[in(:,1);in(:,2)]};
hold off;
end
else
if vi(2)< in1(1,2)& isempty(in1)
ptx=[vi(1) vi(1)];
pty=[vi(1,2) in1(2,2)];

hold on;
plot(ptx,pty,'b')
line(i,:)={[ptx;pty]};
hold off;
elseif vi(2)> in1(2,2) & isempty(in1)
ptx=[vi(1) vi(1)];
pty=[vi(1,2) in1(2,2)];
hold on;
plot(ptx,pty,'b')
line(i,:)={[ptx;pty]};
hold on;
elseif vi(2)==in1(2,2)
ptx=[vi(1) vi(1)];
pty=[vi(1,2) in(2,2)];
hold on;
plot(ptx,pty,'b')
line(i,:)={[ptx;pty]};
hold on;
elseif ~isempty(in1)
if svertices(i-1,2)>vi(2)
ptx=[vi(1) vi(1)];
pty=[vi(1,2) in(1,2)];
hold on;
plot(ptx,pty,'b')
if pty(1)> pty(2)
swap=pty(1);
pty(1)=pty(2);
pty(2)=swap;
end
line(i,:)={[ptx;pty]};
hold on;
else
ptx=[vi(1) vi(1)];
pty=[vi(1,2) in(2,2)];
hold on;
plot(ptx,pty,'b')
if pty(1)> pty(2)
swap=pty(2);
pty(1)=pty(2);
pty(2)=swap;
end
line(i,:)={[ptx;pty]};
hold on;
end
else
ptx=[vi(1) vi(1)];
pty=[vi(1,2) in(1,2)];
hold on;
plot(ptx,pty,'b')
line(i,:)={[ptx;pty]};
hold on;



end


end













%line([vi(1) vi(1)], [-10 10]);

end
line(size(svertices,1),:)={[]};

if size(incell1,1)~=size(svertices,1)
for j=size(incell1,1)+1:size(svertices,1)
incell1(j,:)={[]};
end
end
k=1;
n=[];
polyStoreCell=[];
flag3=0;
flag4=0;
flag8=0;
flag10=0;
for i=1:size(svertices,1)-1

w =line{(i+1),:};
w = reshape(w',[],1);

if i<(size(svertices,1)-1)
z =line{i+2,:};
end
v =line{i,:};
n=incell1{i,:};
flag=0;
flag1=0;
txt=string(k);

u = svertices(i,:);
for j=1:size(ov,1)
if ov(j,1)==u(1)
ovpoint=ov(j,:);
flag1=1;
end
end
if isempty(v)
w=line{i+1,:};
flag=1;
pgon= polyshape([u(1) w(1) w(2)],[u(2) w(3) w(4)]);
polyStoreCell=[polyStoreCell;NaN NaN];
polyStoreCell=[polyStoreCell;pgon.Vertices];
hold on;
plot(pgon)
[q,r] = centroid(pgon);
text(q,r,txt)
hold off;

elseif i==size(svertices,1)-1
pgon1= polyshape([v(1) v(2) svertices(i+1,1)],[v(3) v(4) svertices(i+1,2)]);
polyStoreCell=[polyStoreCell;NaN NaN];
polyStoreCell=[polyStoreCell;pgon1.Vertices];
hold on;
plot(pgon1)
[q,r] = centroid(pgon1);
txt=string(k);
text(q,r,txt)
hold off;
break;

elseif ismember(u(1), ove(:,1))
pgon1= polyshape([v(1) v(2) w(1) w(2)],[v(3) v(4) w(4) w(3)]);
polyStoreCell=[polyStoreCell;NaN NaN];
polyStoreCell=[polyStoreCell;pgon1.Vertices];
hold on;
plot(pgon1)
[q,r] = centroid(pgon1);
txt=string(k);
text(q,r,txt)
hold off;
elseif flag1==1 && isempty(n) && u(1)~=max(ov(:,1))
if u(2)<svertices(i+1,2)
a=size(w);
a=a(1,1);
if a==2
small=min(w(2),w(4));
large=max(w(2),w(4));
xc=w(1);
x2c=w(3);
end
xc=w(1);
xc1=w(2);
small=w(3);
large=w(4);
pgon1= polyshape([u(1) v(2) w(1) w(2)],[u(2) v(4) large small]);
polyStoreCell=[polyStoreCell;NaN NaN];
polyStoreCell=[polyStoreCell;pgon1.Vertices];
hold on;
plot(pgon1)
[q,r] = centroid(pgon1);
text(q,r,txt)
hold off;
pgon2= polyshape([u(1) v(2) z(1) z(2)],[u(2) v(3) z(3) svertices(i+2,2)]);
polyStoreCell=[polyStoreCell;NaN NaN];
polyStoreCell=[polyStoreCell;pgon2.Vertices];
hold on;
plot(pgon2)
[q,r] = centroid(pgon2);
k=k+1;
txt=string(k);
text(q,r,txt)
hold off;
flag1=0;
flag8=1;
else
a=size(w);
a=a(1,1);
if a==2
small=min(w(2),w(4));
large=max(w(2),w(4));
xc=w(1);
x2c=w(3);
end
xc=w(1);
xc1=w(2);
small=min(w(3),w(4));
large=max(w(3),w(4));


pgon1= polyshape([u(1) v(2) w(1) w(2)],[u(2) v(3) small large]);
polyStoreCell=[polyStoreCell;NaN NaN];
polyStoreCell=[polyStoreCell;pgon1.Vertices];
hold on;
plot(pgon1)
[q,r] = centroid(pgon1);
text(q,r,txt)
hold off;
pgon2= polyshape([u(1) v(2) z(1) z(2)],[u(2) v(4) z(4) svertices(i+2,2)]);
polyStoreCell=[polyStoreCell;NaN NaN];
polyStoreCell=[polyStoreCell;pgon2.Vertices];
hold on;
plot(pgon2)
[q,r] = centroid(pgon2);
k=k+1;
txt=string(k);
text(q,r,txt)
hold off;
flag3=1;
flag1=0;
end


elseif ~isempty(n) && flag8==1


a=size(w);
a=a(1,1);
if a==2
small=min(w(2),w(4));
large=max(w(2),w(4));
xc=w(1);
x2c=w(3);
end
xc=w(1);
xc1=w(2);
small=w(3);
large=w(4)

small=min(v(2),v(4));
large=max(v(2),v(4));
pgon=polyshape([u(1) v(3) w(1) w(2)],[u(2) large w(4) svertices(i+1,2)]);
polyStoreCell=[polyStoreCell;NaN NaN];
polyStoreCell=[polyStoreCell;pgon.Vertices];
hold on;
plot(pgon)
[q,r] = centroid(pgon);
text(q,r,txt)
hold off;
flag10=1;
flag8=0;


elseif ~isempty(n) && flag3==1
small=min(v(2),v(4));
large=max(v(2),v(4));

pgon=polyshape([u(1) v(3) w(1) w(2)],[u(2) small w(3) svertices(i+1,2)]);
polyStoreCell=[polyStoreCell;NaN NaN];
polyStoreCell=[polyStoreCell;pgon.Vertices];
hold on;
plot(pgon)
[q,r] = centroid(pgon);
text(q,r,txt)
hold off;
flag4=1;
flag3=0;

elseif flag4==1 || flag10==1

flag4=0;
flag10=0;

pgon=polyshape([v(1) v(2) w(1) w(2)],[v(4) v(3) w(3) w(4)]);
polyStoreCell=[polyStoreCell;NaN NaN];
polyStoreCell=[polyStoreCell;pgon.Vertices];
hold on;
plot(pgon)
[q,r] = centroid(pgon);
text(q,r,txt)
hold off;
elseif i<size(svertices,1)-1
flag=0;
pgon=polyshape([v(1) v(2) w(1) w(2)],[v(4) v(3) w(3) w(4)]);
polyStoreCell=[polyStoreCell;NaN NaN];
polyStoreCell=[polyStoreCell;pgon.Vertices];
hold on;
plot(pgon)
[q,r] = centroid(pgon);
text(q,r,txt)
hold off;
else
flag=0;
small=min(v(2),v(4));
large=max(v(2),v(4));
pgon=polyshape([v(1) v(2) svertices(i+1,1)],[v(3) v(4) svertices(i+1,2)]);
polyStoreCell=[polyStoreCell;NaN NaN];
polyStoreCell=[polyStoreCell;pgon.Vertices];
hold on;
plot(pgon)
[q,r] = centroid(pgon);
text(q,r,txt)
hold off;
end
k=k+1;
end
toc

