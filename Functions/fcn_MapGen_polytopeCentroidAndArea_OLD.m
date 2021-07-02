function [Centroid,Area] = fcn_MapGen_polytopeCentroidAndArea(vertices)
% FCN_POLYTOPE_CALCULATION_CENTROID_AND_AREA calculates the centroid and
% area of a closed polytope.
%
% For details on the calculation, see:
% Centroid & Area equations: https://en.wikipedia.org/wiki/Centroid
%
% [CENTROID,AREA]=FCN_POLYTOPE_CALCULATION_CENTROID_AND_AREA(VERTICES)
% returns:
% CENTROID: the centroid given as [x-coordinate y_coordinate]
% AREA: the unsigned area enclosed by the polytope
%
% with inputs:
% vertices: [x y] where x and y are column vectors
% X: x coordinates of the polytope (with the same first and last point)
% Y: y coordinates of the polytope (with the same first and last point)
%
% The function outputs:
%
% Example:
% x = [3; 4; 2; -1; -2; -3; -4; -2; 1; 2; 3];
% y = [1; 2; 2; 3; 2; -1; -2; -3; -3; -2; 1];
% [Centroid,Area] = fcn_polytope_calculation_centroid_and_area([x,y])
% plot(x,y,'g-','linewidth',2)
% hold on
% plot(Centroid(:,1),Centroid(:,2),'kx','linewidth',1)
%
% This function was written on 2019_04_02 by Seth Tau, modifications by S.
% Brennan
% Questions or comments? sbrennan@psu.edu

% Revision History:
% 2021_02_23 by Seth Tau
% -- Added comments
% 2021_03_02 by Seth Tau
% -- Removed old add path stuff on 
% 2021_07_02
% -- Cleaned up arguments a bit to compactify x,y coordinate convention


% current points
xi = vertices(1:end-1,1); 
yi = vertices(1:end-1,2);

% next points
xip1 = vertices(2:end,1); 
yip1 = vertices(2:end,2);

% signed area
A = sum(xi.*yip1 - xip1.*yi)/2; 

% Centroid calculation
Cx = sum((xi+xip1).*(xi.*yip1 - xip1.*yi))/(6*A); % centroid x coordinate
Cy = sum((yi+yip1).*(xi.*yip1 - xip1.*yi))/(6*A); % centroid x coordinate
Centroid = [Cx, Cy];

Area = abs(A); % unsigned area

% Plotting
fig_num = 2222;
figure(fig_num)
clf;
hold on

plot(vertices(:,1),vertices(:,2),'b-','linewidth',2)

plot(Cx,Cy,'go','Markersize',10)