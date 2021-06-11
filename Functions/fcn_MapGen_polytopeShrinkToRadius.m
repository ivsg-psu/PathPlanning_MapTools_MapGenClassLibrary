function [shrunk_polytope] = ...
    fcn_MapGen_polytopeShrinkToRadius(...
    shrinker,...
    new_radius,...
    tolerance,...
    varargin)
% fcn_MapGen_polytopesShrinkToRadius shrinks the polytopes to achieve the
% desired mean radius and specified variance
%
% FORMAT:
% 
% [shrunk_polytopes,mu_final,sigma_final] = ...
%     fcn_MapGen_polytopesShrinkToRadius(...
%     shrinker,...
%     new_radius,...
%     tolerance,...
%     (fig_num))
%
% INPUTS:
%
%     shrinker: original polytope with same fields as shrunk_polytopes
%     below
%
%     new_radius: desired polytope radius
%
%     tolerance: distance tolerance below which points of a polytope are
%     merged together
%     
%    (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%     SHRUNK_POLYTOPES: a 1-by-n seven field structure of shrunken polytopes, 
%     where n <= number of polytopes with fields:
%       vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%         the number of the individual polytope vertices
%       xv: a 1-by-m vector of vertice x-coordinates
%       yv: a 1-by-m vector of vertice y-coordinates
%       distances: a 1-by-m vector of perimeter distances from one point to the
%         next point, distances(i) = distance from vertices(i) to vertices(i+1)
%       mean: centroid xy coordinate of the polytope
%       area: area of the polytope
%       max_radius: distance from the mean to the farthest vertex
%   
% EXAMPLES:
%      
%
% For additional examples, see: script_test_fcn_MapGen_polytopeShrinkToRadius
%
% This function was written on 2019_08_29 by Seth Tau
% Questions or comments? sat5340@psu.edu 
%

% Revision History:
% 2021-06-08 - S. Brennan
% -- revised function to prep for MapGen class 
% -- added plotting option
% -- added comments, added debugging option

% TO DO
% -- Vectorize the for loop if possible
% -- check inputs are positive numbers where appropriate (e.g. make a
% "positive number" check

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    fig_for_debug = 9453;
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end

%% check input arguments?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if flag_check_inputs
    % Are there the right number of inputs?
    if nargin < 3 || nargin > 4
        error('Incorrect number of input arguments')
    end
    
    % Check the shrinker input
    fcn_MapGen_checkInputsToFunctions(...
        shrinker, 'polytopes');
    
    % Check the new_radius input
    fcn_MapGen_checkInputsToFunctions(...
        new_radius, 'column_of_numbers',1);
    
    % Check the tolerance input
    fcn_MapGen_checkInputsToFunctions(...
        tolerance, 'column_of_numbers',1);  
    
end
    

% Does user want to show the plots?
if  4== nargin
    fig_num = varargin{end};
    figure(fig_num);
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
        flag_do_plot = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the result:
shrunk_polytope = shrinker;

% pull values
vertices = shrinker.vertices;
centroid = shrinker.mean;
rad = shrinker.max_radius;

% determine scale factor
scale = new_radius/rad;

if scale < 1 % calculation error can sometimes make this greater than 1
    % find new vertices
    new_vert = centroid + scale*(vertices-centroid);
    
    % Sometimes, when shrinking, the new verticies are particularly
    % close to each other to where an edge "disappears". To prevent
    % this, we get rid of one of any vertices that are too close to
    % each other. This is set by a tolerance.
    tol_squared = tolerance.^2;
    
    % find all vertices that are larger away from next one, relative to
    % tolerance
    good_ind = ...
        sum((new_vert(1:end-1,:)-new_vert(2:end,:)).^2,2)>(tol_squared);
    
    % Check how many are good. If only 2 points are left, then the
    % result is just a line. If 1 or zero, then result is just a point.
    if sum(good_ind)>2 % sufficient good points to make a shape
        new_vert = new_vert(good_ind,:);
    elseif sum(good_ind)==2 % line shape
        new_vert = new_vert(good_ind,:);
        
        % The line may not go through the centroid, which is odd. We force
        % this by removing the point closest to the centroid
        distances_to_centroid = sum((new_vert-centroid).^2,2).^0.5;        
        if distances_to_centroid(1,1)>distances_to_centroid(2,1)
            new_vert(2,:)=centroid;
        else
            new_vert(1,:)=centroid;
        end
        
        new_vert = [new_vert; flipud(new_vert)];
    else % singular shape (i.e. point) or no shape
        new_vert = [centroid; centroid; centroid];
    end
    
    % adjust polytopes
    shrunk_polytope.vertices = [new_vert; new_vert(1,:)];
    shrunk_polytope.xv = new_vert(:,1)';
    shrunk_polytope.yv = new_vert(:,2)';
    shrunk_polytope.distances = INTERNAL_fcn_geometry_euclideanPointsToPointsDistance(new_vert(1:end-1,:),new_vert(2:end,:));
    shrunk_polytope.area = shrinker.area*scale^2;
    shrunk_polytope.max_radius = shrinker.max_radius*scale;
end
    


%% Plot results?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_do_plot
    figure(fig_num);
    hold on
       
    % Plot the input shrinker in red
    fcn_MapGen_plotPolytopes(shrinker,fig_num,'r',2);
    
    % plot the output polytope in blue
    fcn_MapGen_plotPolytopes(shrunk_polytope,fig_num,'b',2);

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _                 
%  |  ____|              | | (_)                
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___ 
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                               


function [dist] = ...
    INTERNAL_fcn_geometry_euclideanPointsToPointsDistance(...
    points1,...
    points2,...
    varargin)
% fcn_geometry_euclideanPointsToPointsDistance calculates the 
% distance(s) between a vector of points, POINTS1, and another vector of
% points, POINTS2.
%
% FORMAT:
%
% [DIST] = fcn_geometry_euclideanPointsToPointsDistance(POINTS1,POINTS2,(fig_num))
%
% INPUTS:
%
%      POINTS1: an Nx2 or Nx3 series of xy or xyz points 
%      in the form: [x1 y1 z1; x2 y2 z2; ... ; xn yn zn]
%
%      POINTS2: an Nx2 or Nx3 series of xy or xyz points 
%      in the form: [x1 y1 z1; x2 y2 z2; ... ; xn yn zn]
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      DIST: an N x  1 vector of distances [d1; d2; ... ; dn], where N is
%      the number of point sets
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%
% EXAMPLES:
%
%         pt1 = [1 1 5; 5 3 64; 7 2 -2];
%         pt2 = [0 -3 -6; 34 1 17; 18 7 0];
%         dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2);
%
% See the script: script_test_fcn_geometry_euclideanPointsToPointsDistance
% for a full test suite.
%
% This function was written on 2018_11_17 by Seth Tau
% Questions or comments? sat5340@psu.edu 

% Revision History:
% 2021-05-28 - S. Brennan
% -- revised function to prep for geometry class 
% -- rewrote function to use vector sum
% -- added plotting option
% 2021-06-05
% -- fixed comments, added debugging option


%% Debugging and Input checks
flag_check_inputs = 0; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end

%% check input arguments?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_check_inputs    
    % Are there the right number of inputs?    
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end
    
    % Check the points1 input
    fcn_geometry_checkInputsToFunctions(...
        points1, '2or3column_of_numbers');
end

% Use number of rows in points1 to calculate Npoints
Npoints = length(points1(:,1));

if flag_check_inputs
    
    % Check the points2 input, forcing length to match points1
    fcn_geometry_checkInputsToFunctions(...
        points2, '2or3column_of_numbers',Npoints);
    
end


% Does user want to show the plots?
if 3 == nargin
    fig_num = varargin{end};
    figure(fig_num);
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plot = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dist = sum((points1-points2).^2,2).^0.5;

%% Plot results?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_do_plot
    % Set up the figure
    figure(fig_num);
    clf
    hold on;
    grid on; grid minor;
        
    midpoints = (points1+points2)/2;
    for ith_point=1:Npoints
        % 2D plot?
        if length(midpoints(1,:))==2
            % Plot the points
            xdata = [points1(ith_point,1) points2(ith_point,1)];
            ydata = [points1(ith_point,2) points2(ith_point,2)];
            plot(xdata,ydata,'.-','Linewidth',3,'Markersize',20);
            
            % Label the midpoints
            text(midpoints(ith_point,1),midpoints(ith_point,2),sprintf('d - %.1f',dist(ith_point,1)));
        else
            % Plot the points
            xdata = [points1(ith_point,1) points2(ith_point,1)];
            ydata = [points1(ith_point,2) points2(ith_point,2)];
            zdata = [points1(ith_point,3) points2(ith_point,3)];
            plot3(xdata,ydata,zdata,'.-','Linewidth',3,'Markersize',20);

            % Label the midpoints
            text(midpoints(ith_point,1),midpoints(ith_point,2),midpoints(ith_point,3),sprintf('d - %.1f',dist(ith_point,1)));
            
            % Set to 3D view
            view(3);
        end
        
    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end


end % Ends the function

