function [shrunk_polytopes,mu_final,sigma_final] = ...
    fcn_MapGen_polytopeShrinkToRadius(...
    polytopes,...
    des_radius,...
    sigma_radius,...
    min_rad,...
    varargin)
% fcn_MapGen_polytopeShrinkToRadius shrinks the polytopes to achieve the
% desired mean radius and specified variance
%
% FORMAT:
% 
% [shrunk_polytopes,mu_final,sigma_final] = ...
%     fcn_MapGen_polytopeShrinkToRadius(...
%     polytopes,...
%     des_radius,...
%     sigma_radius,...
%     min_rad,...
%     (fig_num))
%
% INPUTS:
%
%     POLYTOPES: original polytopes with same fields as shrunk_polytopes
%
%     DES_RAD: desired average max radius   
%
%     SIGMA_RADIUS: desired variance in the radii 
%
%     MIN_RAD: minimum acceptable radius
%
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
%     MU_FINAL: final average maximum radius achieved
%
%     SIGMA_FINAL: final variance achieved
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
    if nargin < 4 || nargin > 5
        error('Incorrect number of input arguments')
    end
    
    % Check the polytopes input
    fcn_MapGen_checkInputsToFunctions(...
        polytopes, 'polytopes');
    
    % Check the des_radius input
    fcn_MapGen_checkInputsToFunctions(...
        des_radius, 'column_of_numbers',1);
    
    % Check the sigma_radius input
    fcn_MapGen_checkInputsToFunctions(...
        sigma_radius, 'column_of_numbers',1);
 
    % Check the min_rad input
    fcn_MapGen_checkInputsToFunctions(...
        min_rad, 'column_of_numbers',1);
    
    
end
    

% Does user want to show the plots?
if  5== nargin
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


%% find current distribution
radii = [polytopes.max_radius];

if flag_do_debug
    fcn_MapGen_plotPolytopes(polytopes,fig_for_debug,'b',2);
    
    figure(fig_for_debug+1);
    histogram(radii,20)
    title('Histogram of input radii');
    r_sigma = std(radii);
    fprintf(1,'Standard deviation in r is: %.2f \n', r_sigma);
end

r_mu = mean(radii);
r_size = length(radii);

if r_mu < des_radius
    error('cannot achieve the desired radius by shrinking because average radius is already smaller than desired radius')
end

%% determine desired distribution
r_dist = normrnd(des_radius,sigma_radius,[r_size,1]);

if flag_do_debug
   figure(fig_for_debug+2);  
   histogram(r_dist,20)
end

r_dist = r_dist + (des_radius-mean(r_dist)); % adjust to ensure the mean value is mu
max_r_dist = max(radii); % largest possible radius
min_r_dist = min_rad; % smallest possible radius
if sum((r_dist>max_r_dist)+(r_dist<min_r_dist)) > 0
    warning('standard deviation skewed due to truncated distribution')
end
r_dist(r_dist>max_r_dist) = max_r_dist; % truncate any values that are too large
r_dist(r_dist<min_r_dist) = min_r_dist; % truncate any values that are too small
while abs(mean(r_dist) - des_radius) > 1e-10
    r_dist = r_dist + (des_radius-mean(r_dist));
    r_dist(r_dist>max_r_dist) = max_r_dist;
    r_dist(r_dist<min_r_dist) = min_r_dist;
end

mu_final = mean(r_dist);
sigma_final = std(r_dist);


if flag_do_debug
   figure(fig_for_debug+2);  
   histogram(r_dist,20)
end

%% shrink polytopes to achieve the distribution
[new_rads,ob_ind] = sort(r_dist);
if sum((sort(radii)'-sort(r_dist))>=-2*min_rad) < r_size
    error('distribution is unachievable with generated map')
end

shrunk_polytopes = polytopes;
for idx = 1:length(new_rads)
    shrinker = polytopes(ob_ind(idx)); % obstacle to be shrunk
    
    % pull values
    vertices = shrinker.vertices;
    centroid = shrinker.mean;
    rad = shrinker.max_radius;
    
    % determine scale factor
    scale = new_rads(idx)/rad;
    
    if scale < 1 % calculation error can sometimes make this greater than 1
        % find new vertices
        new_vert = centroid + scale*(vertices-centroid);

        % adjust polytopes
        shrinker.vertices = new_vert;
        shrinker.xv = new_vert(1:end-1,1)';
        shrinker.yv = new_vert(1:end-1,2)';
        shrinker.distances = INTERNAL_fcn_geometry_euclideanPointsToPointsDistance(new_vert(1:end-1,:),new_vert(2:end,:));
        shrinker.area = shrinker.area*scale^2;
        shrinker.max_radius = shrinker.max_radius*scale;
    end
    
    % assign to shrunk_polytopes
    shrunk_polytopes(ob_ind(idx)) = shrinker;
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
    
    % Plot the input polytopes in red
    fcn_MapGen_plotPolytopes(polytopes,fig_num,'r',2,[0 1 0 1]);
    
    % plot the shrunk in blue
    fcn_MapGen_plotPolytopes(shrunk_polytopes,fig_num,'b',2,[0 1 0 1]);

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

