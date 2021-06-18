function [polytopes] = ...
    fcn_MapGen_haltonVoronoiTiling(Halton_range,varargin)
% fcn_MapGen_haltonVoronoiTiling generates a map with
% obstacles perfectly tiled together using the Voronoi cells generated from
% the Halton sequence
%
% FORMAT:
% 
% [polytopes] = ...
%    fcn_MapGen_haltonVoronoiTiling(Halton_range,(stretch),(fig_num))
%
% INPUTS:
%
%    Halton_range: 1 x 2 vector of integers specifying the [low high] range
%    of Halton point indices to use to generate the tiling, where  
%    low: the lowest point index to use in the halton sequence
%    high: the highest point index to use in the halton sequence
%     
%    (OPTIONAL INPUTS)
%
%    stretch: [x,y] scaling factors to allow the values from the halton set
%    to be scaled to fit maps shapes other than square
%
%    fig_num: a figure number to plot results.
%
% OUTPUTS:
%
% POLYTOPES: a 1-by-n seven field structure, where n <= number of polytopes
%   with fields:
%   vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%     the number of the individual polytope vertices
%   xv: a 1-by-m vector of vertice x-coordinates
%   yv: a 1-by-m vector of vertice y-coordinates
%   distances: a 1-by-m vector of perimeter distances from one point to the
%     next point, distances(i) = distance from vertices(i) to vertices(i+1)
%   mean: centroid xy coordinate of the polytope
%   area: area of the polytope
%
% DEPENDENCIES:
%
%      fcn_MapGen_checkInputsToFunctions
%      fcn_MapGen_polytopeCentroidAndArea
%      fcn_MapGen_plotPolytopes
%
% EXAMPLES:
%      
% See the script: script_test_fcn_MapGen_haltonVoronoiTiling
% for a full test suite.
%
% This function was written on 2019_06_13 by Seth Tau
% Questions or comments? sbrennan@psu.edu 

% REVISION HISTORY:
% 2021_06_06 
% -- edited by S. Brennan to put it into MapGen format

% TO DO:
% -- check cross product around entire polytope
% -- add unit normal vectors for each edge
% -- buffer the Halton set to be sure polytopes are well formed on edges
% -- force correct number of polytopes?

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    fig_for_debug = 9993;
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
    if nargin < 1 || nargin > 3
        error('Incorrect number of input arguments')
    end
    
    % Check the Halton_range input
    fcn_MapGen_checkInputsToFunctions(...
        Halton_range, '2column_of_integers');

end


% check variable argument
stretch = [1 1]; % default stretch value
if 2 == nargin
    stretch = varargin{1};
    
    % Check the stretch input
    fcn_MapGen_checkInputsToFunctions(...
        stretch, '2column_of_numbers',1);
       
end


% Does user want to show the plots?
if 3 == nargin
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
%% pull halton set
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values

%% pick values from halton set
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
seed_points = points_scrambled(low_pt:high_pt,:);
seed_points = seed_points.*stretch;
[V,C] = voronoin(seed_points);
% V = V.*stretch;

%% fill polytopes from tiling
polytopes = INTERNAL_fcn_MapGen_createPolysFromTiling(seed_points,V,C);


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
 
    % plot the polytopes
    fcn_MapGen_plotPolytopes(polytopes,fig_num,'b',2,[0 1 0 1]);

    % plot the seed points
    plot(seed_points(:,1),seed_points(:,2),'r.','Markersize',10);
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end


end % Ends the function


function polytopes = INTERNAL_fcn_MapGen_createPolysFromTiling(seed_points,V,C)
%% create tiling
num_poly = size(seed_points,1);
polytopes(num_poly) = ...
    struct(...
    'vertices',[],...
    'xv',[],...
    'yv',[],...
    'distances',[],...
    'mean',[],...
    'area',[],...
    'max_radius',[]);


remove = 0; % keep track of how many cells to be removed
for poly = 1:num_poly % pull each cell from the voronoi diagram
    % x and y values from this cell
    xv = V(C{poly},1)'; 
    yv = V(C{poly},2)';
    
    verticies = V(C{poly},:);
    interior_point = seed_points(poly,:);
    
    %     if interior_point(1,1)>0.97851 && interior_point(1,1)<0.97852
    %         disp('Stop here')
    %     end
    
    % FOR DEBUGGING
    %         plotting_vertices = [vertices; vertices(1,:)];
    %         figure(999);
    %         hold on;
    %         plot(interior_point(1,1),interior_point(1,2),'r.');
    %         plot(plotting_vertices(:,1),plotting_vertices(:,2),'g-');

    % Are any vertices outside the range?
    if any(xv>1) || any(yv>1) || any(xv<0) || any(yv<0)
        
        % Are any verticies infinite? If so, we need to replace the
        % infinite terms with real data. A way to do this is to fill the
        % infinite point with "snapped" points created by projecting the
        % adjacent vertices to the infinite one toward the nearest wall.
        if any(isinf(verticies),'all')
            Npoints = length(xv);            
            bad_indices = find(any(isinf(verticies),2));
            if length(bad_indices)>1
                warning('More than 2 infinities found in one vector. Code may break');
            end
            for ith_index = 1:length(bad_indices)
                bad_index = bad_indices(ith_index);

                % Grab prior and next indices, being careful to check for
                % situations where the infinite index is at start or end.
                prior_index = bad_index-1;
                next_index  = bad_index+1;
                if bad_index == 1
                    prior_index = Npoints;
                    start_data = [];
                else
                    start_data = verticies(1:prior_index,:);
                end
                if bad_index == Npoints
                    next_index = 1;
                    end_data = [];
                else
                    end_data = verticies(next_index:end,:);                    
                end
                
                % snap prior to closest wall
                snap_prior = round(verticies(prior_index,:));
                differences = (verticies(prior_index,:)-snap_prior).^2;
                if differences(1,1)<differences(1,2)
                    new_prior = [snap_prior(1,1) verticies(prior_index,2)];
                else
                    new_prior = [verticies(prior_index,1) snap_prior(1,2)];
                end
                
                % snap next to closest wall
                snap_next = round(verticies(next_index,:));
                differences = (verticies(next_index,:)-snap_next).^2;
                if differences(1,1)<differences(1,2)
                    new_next = [snap_next(1,1) verticies(next_index,2)];
                else
                    new_next = [verticies(next_index,1) snap_next(1,2)];
                end
                
                % Substitute data in
                verticies = [start_data; new_prior; new_next; end_data];                
  
            end
        end
            
        %verticies = verticies(flag_is_finite>0,:);

        [cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies, interior_point);
        xv = cropped_vertices(:,1)';
        yv = cropped_vertices(:,2)';        
    end

    % Are polytopes not trivial in length?
    if length(xv)>2                
    
        % make sure cw
        vec1 = [xv(2)-xv(1),yv(2)-yv(1),0]; % vector leading into point
        vec2 = [xv(3)-xv(2),yv(3)-yv(2),0]; % vector leading out of point
        xing = cross(vec1,vec2); % cross product of two vectors
        if sign(xing(3)) == -1 % points ordered in wrong direction
            xv = fliplr(xv);
            yv = fliplr(yv);
        end
        
        % enter info into polytope structure
        polytopes(poly-remove).xv = xv;
        polytopes(poly-remove).yv = yv;
        polytopes(poly-remove).vertices = [[xv xv(1)]' [yv yv(1)]']; % repeat first vertice for easy plotting
        
        [Cx,Cy,polytopes(poly-remove).area] = ...
            fcn_MapGen_polytopeCentroidAndArea([xv xv(1)]',[yv yv(1)]');
        
        polytopes(poly-remove).mean = [Cx, Cy]; % enter polytope centroid
        % calculate perimeter distances around the polytope
        polytopes(poly-remove).distances = INTERNAL_fcn_geometry_euclideanPointsToPointsDistance(polytopes(poly-remove).vertices(1:end-1,:),polytopes(poly-remove).vertices(2:end,:));
        % calculate the maximum distance from center to a vertex
        polytopes(poly-remove).max_radius = max(INTERNAL_fcn_geometry_euclideanPointsToPointsDistance(polytopes(poly-remove).vertices(1:end-1,:),ones(length(xv),1)*polytopes(poly-remove).mean));


    else % if 2 or less points in cell 
        remove = remove+1; % skip cell and remove later
    end
end

%% remove extra empty polytopes
polytopes = polytopes(1:num_poly);
end


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
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
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

% if flag_check_inputs    
%     % Are there the right number of inputs?    
%     if nargin < 2 || nargin > 3
%         error('Incorrect number of input arguments')
%     end
%     
%     % Check the points1 input
%     fcn_geometry_checkInputsToFunctions(...
%         points1, '2or3column_of_numbers');
%     
%     % Use number of rows in points1 to calculate Npoints
%     Npoints = length(points1(:,1));
%     
%     % Check the points2 input, forcing length to match points1
%     fcn_geometry_checkInputsToFunctions(...
%         points2, '2or3column_of_numbers',Npoints);       
% end
    

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













