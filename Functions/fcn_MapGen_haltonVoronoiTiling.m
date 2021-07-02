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
% 2021_07_02
% -- fixed boundary edge case
% -- buffer the Halton set to be sure polytopes are well formed on edges
% -- force correct number of polytopes!

% TO DO:
% -- check cross product around entire polytope
% -- add unit normal vectors for each edge

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
if 2 <= nargin
    stretch = varargin{1};
    
    % Check the stretch input
    fcn_MapGen_checkInputsToFunctions(...
        stretch, '2column_of_numbers',1);
       
end


% Does user want to show the plots?
if 3 == nargin
    fig_num = varargin{end};
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
[V,C] = voronoin(seed_points);
% V = V.*stretch;

%% fill polytopes from tiling
polytopes = INTERNAL_fcn_MapGen_generatePolysFromTiling(seed_points,V,C,stretch);


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
    fcn_MapGen_plotPolytopes(polytopes,fig_num,'b',2);

    % plot the seed points
    plot(seed_points(:,1),seed_points(:,2),'r.','Markersize',10);
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end


end % Ends the function

%% create tiling
function polytopes = INTERNAL_fcn_MapGen_generatePolysFromTiling(seed_points,V,C,stretch)

flag_do_debug = 1;

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
    
    % Are any vertices outside the [0,1] range?
    if any(xv>1) || any(yv>1) || any(xv<0) || any(yv<0)

        % Crop vertices to allowable range
        [cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies, interior_point);
        xv = cropped_vertices(1:end-1,1)';
        yv = cropped_vertices(1:end-1,2)'; 
    else
        
    end 

    % Are polytopes not trivial in length? (This may not be needed)
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
        polytopes(poly-remove).distances = sum((polytopes(poly-remove).vertices(1:end-1,:)-polytopes(poly-remove).vertices(2:end,:)).^2,2).^0.5;
        
        % calculate the maximum distance from center to a vertex
        polytopes(poly-remove).max_radius = max(sum((polytopes(poly-remove).vertices(1:end-1,:)-ones(length(xv),1)*polytopes(poly-remove).mean).^2,2).^0.5);


    else % if 2 or less points in cell 
        remove = remove+1; % skip cell and remove later
    end
end

% remove extra empty polytopes
polytopes = polytopes(1:(num_poly-remove));

% Apply the stretch
num_poly = length(polytopes);
for poly = 1:num_poly % pull each cell from the voronoi diagram
    
    polytopes(poly).vertices  = polytopes(poly).vertices.*stretch;
    polytopes(poly).xv        = (polytopes(poly).vertices(1:end-1,1)').*stretch(1,1);
    polytopes(poly).yv        = (polytopes(poly).vertices(1:end-1,2)').*stretch(1,2);
    polytopes(poly).distances = sum((polytopes(poly).vertices(1:end-1,:)-polytopes(poly).vertices(2:end,:)).^2,2).^0.5;
    [Cx,Cy,polytopes(poly).area] = ...
        fcn_MapGen_polytopeCentroidAndArea([xv xv(1)]',[yv yv(1)]');
    
    polytopes(poly).mean = [Cx, Cy]; % enter polytope centroid
    polytopes(poly).max_radius = max(sum((polytopes(poly).vertices(1:end-1,:)-ones(length(polytopes(poly).xv),1)*polytopes(poly).mean).^2,2).^0.5);
end % Ends for loop for stretch

end


