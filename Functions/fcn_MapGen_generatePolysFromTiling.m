function [ ...
polytopes ...
] = ...
fcn_MapGen_generatePolysFromTiling( ...
seed_points, ...
V, ...
C, ...
AABB, ...
stretch, ...
varargin...
)
% fcn_MapGen_generatePolysFromTiling
% creates polytopes given seed points, V and C matrices from Voronoi 
% tiling, and stretch matrix
% 
% 
% 
% FORMAT:
% 
%    [ ...
%    polytopes ...
%    ] = ...
%    fcn_MapGen_generatePolysFromTiling( ...
%    seed_points, ...
%    V, ...
%    C, ...
%    AABB, ...
%    stretch, ...
%    (fig_num) ...
%    )
% 
% INPUTS:
% 
%     seed_points: the list of seed points in [x y] format, where x and y 
%     are columns
% 
%     V: the V matrix resulting from Voronoi calculations
% 
%     C: the C matrix resulting from Voronoi calculations
%
%     AABB: the axis-aligned bounding box, in format of 
%     [xmin ymin xmax ymax], wherein the resulting polytopes must be
%     bounded.
% 
%     stretch: the factor to stretch the polytopes after they are
%     calculated, in format of [xmagnification ymagnification]
% 
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a 
%     figure to be drawn showing results.
% 
% 
% OUTPUTS:
% 
%     polytopes: the resulting polytopes after converting to polytope form.
% 
% 
% DEPENDENCIES:
% 
%     fcn_MapGen_checkInputsToFunctions
% 
%     fcn_MapGen_cropPolytopeToRange
% 
%     fcn_MapGen_fillPolytopeFieldsFromVertices
% 
% 
% EXAMPLES:
% 
% See the script: script_test_fcn_MapGen_generatePolysFromTiling
% for a full test suite.
% 
% This function was written on 2021_07_02 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

% 
% REVISION HISTORY:
% 
% 2021_07_02 by Sean Brennan
% -- first write of function

% 
% TO DO:
% 
% -- fill in to-do items here.

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments 
flag_do_plot = 0;      % Set equal to 1 for plotting 
flag_do_debug = 0;     % Set equal to 1 for debugging 

if flag_do_debug
    fig_for_debug = 846;
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


if 1 == flag_check_inputs

    % Are there the right number of inputs?
    if nargin < 5 || nargin > 6
        error('Incorrect number of input arguments')
    end

    % Check the seed_points input, make sure it is '2column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        seed_points, '2column_of_numbers');
 
    % Check the stretch input, make sure it is '2column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        stretch, '2column_of_numbers',1);

end

% Does user want to show the plots?
if  6== nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
        flag_do_plot = 1;
    end
end

%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%See: http://patorjk.com/software/taag/#p=display&f=Big&t=Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§




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
    % Grab verticies    
    verticies = V(C{poly},:);    
    interior_point = seed_points(poly,:);

    %     if (interior_point(1,1)<0.0079/2) && (interior_point(1,1)>0.0078/2) && (interior_point(1,2)<0.754/2) && (interior_point(1,2)>0.753/2)
    %         disp('stop here');
    %     end
    
    %     % Are any vertices outside the [0,1] range?
    %     if any(xv>1) || any(yv>1) || any(xv<0) || any(yv<0)
    
    % Are any vertices outside the AABB?
    if ~all(fcn_MapGen_isWithinABBB(AABB,verticies)==1)
        
        % Crop vertices to allowable range        
        [cropped_vertices] = ...
            fcn_MapGen_cropPolytopeToRange(verticies, interior_point, AABB);
        %         xv = cropped_vertices(1:end-1,1)';
        %         yv = cropped_vertices(1:end-1,2)';
    else
        cropped_vertices = [verticies; verticies(1,:)]; % Close off the vertices
        %         % x and y values from this cell
        %         xv = verticies(:,1)';
        %         yv = verticies(:,2)';
    end 

    % Are polytopes not trivial in length? (This may not be needed)
    if length(cropped_vertices(:,1))>2                
    
        % make sure cw
        vec1 = [cropped_vertices(2,:)-cropped_vertices(1,:),0]; % vector leading into point
        vec2 = [cropped_vertices(3,:)-cropped_vertices(2,:),0]; % vector leading out of point
        xing = cross(vec1,vec2); % cross product of two vectors
        if sign(xing(3)) == -1 % points ordered in wrong direction
            cropped_vertices = flipud(cropped_vertices);
        end
        
        % enter info into polytope structure
        polytopes(poly-remove).vertices = cropped_vertices; 

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
end % Ends for loop for stretch

% Fill in all the other fields
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes);
    
%§
%% Plot the results (for debugging)?
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
    
    % plot the seed points in red
    plot(seed_points(:,1),seed_points(:,2),'r.','Markersize',10);
    
    % plot the means in black
    temp = zeros(length(polytopes),2);
    for ith_poly = 1:length(polytopes)
        temp(ith_poly,:) = polytopes(ith_poly).mean;
    end
    plot(temp(:,1),temp(:,2),'ko','Markersize',3);
    
end % Ends the flag_do_plot if statement    

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends the function

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _                 
%  |  ____|              | | (_)                
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___ 
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%                                               
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§




