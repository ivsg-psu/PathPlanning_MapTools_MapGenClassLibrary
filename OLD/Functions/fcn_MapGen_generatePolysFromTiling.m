function [ ...
polytopes ...
] = ...
fcn_MapGen_generatePolysFromTiling( ...
seed_points, ...
V, ...
C, ...
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
%     stretch: the list of seed points in [x y] format, where x and y are 
%     columns
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
%     fcn_MapGen_fillPolytopeFieldsFromVerticies
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
    fig_for_debug = 404;
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
    if nargin < 4 || nargin > 5
        error('Incorrect number of input arguments')
    end

    % Check the seed_points input, make sure it is '2column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        seed_points, '2column_of_numbers');
 
    % Check the stretch input, make sure it is '2column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        stretch, '2column_of_numbers',[1 1 ]);
 
end

% Does user want to show the plots?
if  5== nargin
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
        polytopes(poly-remove).vertices = [[xv xv(1)]' [yv yv(1)]']; % repeat first vertice for easy plotting
        
        polytopes(poly-remove) = fcn_MapGen_fillPolytopeFieldsFromVerticies(polytopes(poly-remove));

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
    polytopes(poly) = fcn_MapGen_fillPolytopeFieldsFromVerticies(polytopes(poly));
    
end % Ends for loop for stretch

% 

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
    % Nothing to plot here
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

end



