function [angles, unit_in_vectors, unit_out_vectors] = ...
    fcn_MapGen_polytopeFindVertexAngles(...
    vertices,...
    varargin)
% fcn_MapGen_polytopeFindVertexAngles finds the angles, in radians, of each
% vertex relative to each other within a polytope. Note that the 1st angle
% corresponds to the 1st vertex, which is located at the 1nd point in the
% vertex list.
%
% This operation is done by doing the cross and dot products. The cross
% product determines the sign of the angle (+ or -) whereas the dot product
% determines the value. This, the angle returned is the OUTSIDE angle
% between the incoming line and the outgoing line. For the INSIDE angle,
% this is given by INSIDE = pi - OUTSIDE.
%
% FORMAT:
% 
% [angles, unit_in_vectors, unit_out_vectors] = ...
%     fcn_MapGen_polytopeFindVertexAngles(...
%     vertices,...
%     (fig_num))
%
% INPUTS:
%
%     vertices: N x 2 vector of vertices, in [x y] format
%   
%    (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%     angles: an Nx1 column of angles, one for each vertex
%  
%     unit_in_vectors: the Nx2 matrix of unit vectors leading into each vertex
%  
%     unit_out_vectors: the Nx2 matrix of unit vectors leading out of each vertex
%   
% DEPENDENCIES:
% 
%     fcn_MapGen_checkInputsToFunctions
% 
% EXAMPLES:
%
% For additional examples, see: script_test_fcn_MapGen_polytopeFindVertexAngles
%
% This function was written on 2021_08_01 by S. Brennan
% Questions or comments? sbrennan@psu.edu 
%

% Revision History:
% 2021_08_01 - S. Brennan
% -- first write of the code

% TO DO
% -- (none)

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
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end
    
    % Check the vertices input
    fcn_MapGen_checkInputsToFunctions(...
        vertices, '2column_of_numbers');
    
end
    

% Does user want to show the plots?
if  2== nargin
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

% Start by defining the vectors - this is done by using sequences of
% three points to define two vectors.
Nvertices = length(vertices(:,1));
Nangles = Nvertices-1;

% Create a set of calculation vertices which are all the verticies,
% plus 2 wrap-arounds (1 is already there) so the difference
% calculations for start/end vectors are easy to represent
calculation_vertices = [vertices(end-1,:); vertices];

% Determine starting and ending vectors. Start is created from first
% points, the end is created from last 2 points
in_vectors = calculation_vertices(2:end-1,:)- calculation_vertices(1:end-2,:);
out_vectors = calculation_vertices(3:end,:)- calculation_vertices(2:end-1,:);

% Convert both to unit vectors
unit_in_vectors = in_vectors./(sum(in_vectors.^2,2).^0.5);
unit_out_vectors   = out_vectors./(sum(out_vectors.^2,2).^0.5);

% Do the cross product, and pull out only z-component result
% cross_products = cross([unit_starting_vector zeros(Nangles,1)],[unit_ending_vector zeros(Nangles,1)]);
% cross_result = cross_products(:,3);
cross_result = crossProduct(unit_in_vectors,unit_out_vectors);


% Grab the angle via cross product
angles_cross = asin(cross_result);

% Do the dot product, and grab angle via dot product
dot_products = dot(unit_in_vectors,unit_out_vectors,2);
angles_dot = acos(dot_products);

% Calculate the angles
angles(1:Nangles,1)  = angles_dot.*sign(angles_cross);

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
    clf;
    grid on;
    hold on
    axis equal;
    
    % Plot the polytope in red
    plot(vertices(:,1),vertices(:,2),'r-','Linewidth',2);

    % Find size of vertices
    size = max(max(vertices)) - min(min(vertices));
    nudge = size*0.01;
    
    % Label the vertices
    for ith_angle = 1:Nangles
        ith_vertex = ith_angle;
        text(vertices(ith_vertex,1)+nudge,vertices(ith_vertex,2),...
            sprintf('%.0f deg',180-angles(ith_angle,1)*180/pi));
    end

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


%% Calculate cross products
function result = crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end
