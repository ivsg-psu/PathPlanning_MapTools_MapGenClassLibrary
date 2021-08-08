function [cropped_vertices] = ...
    fcn_MapGen_polytopeRemoveColinearVertices(...,
    vertices,...
    varargin)

% Given a set of vertices that may have multiple points in a line on the
% same edge, this function returns an equivalent polytope with no repeats
% on any edge.
%
% FORMAT:
% 
% [cropped_vertices] = ...
%    fcn_MapGen_polytopeRemoveColinearVertices(...
%     vertices,...
%     (fig_num))
%
% INPUTS:
%
%     vertices: an Nx2 matrix of [x y] vertices
%
%    (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%     cropped_vertices: an Mx2 matrix of [x y] points, where no vertices
%     are repeated
%   
% DEPENDENCIES:
% 
%     fcn_MapGen_checkInputsToFunctions
% 
% EXAMPLES:
%      
% For additional examples, see:
% script_test_fcn_MapGen_polytopeRemoveColinearVertices
%
% This function was written on 2021_08_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu 
%

% Revision History:
% 2021_08_03 - S. Brennan
% -- first write of code

% TO DO
% -- none

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    fig_for_debug = 4564;
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
if  2 == nargin
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


% Remove repeats
[cropped_vertices,~,~] = unique(vertices,'rows','stable');

% Use the cross-product to eliminate co-linear points
Npoints = length(cropped_vertices(:,1));
good_indices = zeros(Npoints,1);

for ith_point = 1:Npoints
    
    if ith_point>1
        previous_vector = [cropped_vertices(ith_point,:)-cropped_vertices(ith_point-1,:), 0];
    else % must wrap around
        previous_vector = [cropped_vertices(end,:)-cropped_vertices(1,:), 0];
    end
    
    if ith_point<Npoints
        subsequent_vector = [cropped_vertices(ith_point+1,:)-cropped_vertices(ith_point,:), 0];
    else % must wrap around
        subsequent_vector = [cropped_vertices(1,:)-cropped_vertices(Npoints,:), 0];
    end
    
    cross_result = cross(previous_vector,subsequent_vector);
    
    if cross_result(1,3)~=0
        good_indices(ith_point,:) = 1;
    end
end

% Final polytope
cropped_vertices = cropped_vertices(good_indices>0,:);



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
    grid on
    grid minor
    hold on
    axis equal
    
    % Plot the input vertices
    plot(vertices(:,1),vertices(:,2),'r-');
    
    % Plot the cropped_vertices
    plot(cropped_vertices(:,1),cropped_vertices(:,2),'ko');    

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function


