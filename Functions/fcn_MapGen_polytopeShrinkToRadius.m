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

% calculation error can sometimes make the scale greater than 1, so if we
% are doing shrinking, check that actual shrinking is happening!
if scale < 1 
    % find new vertices
    new_vert = centroid + scale*(vertices-centroid);
    shrunk_polytope.vertices = [new_vert; new_vert(1,:)];
end

% adjust polytopes
shrunk_polytope = fcn_MapGen_fillPolytopeFieldsFromVerticies(shrunk_polytope);


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
    
    % Plot the cetroid in black
    plot(centroid(:,1),centroid(:,2),'ko','Markersize',10);
    
    % Plot the input shrinker in red
    fcn_MapGen_plotPolytopes(shrinker,fig_num,'r',2);
    
    % plot the output polytope in blue
    fcn_MapGen_plotPolytopes(shrunk_polytope,fig_num,'b',2);

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function




