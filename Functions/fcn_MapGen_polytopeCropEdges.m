function [trim_polytopes] = ...
    fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,varargin)
% fcn_MapGen_polytopeCropEdges removes polytopes that extend 
% beyond the boundaries specified
%
% FORMAT:
% 
% [trim_polytopes] = ...
%     fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,(fig_num))
%
% INPUTS:
%
%     polytopes: the original polytopes with the same fields as trim_polytopes
% 
%     bounding_box: a 2 x 2 matrix of [xlow ylow; xhigh yhigh] in which all
%     the polytopes must exist, e.g. the corner coordinates of the
%     axis-aligned bounding box.
%     
%    (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results.
%
%
% OUTPUTS:
%
%     TRIM_POLYTOPES: a 1-by-n seven field structure of polytopes within the 
%     boundaries, where n <= number of polytopes with fields:
%     vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%     the number of the individual polytope vertices
%     xv: a 1-by-m vector of vertice x-coordinates
%     yv: a 1-by-m vector of vertice y-coordinates
%     distances: a 1-by-m vector of perimeter distances from one point to the
%     next point, distances(i) = distance from vertices(i) to vertices(i+1)
%     mean: centroid xy coordinate of the polytope
%     area: area of the polytope
%
%   
% EXAMPLES:
%      
%   polytopes = fcn_MapGen_haltonVoronoiTiling([1 1000]);
%   fig = fcn_plot_polytopes(polytopes,[],'b',2,[0 1 0 1]);
%
%   bounding_box = [0,1; 0,1];
%   trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box);
%   fcn_plot_polytopes(trim_polytopes,fig,'g',2,[0 1 0 1]);
%
% For additional examples, see: script_test_fcn_MapGen_polytopeCropEdges
%
% This function was written on 2019_06_13 by Seth Tau
% Questions or comments? sat5340@psu.edu 
%

% REVISIONS
% 2021-06-06
% -- rewrote function from fc_polytope_editing_remove_edge_polytopes
% -- added a test script
% 2021-06-10
% -- updated comments

% TO DO
% -- Vectorize the for loop if possible

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
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end
    
    % Check the polytopes input
    fcn_MapGen_checkInputsToFunctions(...
        polytopes, 'polytopes');
    
    % Check the bounding_box input
    fcn_MapGen_checkInputsToFunctions(...
        bounding_box, '2column_of_numbers',2);

end
    

% Does user want to show the plots?
if  3== nargin
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

xlow = bounding_box(1,1);
ylow = bounding_box(1,2);
xhigh = bounding_box(2,1);
yhigh = bounding_box(2,2);

Npolys = length(polytopes);

% Preallocate the polytopes
trim_polytopes(Npolys) = struct('vertices',[],'xv',[],'yv',[],'distances',[],'mean',[],'area',[],'max_radius',[]);

keep = 0; % number of polytopes to keep
for poly = 1:Npolys % check each polytope within polytopes
    xv = polytopes(poly).xv;
    yv = polytopes(poly).yv;
    if sum((xv<xlow)+(xv>xhigh)+(yv<ylow)+(yv>yhigh))==0 % if the x or y vertices are inside of the bounds
        keep = keep + 1;
        trim_polytopes(keep) = polytopes(poly);
    end
end
trim_polytopes = trim_polytopes(1:keep); % Save only the ones that were filled

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
    fcn_MapGen_plotPolytopes(polytopes,fig_num,'r',2,[xlow xhigh ylow yhigh]);
    
    % plot the tiled_polytopes in blue
    fcn_MapGen_plotPolytopes(trim_polytopes,fig_num,'b',2,[xlow xhigh ylow yhigh]);

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end


end % Ends the function



