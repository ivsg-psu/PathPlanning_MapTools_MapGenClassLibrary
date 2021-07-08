function [polytopes] = ...
    fcn_MapGen_latinVoronoiTiling(latin_range,varargin)
% fcn_MapGen_latinVoronoiTiling generates a map with
% obstacles perfectly tiled together using the Voronoi cells generated from
% the Latin Hypercube sequence. See more about this at:
% https://www.mathworks.com/help/stats/generating-quasi-random-numbers.html
%
% FORMAT:
% 
% [polytopes] = ...
%    fcn_MapGen_latinVoronoiTiling(latin_range,(stretch),(fig_num))
%
% INPUTS:
%
%    latin_range: 1 x 2 vector of integers specifying the [low high] range
%    of latin point indices to use to generate the tiling, where  
%    low: the lowest point index to use in the latin sequence
%    high: the highest point index to use in the latin sequence. NOTE: for
%    latin sequences, only the range of numbers is used, e.g. high minus
%    low. So an input of [1 100] is the same as [2 101].
%     
%    (OPTIONAL INPUTS)
%
%    stretch: [x,y] scaling factors to allow the values from the latin set
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
% See the script: script_test_fcn_MapGen_latinVoronoiTiling
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
% -- buffer the latin set to be sure polytopes are well formed on edges
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
    
    % Check the latin_range input
    fcn_MapGen_checkInputsToFunctions(...
        latin_range, '2column_of_integers');

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
%% pull latin set
Npoints = max(latin_range) - min(latin_range) + 1;
latin_points = lhsdesign(Npoints,2);

%% pick values from latin set
seed_points = latin_points;
[V,C] = voronoin(seed_points);
% V = V.*stretch;

%% fill polytopes from tiling
polytopes = fcn_MapGen_generatePolysFromTiling(seed_points,V,C,stretch);



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















