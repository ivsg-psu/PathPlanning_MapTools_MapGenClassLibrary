function [ ...
exp_polytopes ...
] = ...
fcn_MapGen_polytopesExpandEvenlyForConcave( ...
polytopes, ...
exp_dist, ...
varargin...
)
% fcn_MapGen_polytopesExpandEvenlyForConcave
% Expands an obstacle out by exp_dist on all sides using matlab's polyshape object
% and the polybuffer object function (method). The utility of this is that this
% method is compatible with nonconvex polytopes while the implementation of
% fcn_MapGen_polytopesExpandEvenly is not.
%
%
%
% FORMAT:
%
%    [ ...
%    exp_polytopes ...
%    ] = ...
%    fcn_MapGen_polytopesExpandEvenlyForConcave( ...
%    polytopes, ...
%    delta, ...
%    exp_dist, ...
%    (fig_num) ...
%    )
%
% INPUTS:
%
%     polytopes: the structure of 'polytopes' type that stores the
%     polytopes to be expanded
%
%     exp_dist: distance to expand the obstacle
%
%     (optional inputs)
%
%     fig_num: any number that acts somewhat like a figure number output.
%     If given, this forces the variable types to be displayed as output
%     and as well makes the input check process verbose.
%
%
% OUTPUTS:
%
%     exp_polytopes: structure of expanded polytopes
%
%
% DEPENDENCIES:
%
%     fcn_MapGen_checkInputsToFunctions
%     fcn_MapGen_fillPolytopeFieldsFromVertices
%     fcn_MapGen_plotPolytopes
%     MATLAB's polyshape object and polybuffer object function (method)
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_polytopesExpandEvenlyForConcave
% for a full test suite.
%
% This function was written 5 Feb. 2024 by Steve Harnett
% Questions or comments? contact sjharnett@psu.edu

%
% REVISION HISTORY:
%
% 2024_02_05, Steve Harnett
% -- first write of script

%
% TO DO:
%
% -- fill in to-do items here.

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    fig_for_debug = 681;
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
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end

    % Check the polytopes input, make sure it is 'polytopes' type
    fcn_MapGen_checkInputsToFunctions(...
        polytopes, 'polytopes');


    % Check the exp_dist input, make sure it is 'positive_column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        exp_dist, 'positive_1column_of_numbers',1);

end

% Does user want to show the plots?
if  3== nargin
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

clear exp_polytopes;
for p = 1:length(polytopes)
    this_polytope = polytopes(p); % look at one polytope
    this_polyshape = polyshape(this_polytope.vertices); % convert it to matlab polyshape
    scaled_polyshape = polybuffer(this_polyshape,exp_dist,'JointType','miter','MiterLimit',2); % use polyshape to enlarge it by a buffer
    new_vertices = scaled_polyshape.Vertices; % extract the vertices from the polyshape
    new_vertices = [new_vertices; new_vertices(1,:)]; % duplicate first vertex at end of array
    exp_polytopes(p).vertices = new_vertices; % store vertices in expanded poly struct array
end
exp_polytopes= fcn_MapGen_fillPolytopeFieldsFromVertices(exp_polytopes,1009,1); % fill polytopes from vertices

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
    figure(fig_num)
    clf;

    LineWidth = 2;
    fcn_MapGen_plotPolytopes(polytopes,fig_num,'r-',LineWidth);
    fcn_MapGen_plotPolytopes(exp_polytopes,fig_num,'b-',LineWidth,'square');
    legend('Original','Expanded')
    box on
    xlabel('X Position')
    ylabel('Y Position')

end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends the function
