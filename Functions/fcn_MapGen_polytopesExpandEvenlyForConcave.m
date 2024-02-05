function [ ...
exp_polytopes ...
] = ...
fcn_MapGen_polytopesExpandEvenlyForConcave( ...
polytopes, ...
exp_dist, ...
varargin...
)
% fcn_MapGen_polytopesExpandEvenly
% Expands an obstacle out by exp_dist on all sides.
%
%
%
% FORMAT:
%
%    [ ...
%    exp_polytopes ...
%    ] = ...
%    fcn_MapGen_polytopesExpandEvenly( ...
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
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_polytopesExpandEvenly
% for a full test suite.
%
% This function was written on 2018_11_17, Adjusted example code on 2021_04_28 by Seth Tau, Rebased on 2021_06_26 by S. Brennan by Seth Tau
% Questions or comments? contact sbrennan@psu.edu and sat5340@psu.edu

%
% REVISION HISTORY:
%
% 2018_11_17, Seth Tau
% -- first write of script
% 2021_04_28, Seth Tau
% -- Adjusted example code ,
% 2021_06_26 S. Brennan
% -- Rebased code
% -- Rewrote for clarity
% 2021_07_06 S. Brennan
% -- Vectorized plotting into array structure, to better support legends
% (rather than plotting all polytopes individually)

%
% TO DO:
%
% -- fill in to-do items here.

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    fig_for_debug = 680;
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
    % if nargin < 3 || nargin > 4
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end

    % Check the polytopes input, make sure it is 'polytopes' type
    fcn_MapGen_checkInputsToFunctions(...
        polytopes, 'polytopes');

    %     % Check the delta input, make sure it is 'positive_column_of_numbers' type
    %     fcn_MapGen_checkInputsToFunctions(...
    %         delta, 'positive_column_of_numbers',1);

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
    this_polytope = polytopes(p);
    this_polyshape = polyshape(this_polytope.vertices);
    scaled_polyshape = polybuffer(this_polyshape,exp_dist,'JointType','miter','MiterLimit',2);
    new_vertices = scaled_polyshape.Vertices;
    new_vertices = [new_vertices; new_vertices(1,:)];
    exp_polytopes(p).vertices = new_vertices;
end
exp_polytopes= fcn_MapGen_fillPolytopeFieldsFromVertices(exp_polytopes);
