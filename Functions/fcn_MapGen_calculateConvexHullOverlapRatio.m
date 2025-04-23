function [ ...
convex_hull_overlap_ratio,...
A_overlap,...
A_occupied...
] = ...
fcn_MapGen_calculateConvexHullOverlapRatio( ...
polytopes, ...
varargin...
)
% fcn_MapGen_calculateConvexHullOverlapRatio
% calculates the convex hull of every obstacle.  The area of the overlap between
% these hulls relative to the total occupied area
%
%
%
% FORMAT:
%
% function [ ...
% convex_hull_overlap_ratio...
% ] = ...
% fcn_MapGen_calculateConvexHullOverlapRatio( ...
% polytopes, ...
% varargin...
% )
%
% INPUTS:
%
%     polytopes: the structure of 'polytopes' type that stores the
%     polytopes to be expanded
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
%     covex_hull_overlap_ratio: portion of overlapping convex hull area to total obstacle area
%
%
% DEPENDENCIES:
%
%     MATLAB's polyshape object and union object function (method)
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_polytopesExpandEvenlyForConcave
% for a full test suite.
%
% This function was written 23 Feb. 2024 by Steve Harnett
% Questions or comments? contact sjharnett@psu.edu

%
% REVISION HISTORY:
%
% 2024_02_23, Steve Harnett
% -- first write of function
% 2025_04_16, Steve Harnett
% -- add legend to plotting

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
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end

    % Check the polytopes input, make sure it is 'polytopes' type
    fcn_MapGen_checkInputsToFunctions(...
        polytopes, 'polytopes');


end

% Does user want to show the plots?
if  2== nargin
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
A_overlap = 0;
A_occupied = 0;
conv_hull_polyshapes = [];
for p = 1:length(polytopes)
    this_polytope = polytopes(p); % look at one polytope
    if flag_do_plot
        figure(fig_num); hold on; box on; fill(this_polytope.vertices(:,1)',this_polytope.vertices(:,2),[0 0 1],'FaceAlpha',1);
        if p == 1
            leg_str = {'obstacles'};
        else
            leg_str{end+1} = '';
        end
    end
    these_vertices = this_polytope.vertices(1:(end-1),:); % grab only non-repeating vertices
    k = convhull(these_vertices); % find convex hull of vertices
    convex_hull_vertices = [these_vertices(k,1),these_vertices(k,2)];
    this_conv_hull_polyshape = polyshape(convex_hull_vertices); % convert it to matlab polyshape
    if flag_do_plot
        figure(fig_num); hold on; box on; plot(this_conv_hull_polyshape,'FaceColor','green','FaceAlpha',0.2);
        if p == 1
            leg_str{end+1} = 'obs. convex hull';
        else
            leg_str{end+1} = '';
        end
    end
    conv_hull_polyshapes = [conv_hull_polyshapes; this_conv_hull_polyshape];
    A_occupied = A_occupied+this_polytope.area;
end
flag_havent_plotted = 1;
for i = 1:length(polytopes)
    j = 1;
    while j <= length(polytopes)
        if j<=i
            % if we checked 1,2 we don't need to check 2,1 so ignore j<i
            % also don't want to check j=i because there is meaningless overlap
            j = j+1;
            continue;
        end
        overlap_polyshape = intersect(conv_hull_polyshapes(i),conv_hull_polyshapes(j));
        if flag_do_plot
            figure(fig_num); hold on; box on; plot(overlap_polyshape,'FaceColor','red');
            if flag_havent_plotted
                leg_str{end+1} = 'overlap';
                flag_havent_plotted = 0;
            end
        end
        A_overlap = A_overlap + area(overlap_polyshape);
        j = j+1;
    end
end

convex_hull_overlap_ratio = A_overlap/A_occupied;
if flag_do_plot
    title_str = sprintf('total obstacle area: %.3f\noverlapping convex hull area: %.3f\nconvex hull overlap ratio: %.3f',A_occupied,A_overlap,convex_hull_overlap_ratio);
    figure(fig_num); title(title_str);
    legend(leg_str)
end
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


end % Ends the function
