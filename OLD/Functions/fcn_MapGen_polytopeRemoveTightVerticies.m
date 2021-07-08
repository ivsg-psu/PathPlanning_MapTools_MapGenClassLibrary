function [ ...
cleaned_polytope ...
] = ...
fcn_MapGen_polytopeRemoveTightVerticies( ...
polytopes, ...
tolerance, ...
varargin...
)
% fcn_MapGen_polytopeRemoveTightVerticies
% removes verticies of polytopes that are too close to each other, 
% measured by a tolerance
% 
% Sometimes, when shrinking, the new verticies are particularly close to 
% each other to where an edge has a trivial length. To prevent this, we 
% get rid of one of any vertices that are too close to each other. This 
% proximity is set by a user-defined tolerance.
% 
% FORMAT:
% 
%    [ ...
%    cleaned_polytope ...
%    ] = ...
%    fcn_MapGen_polytopeRemoveTightVerticies( ...
%    polytopes, ...
%    tolerance, ...
%    (fig_num) ...
%    )
% 
% INPUTS:
% 
%     polytopes: an individual structure or structure array of 'polytopes' 
%     type that stores the polytopes to be evaluated
% 
%     tolerance: a numeric value that defines how close points should be 
%     to be removed
% 
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a 
%     figure to be drawn showing results.
% 
% 
% OUTPUTS:
% 
%     cleaned_polytope: the resulting polytope after close edges are 
%     removed.
% 
% 
% DEPENDENCIES:
% 
%     fcn_MapGen_checkInputsToFunctions
% 
%     fcn_MapGen_fillPolytopeFieldsFromVerticies
% 
% 
% EXAMPLES:
% 
% See the script: script_test_fcn_MapGen_polytopeRemoveTightVerticies
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
    fig_for_debug = 380;
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
 
    % Check the tolerance input, make sure it is 'numeric' type
    fcn_MapGen_checkInputsToFunctions(...
        tolerance, 'numeric');
 
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





% pull values
vertices = polytope.vertices;
centroid = polytope.mean;


% Work with tolerance squared, since it avoids square-root calculations
tol_squared = tolerance.^2;

% find all vertices that are larger away from next one, relative to
% tolerance
good_ind = ...
    sum((vertices(1:end-1,:)-vertices(2:end,:)).^2,2)>(tol_squared);

% Check how many are good. If only 2 points are left, then the
% result is just a line. If 1 or zero, then result is just a point.
if sum(good_ind)>2 % sufficient good points to make a shape
    new_vert = vertices(good_ind,:);
elseif sum(good_ind)==2 % line shape
    new_vert = vertices(good_ind,:);
    
    % The line may not go through the centroid, which is odd. We force
    % this by removing the point closest to the centroid
    distances_to_centroid = sum((new_vert-centroid).^2,2).^0.5;
    if distances_to_centroid(1,1)>distances_to_centroid(2,1)
        new_vert(2,:)=centroid;
    else
        new_vert(1,:)=centroid;
    end
    
    new_vert = [new_vert; flipud(new_vert)];
else % singular shape (i.e. point) or no shape
    new_vert = [centroid; centroid; centroid];
end

% adjust polytopes
cleaned_polytope.vertices = [new_vert; new_vert(1,:)];
cleaned_polytope = fcn_MapGen_fillPolytopeFieldsFromVerticies(cleaned_polytope);


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

