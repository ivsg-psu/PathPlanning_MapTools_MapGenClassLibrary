function [ ...
exp_polytopes ...
] = ...
fcn_MapGen_polytopesExpandEvenly( ...
polytopes, ...
delta, ...
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
%     delta: a small number relative to vehicle size to determine the 
%     inside of an obstacle
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
%     fcn_MapGen_polytopeCentroidAndArea
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
% 2018_11_17, Adjusted example code on 2021_04_28 by Seth Tau, Rebased on 2021_06_26 by S. Brennan by Seth Tau
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
    fig_for_debug = 486;
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
    if nargin < 3 || nargin > 4
        error('Incorrect number of input arguments')
    end

    % Check the polytopes input, make sure it is 'polytopes' type
    fcn_MapGen_checkInputsToFunctions(...
        polytopes, 'polytopes');
 
    % Check the delta input, make sure it is 'column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        delta, 'column_of_numbers',1);
 
    % Check the exp_dist input, make sure it is 'column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        exp_dist, 'column_of_numbers',1);
 
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





exp_polytopes = polytopes; % both structures will be the same size
for polys = 1:size(polytopes,2) % check each obstacle
    %% pull values
    xvert = polytopes(polys).xv; % pull vertice values
    yvert = polytopes(polys).yv;
    numverts = length(polytopes(polys).xv); % find number of vertices
    xv = zeros(1,numverts); % pre-allocate xv and yv 
    yv = xv;
    for vert1 = 1:numverts % repeat for each vert
        %% assign 3 vertex indices
        if vert1 < numverts-1 
            vert2 = vert1 + 1;
            vert3 = vert2 + 1;
        elseif vert1 < numverts
            vert2 = vert1 + 1;
            vert3 = 1;
        else
            vert2 = 1;
            vert3 = vert2 + 1;
        end
        %% find which side of the line is in the obstacle
        ang1 = atan2(yvert(vert2)-yvert(vert1),xvert(vert2)-xvert(vert1)); % angle of first segment
        if vert1 == 1 % only execute once per obstacle
            mid_point = [(xvert(vert2)+xvert(vert1))/2; (yvert(vert2)+yvert(vert1))/2; 0]; % segment midpoint
            % check if turning ccw puts the point in the obstacle
            if inpolygon(mid_point(1)+delta*cos(ang1+pi/2),mid_point(2)+delta*sin(ang1+pi/2),xvert,yvert)
                turn = -pi/2; % turn the other way
            else % if puts point out of the obstacle
                turn = pi/2; % keep turning that way
            end
        end
        ang2 = atan2(yvert(vert3)-yvert(vert2),xvert(vert3)-xvert(vert2)); % angle of second segment
        %% force both angles to be positive
        if ang1 < 0
            ang1 = ang1 + 2*pi;
        end
        if ang2 < 0
            ang2 = ang2 + 2*pi;
        end
        %% move vertices
        exp_ang = (ang1 + ang2 + 2*turn)/2; % direction to move vertex
        inner_ang = exp_ang - ang2 - turn; % angle for determining distance
        dist = exp_dist/cos(inner_ang); % how far to move to ensure both sides move out by exp_dist
        xv(vert2) = polytopes(polys).xv(vert2) + dist*cos(exp_ang); % change x of vert2
        yv(vert2) = polytopes(polys).yv(vert2) + dist*sin(exp_ang); % change y of vert2
    end
    %% assign new exp_polytopes values
    exp_polytopes(polys).vertices = [[xv xv(1)]' [yv yv(1)]'];
    exp_polytopes(polys).xv = xv;
    exp_polytopes(polys).yv = yv;
    exp_polytopes(polys).distances = sum((exp_polytopes(polys).vertices(1:end-1,:)-exp_polytopes(polys).vertices(2:end,:)).^2,2).^0.5;
    [Cx,Cy,exp_polytopes(polys).area] = fcn_MapGen_polytopeCentroidAndArea([xv xv(1)]', [yv yv(1)]');
    exp_polytopes(polys).mean = [Cx, Cy];        
    exp_polytopes(polys).max_radius = max(sum((exp_polytopes(polys).vertices(1:end-1,:)-ones(length(xv),1)*exp_polytopes(polys).mean).^2,2).^0.5);
end

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

