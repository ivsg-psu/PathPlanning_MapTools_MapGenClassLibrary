function [polytopes, allVertices] = fcn_MapGen_generatePolysFromVoronoiAABBWithTiling(seedPoints, AABB, stretch, varargin)
% fcn_MapGen_generatePolysFromVoronoiAABBWithTiling
% creates polytopes given seed points, V and C matrices from Voronoi
% tiling, and stretch matrix. This function is very similar to:
% fcn_MapGen_generatePolysFromVoronoiAABB except that it tiles the points
% such that the resulting polytope map is itself completely tile-symmetric,
% e.g. it also tiles correctly with itself. This avoids particular edge
% issues that occur when using fcn_MapGen_generatePolysFromVoronoiAABB.
%
% FORMAT:
%
%    [ ...
%    polytopes ...
%    ] = ...
%    fcn_MapGen_generatePolysFromVoronoiAABBWithTiling( ...
%    seedPoints, ...
%    AABB, ...
%    stretch, ...
%    (fig_num) ...
%    )
%
% INPUTS:
%
%     seedPoints: the list of seed points in [x y] format, where x and y
%     are columns
%
%     AABB: the axis-aligned bounding box, in format of
%     [xmin ymin xmax ymax], wherein the resulting polytopes must be
%     bounded.
%
%     stretch: the factor to stretch the polytopes after they are
%     calculated, in format of [xmagnification ymagnification]
%
%     (optional inputs)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
% OUTPUTS:
%
%     polytopes: the resulting polytopes after converting to polytope form.
%
%     allVertices: a Mx3 array with M>N where N is number of polytopes
%     containing NaN values separating rows of real data. Each column
%     consists of [polyID vertexX vertexY] data, when not nan values.
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_tilePoints
%     fcn_MapGen_polytopeFillEmptyPoly
%     fcn_MapGen_plotPolytopes
%     fcn_MapGen_fillPolytopeFieldsFromVertices
%     fcn_MapGen_polytopesShrinkFromEdges
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_generatePolysFromVoronoiAABBWithTiling
% for a full test suite.
%
% This function was written on 2021_07_02 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

%
% REVISION HISTORY:
%
% 2021_07_02 by Sean Brennan
% -- first write of function
% 2021_07_30 by Sean Brennan
% -- fixed errors due to corners being omitted
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% -- fixed docstrings to correct argument listing
% 2025_07_15 by Sean Brennan
% -- cleaned variable naming to remove underscores
% -- turned on fast mode for subfunctions

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS");
    MATLABFLAG_MAPGEN_FLAG_DO_DEBUG = getenv("MATLABFLAG_MAPGEN_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_MAPGEN_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_MAPGEN_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
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

if (0==flag_max_speed)
    if 1 == flag_check_inputs

        % Are there the right number of inputs?
        narginchk(3,4);

        % Check the seedPoints input, make sure it is '2column_of_numbers' type
        fcn_DebugTools_checkInputsToFunctions(...
            seedPoints, '2column_of_numbers');

        % Check the AABB input, make sure it is '4column_of_numbers' type
        fcn_DebugTools_checkInputsToFunctions(...
            AABB, '4column_of_numbers',1);

        % Check the stretch input, make sure it is '2column_of_numbers' type
        fcn_DebugTools_checkInputsToFunctions(...
            stretch, '2column_of_numbers',1);

    end
end


% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  4 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
    temp = varargin{end}; % Last argument is always figure number
    if ~isempty(temp) % Make sure the user is not giving empty input
        fig_num = temp;
        flag_do_plot = 1; % Set flag to do plotting
    end
else
    if flag_do_debug % If in debug mode, do plotting but to an arbitrary figure number
        fig = figure;
        fig_for_debug = fig.Number; %#ok<NASGU>
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

 % Put one copy of the points in a perimeter around the points, tiled by
 % the AABB
tile_depth = 1;
Nseedpoints = length(seedPoints);
[tiled_original_seedPoints] = fcn_MapGen_tilePoints(seedPoints,tile_depth,AABB, -1);

% Calculate the resulting Voronoi diagram
[V,C] = voronoin(tiled_original_seedPoints);

% Choose the seed points from the middle area
mid_tile_superindex = ceil((2*tile_depth+1)^2/2);
mid_tile_range = ((Nseedpoints*(mid_tile_superindex-1)+1):Nseedpoints*mid_tile_superindex)';
seedPoints_to_use = tiled_original_seedPoints(mid_tile_range,:);

%% for each of the seed points, save the polytope for it
% Start by initializing the data structure
num_poly = size(seedPoints_to_use,1);
polytopes(num_poly) = fcn_MapGen_polytopeFillEmptyPoly((-1));
Npolys = length(polytopes);
Nvertices_per_poly = 20; % Maximum estimate
Nvertices_per_map = Npolys*Nvertices_per_poly;
allVertices = nan(Nvertices_per_map,3);
% all_neighbors = nan(Nvertices_per_map,1);


%% Loop through the polytopes, filling seed point, verticies
% (and neighbors matrix?)

% Prep a plot?
if flag_do_debug
    fig_debug = 12345;
    figure(fig_debug);
    clf;
    hold on;

    scale = max(3*AABB,[],'all') - min(3*AABB,[],'all');
    new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];
    axis(new_axis);
end

% Start the loop
for ith_poly = 1:length(seedPoints)
    % Find which seed index to use
    offset_seed_index = mid_tile_range(ith_poly);

    % Fill in seed_point
    polytopes(ith_poly).seed_point = seedPoints_to_use(ith_poly);


    %     if ith_poly==50
    %         disp('Stop here');
    %     end

    vertices_open = V(C{offset_seed_index},:);

    % Remove ill-conditioned points by setting the ones really, really far
    % away to infinity
    scale = max(AABB,[],'all') - min(AABB,[],'all');
    center = [AABB(1)+AABB(3), AABB(2)+AABB(4)]/2;
    distances_from_center = sum((vertices_open-center).^2,2).^0.5;
    near_infinite = (distances_from_center/scale)>1E7;
    if any(near_infinite)
        warning('on','backtrace');
        warning('Near-infinite vertex found for polytope: %.0d, for seed point: (%.3f, %.3f)', ith_poly,seedPoints_to_use(ith_poly,1),seedPoints_to_use(ith_poly,2))
        vertices_open(near_infinite,:) = inf;
        % Remove repeated infinities
        vertices_open = unique(vertices_open,'rows','stable');
    end

    % Append results to close off the vector loop
    vertices = [vertices_open; vertices_open(1,:)]; % Close off the vertices

    % Save verticies to this polytope
    polytopes(ith_poly).vertices = vertices;

    % Plot the poly?
    if flag_do_debug

        % Replot the poly before this one in blue (to cover it up)
        line_width = 3;
        if ith_poly>1
            % fcn_MapGen_OLD_plotPolytopes(polytopes(ith_poly-1),fig_debug,'b-',line_width);
            plotFormat.LineWidth = line_width;
            plotFormat.MarkerSize = 10;
            plotFormat.LineStyle = '-';
            plotFormat.Color = [0 0 1];
            fillFormat = [];
            h_plot = fcn_MapGen_plotPolytopes(polytopes(ith_poly-1), (plotFormat),(fillFormat),(fig_debug)); %#ok<NASGU>
        end

        % Plot this one in red
        % fcn_MapGen_OLD_plotPolytopes(polytopes(ith_poly),fig_debug,'r-',line_width);
        plotFormat.LineWidth = line_width;
        plotFormat.MarkerSize = 10;
        plotFormat.LineStyle = '-';
        plotFormat.Color = [1 0 0];
        fillFormat = [];
        h_plot = fcn_MapGen_plotPolytopes(polytopes(ith_poly), (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>
        title(sprintf('Polytope: %.0d',ith_poly));

        % Is it the last poly? Plot it in blue
        if ith_poly==Npolys
            % fcn_MapGen_OLD_plotPolytopes(polytopes(ith_poly),fig_debug,'b-',line_width);
            plotFormat.LineWidth = line_width;
            plotFormat.MarkerSize = 10;
            plotFormat.LineStyle = '-';
            plotFormat.Color = [0 0 1];
            fillFormat = [];
            h_plot = fcn_MapGen_plotPolytopes(polytopes(ith_poly), (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        end

        % Plot this one's seed point:
        plot(seedPoints_to_use(ith_poly,1),seedPoints_to_use(ith_poly,2),'r.');

    end

    % Save verticies to the all_verticies array
    Nvertices = length(vertices(:,1));
    if Nvertices>Nvertices_per_poly
        error('Need to resize the number of allowable vertices');
    else
        row_offset = (ith_poly-1)*Nvertices_per_poly;
        allVertices(row_offset+1:row_offset+Nvertices,1) = ith_poly;
        allVertices(row_offset+1:row_offset+Nvertices,2:3) = vertices;
    end
end


%% Apply the stretch
num_poly = length(polytopes);
for poly = 1:num_poly % pull each cell from the voronoi diagram
    try
        polytopes(poly).vertices  = polytopes(poly).vertices.*stretch;
    catch
        error('stop here');
    end
end % Ends for loop for stretch

%% Fill in all the other fields
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes, [], -1);

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

    figure(fig_num);
    clf;
    hold on
    scale = max(AABB,[],'all') - min(AABB,[],'all');
    new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];
    axis(new_axis);

    % plot the polytopes    
    % fcn_MapGen_OLD_plotPolytopes(polytopes,fig_num,'b',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [0 0 1];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

    % plot all vertices
    plot(allVertices(:,2),allVertices(:,3),'-', 'Color',0.8*[1 1 1],'Linewidth',1, 'MarkerSize',20, 'DisplayName','vertices');

    % plot the seed points in red
    plot(seedPoints(:,1),seedPoints(:,2),'r.','Markersize',10, 'DisplayName', 'seedPoints');

    % plot the means in black
    temp = zeros(length(polytopes),2);
    for ith_poly = 1:length(polytopes)
        temp(ith_poly,:) = polytopes(ith_poly).mean;
    end
    plot(temp(:,1),temp(:,2),'ko','Markersize',3,'DisplayName', 'polytope means');

    legend('Interpreter','none','location','best');

    % number the polytopes at seed points?
    if 1== flag_do_plot
        for ith_poly = 1:length(polytopes)
            text_location = seedPoints(ith_poly,:);
            text(text_location(1,1),text_location(1,2),sprintf('%.0d',ith_poly));
        end
    end

    %     % number the polytopes at means
    %     for ith_poly = 1:length(polytopes)
    %         text_location = polytopes(ith_poly).mean;
    %         text(text_location(1,1),text_location(1,2),sprintf('%.0d',ith_poly));
    %     end

    % plot the connections between the polytope neighbors
    if 1==0
        % Clean up and sort the vertices so that we can associate neighbors
        allVertices_no_nan = allVertices(~isnan(allVertices(:,1)),:);
        sorted_allVertices = sortrows(allVertices_no_nan,[2 3]);

        % Remove repeats
        sorted_allVertices = unique(sorted_allVertices,'rows','stable');

        % Remove infinities
        sorted_allVertices = sorted_allVertices(~isinf(sorted_allVertices(:,2)));

        Nrealvertices = floor(length(sorted_allVertices(:,1))/3);
        data = zeros(Nrealvertices*6,2);
        for ith_poly = 1:Nrealvertices
            row_offset = (ith_poly-1)*3;
            neighbors = sorted_allVertices(row_offset+1:row_offset+3,1);

            for jth_neighbor = 2:length(neighbors)
                neigh_offset = (ith_poly-1)*6 + ((jth_neighbor-2)*3);
                data(neigh_offset+1:neigh_offset+3,:) = [seedPoints(neighbors(1),:); seedPoints(neighbors(jth_neighbor),:); nan nan];
            end
        end
        plot(data(:,1),data(:,2),'-','Linewidth',0.5);
    end


    %% Show a detailed step-by-step process behind construction of obstacle map
    if 1==0
        % using fcn_MapGen_generatePolysFromVoronoiAABB
        fig_num = 1010;
        figure(fig_num);
        clf;

        % Calculate the scale
        scale = max(AABB,[],'all') - min(AABB,[],'all');
        new_axis = [AABB(1)-1.5*scale AABB(3)+1.5*scale AABB(2)-1.5*scale AABB(4)+1.5*scale];


        %% plot the seed points in red
        subplot(2,3,1);

        plot(seedPoints(:,1),seedPoints(:,2),'r.','Markersize',10);


        % number the polytopes at seed points
        for ith_poly = 1:length(polytopes)
            text_location = seedPoints(ith_poly,:);
            text(text_location(1,1),text_location(1,2),sprintf('%.0d',ith_poly));
        end
        axis(new_axis);
        title('Seed points');

        %% plot the tiled seed points
        subplot(2,3,2);
        cla;
        plot(tiled_original_seedPoints(:,1),tiled_original_seedPoints(:,2),'b.','Markersize',20);
        hold on;
        plot(seedPoints(:,1),seedPoints(:,2),'r.','Markersize',10);

        axis(new_axis);
        title('Tiled Seed points');


        %% PLOT THE VORONOI lines with the points
        subplot(2,3,3);

        % Start the loop to calculate all the vertices
        num_poly = size(tiled_original_seedPoints,1);
        voronoi_polytopes(num_poly) = fcn_MapGen_polytopeFillEmptyPoly((-1));

        Npolys = length(voronoi_polytopes);
        Nvertices_per_poly = 20; % Maximum estimate
        Nvertices_per_map = Npolys*Nvertices_per_poly;
        all_voronoi_vertices = nan(Nvertices_per_map,3);

        for ith_poly = 1:length(tiled_original_seedPoints(:,1))

            % Fill in seed_point
            voronoi_polytopes(ith_poly).seed_point = tiled_original_seedPoints(ith_poly);

            % Get the verticies
            vertices_open = V(C{ith_poly},:);

            % Append results to close off the vector loop
            vertices = [vertices_open; vertices_open(1,:)]; % Close off the vertices

            % Save verticies to this polytope
            voronoi_polytopes(ith_poly).vertices = vertices;

            % Save verticies to the all_verticies array
            Nvertices = length(vertices(:,1));
            if Nvertices>Nvertices_per_poly
                error('Need to resize the number of allowable vertices');
            else
                row_offset = (ith_poly-1)*Nvertices_per_poly;
                all_voronoi_vertices(row_offset+1:row_offset+Nvertices,1) = ith_poly;
                all_voronoi_vertices(row_offset+1:row_offset+Nvertices,2:3) = vertices;
            end
        end


        % plot the tiled seed points in blue
        plot(tiled_original_seedPoints(:,1),tiled_original_seedPoints(:,2),'b.','Markersize',20);
        hold on;

        % plot all vertices
        plot(all_voronoi_vertices(:,2),all_voronoi_vertices(:,3),'c','Linewidth',1);

        axis(new_axis);
        title('Voronoi boundaries');


        %% PLOT THE VORONOI lines for just the middle
        subplot(2,3,4);

        % plot the polytopes on current axis
        % fcn_MapGen_OLD_plotPolytopes(polytopes,gca,'b',2);
        plotFormat.LineWidth = 2;
        plotFormat.MarkerSize = 10;
        plotFormat.LineStyle = '-';
        plotFormat.Color = [0 0 1];
        fillFormat = [];
        h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        hold on;


        % plot the seed points in red
        plot(seedPoints(:,1),seedPoints(:,2),'r.','Markersize',10);

        axis(new_axis);
        title('Extract middle, tilable polytopes')


        %% show that it's tilable
        subplot(2,3,5);

        % plot the polytopes
        % fcn_MapGen_OLD_plotPolytopes(polytopes,gca,'b',5);
        plotFormat.LineWidth = 5;
        plotFormat.MarkerSize = 10;
        plotFormat.LineStyle = '-';
        plotFormat.Color = [0 0 1];
        fillFormat = [];
        h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>


        hold on;

        xscale = AABB(3)-AABB(1);
        yscale = AABB(4)-AABB(2);

        plotFormat.LineWidth = 2;
        plotFormat.MarkerSize = 10;
        plotFormat.LineStyle = '-';
        fillFormat = [];

        shifted_left_up = fcn_INTERNAL_shift_polys(polytopes,-xscale,yscale);
        % fcn_MapGen_OLD_plotPolytopes(shifted_left_up,gca,'-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_left_up, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_left = fcn_INTERNAL_shift_polys(polytopes,-xscale,0);
        % fcn_MapGen_OLD_plotPolytopes(shifted_left,gca,'-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_left, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_left_down = fcn_INTERNAL_shift_polys(polytopes,-xscale,-yscale);
        % fcn_MapGen_OLD_plotPolytopes(shifted_left_down,gca,'-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_left_down, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_up = fcn_INTERNAL_shift_polys(polytopes,0,yscale);
        % fcn_MapGen_OLD_plotPolytopes(shifted_up,gca,'-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_up, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_down = fcn_INTERNAL_shift_polys(polytopes,0,-yscale);
        % fcn_MapGen_OLD_plotPolytopes(shifted_down,gca,'-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_down, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_right_up = fcn_INTERNAL_shift_polys(polytopes,xscale,yscale);
        % fcn_MapGen_OLD_plotPolytopes(shifted_right_up,gca,'-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_right_up, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_right = fcn_INTERNAL_shift_polys(polytopes,xscale,0);
        % fcn_MapGen_OLD_plotPolytopes(shifted_right,gca,'-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_right, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_right_down = fcn_INTERNAL_shift_polys(polytopes,xscale,-yscale);
        % fcn_MapGen_OLD_plotPolytopes(shifted_right_down,gca,'-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_right_down, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        axis(new_axis);
        title('Polytopes can tile');

        %% plot the shrink to edge
        subplot(2,3,6);

        des_gap_size = 0.05;

        shrunk_polytopes=...
            fcn_MapGen_polytopesShrinkFromEdges(...
            polytopes,des_gap_size, -1);

        % plot the shrunk polytopes
        fcn_MapGen_plotPolytopes(shrunk_polytopes,gca,'r',2);
        hold on;

        plotFormat.LineWidth = 5;
        plotFormat.MarkerSize = 10;
        plotFormat.LineStyle = '-';
        plotFormat.Color = [0 0 1];
        fillFormat = [];        

        shifted_left_up = fcn_INTERNAL_shift_polys(shrunk_polytopes,-xscale,yscale);
        % fcn_MapGen_OLD_plotPolytopes(shifted_left_up,gca,'b-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_left_up, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_left = fcn_INTERNAL_shift_polys(shrunk_polytopes,-xscale,0);
        % fcn_MapGen_OLD_plotPolytopes(shifted_left,gca,'b-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_left, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_left_down = fcn_INTERNAL_shift_polys(shrunk_polytopes,-xscale,-yscale);
        % fcn_MapGen_OLD_plotPolytopes(shifted_left_down,gca,'b-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_left_down, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_up = fcn_INTERNAL_shift_polys(shrunk_polytopes,0,yscale);
        % fcn_MapGen_OLD_plotPolytopes(shifted_up,gca,'b-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_up, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_down = fcn_INTERNAL_shift_polys(shrunk_polytopes,0,-yscale);
        % fcn_MapGen_OLD_plotPolytopes(shifted_down,gca,'b-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_down, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_right_up = fcn_INTERNAL_shift_polys(shrunk_polytopes,xscale,yscale);
        % fcn_MapGen_OLD_plotPolytopes(shifted_right_up,gca,'b-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_right_up, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_right = fcn_INTERNAL_shift_polys(shrunk_polytopes,xscale,0);
        % fcn_MapGen_OLD_plotPolytopes(shifted_right,gca,'b-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_right, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        shifted_right_down = fcn_INTERNAL_shift_polys(shrunk_polytopes,xscale,-yscale);
        % fcn_MapGen_OLD_plotPolytopes(shifted_right_down,gca,'b-',2);
        h_plot = fcn_MapGen_plotPolytopes(shifted_right_down, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

        axis(new_axis);
        title('Shrunk from edge polytopes, with tiling');
    end

    if 1==0
        %% Confirm the tiling
        % if it is a true tiling, then the parts that stick out, e.g. values
        % greater than the AABB boundary, should match points that are in the
        % interior - in other words, there is a seamless connection between the two
        % areas.

        % Find the parts that stick out
        all_points = allVertices(~isnan(allVertices(:,1)),2:3);
        all_points_shifted = all_points - AABB(1:2);
        rounded_points = mod(all_points_shifted,1);
        x_indices_stick_out = find(all_points_shifted(:,1)~=rounded_points(:,1));
        y_indices_stick_out = find(all_points_shifted(:,2)~=rounded_points(:,2));
        xy_indices_stick_out = intersect(x_indices_stick_out,y_indices_stick_out);

        % Find the points that are entirely within the AABB
        x_indices_inside = find(all_points_shifted(:,1)==rounded_points(:,1));
        y_indices_inside = find(all_points_shifted(:,2)==rounded_points(:,2));
        xy_indices_inside = intersect(x_indices_inside,y_indices_inside);
        points_inside = all_points(xy_indices_inside,:);

        % DEBUG AREA
        tiling_fig_num = 12312;
        figure(tiling_fig_num);
        clf;
        hold on

        scale = max(AABB,[],'all') - min(AABB,[],'all');
        new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];
        axis(new_axis);

        % plot the polytopes
        % fcn_MapGen_OLD_plotPolytopes(polytopes,tiling_fig_num,'b',2);
        plotFormat.LineWidth = 2;
        plotFormat.MarkerSize = 10;
        plotFormat.LineStyle = '-';
        plotFormat.Color = [0 0 1];
        fillFormat = [];
        h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat),(fillFormat),(fig_debug)); %#ok<NASGU>

        % plot all vertices
        plot(all_points(:,1),all_points(:,2),'c.','Linewidth',1);

        % Plot all the points that stick out
        plot(all_points(x_indices_stick_out,1),all_points(x_indices_stick_out,2),'r.','Markersize',20);
        plot(all_points(y_indices_stick_out,1),all_points(y_indices_stick_out,2),'m.','Markersize',15);
        plot(all_points(xy_indices_stick_out,1),all_points(xy_indices_stick_out,2),'k.','Markersize',10);
        legend('Polytopes','All points','X-data out of bounds','Y-data out of bounds','Both X and Y out of bounds');

        % For each of the "stick out" points, check to see that there's an inside
        % point that corresponds to the same location, but inside.

        legend off;

        % Start with x
        error_inside_to_outside_x = nan(length(x_indices_stick_out),1);
        for ith_outside_point = 1:length(x_indices_stick_out)
            current_outside_point = all_points(x_indices_stick_out(ith_outside_point),:);
            rounded_current_outside_point = mod(current_outside_point,1);

            % Use a distance metric to find closest actual point
            tolerance = 1E-10;
            distances_to_outside_point = sum((rounded_current_outside_point-points_inside).^2,2).^0.5;
            [min_distance,~] = min(distances_to_outside_point);
            index_min = distances_to_outside_point<tolerance;
            closest_points = points_inside(index_min,:);

            % Show the results
            plot(current_outside_point(:,1),current_outside_point(:,2),'ro','Markersize',20);
            plot(closest_points(:,1),closest_points(:,2),'ro','Markersize',20);
            error_inside_to_outside_x(ith_outside_point)= min_distance;
        end

        % Now with y
        error_inside_to_outside_y = nan(length(y_indices_stick_out),1);
        for ith_outside_point = 1:length(y_indices_stick_out)
            current_outside_point = all_points(y_indices_stick_out(ith_outside_point),:);
            rounded_current_outside_point = mod(current_outside_point,1);

            % Use a distance metric to find closest actual point
            tolerance = 1E-10;
            distances_to_outside_point = sum((rounded_current_outside_point-points_inside).^2,2).^0.5;
            [min_distance,~] = min(distances_to_outside_point);
            index_min = distances_to_outside_point<tolerance;
            closest_points = points_inside(index_min,:);

            % Show the results
            plot(current_outside_point(:,1),current_outside_point(:,2),'go','Markersize',20);
            plot(closest_points(:,1),closest_points(:,2),'go','Markersize',20);
            error_inside_to_outside_y(ith_outside_point)= min_distance;
        end
        figure(1234);
        clf;
        hold on;
        plot(error_inside_to_outside_x);
        plot(error_inside_to_outside_y);
        legend('X errors', 'Y errors');
    end

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

function shifted_polytopes = fcn_INTERNAL_shift_polys(polytopes,xshift,yshift)
shifted_polytopes = polytopes;
for ith_poly = 1:length(polytopes)
    % Shift verticies
    shifted_polytopes(ith_poly).vertices = polytopes(ith_poly).vertices + [xshift yshift];
end
end
