% script_test_fcn_MapGen_polytopeFindVertexSkeleton
% Tests function: fcn_MapGen_polytopeFindVertexSkeleton

% REVISION HISTORY:
% 2022_02_15
% -- first written by S. Brennan

fig_num = 675;
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num);

figure(4747);
grid on
grid minor
hold on
axis equal

for cut = 0:0.5:cut_distance(end)
    
    % Find the shape that is less than or equal to the cut
    shape_index = find(cut_distance<=cut,1,'last');
    
    % Grab vertices to start from, cut to start from
    starting_vertices = new_vertices{shape_index};
    starting_cut = cut_distance(shape_index);
    
    % Calculate projection distance
    projection_distance = cut - starting_cut;
    
    % Determine final vertices
    final_vertices = starting_vertices + new_projection_vectors{shape_index}*projection_distance;
    
    % Plot results
    plot(final_vertices(:,1),final_vertices(:,2),'r.-','Linewidth',2,'Markersize',20);
end



%% Example case 1: non-normal wall shrinking
fig_num = 675;
% this polytope has a vertical wall
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num); %#ok<*ASGLU>


%% Example case 2: square
fig_num = 234343;
vertices = [0 0; 1 0; 1 1; 0 1; 0 0]*10;
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num);


% assert that final vertices are within 5% error of having the same x
% position
final_vertices = new_vertices{end}(1,:);
expected_final_vertices = [5 5];
error_tolerance = 0.05;
vertical_error = sum(sum((final_vertices - expected_final_vertices).^2,1).^0.5,2);
% note this error is currently about 33% as of the commit this line was
% added
assert(vertical_error <= error_tolerance,['Expected and actual final points are: ',...
    '\n\t(%d,%d)  (expected) \n\t(%d,%d) (actual) \nyielding a vertical error of %d. Error tolerance was %d.\n'],...
    expected_final_vertices(1,1), expected_final_vertices(1,2),...
    final_vertices(1,1),final_vertices(1,2),vertical_error,error_tolerance);

%% Example case 3: wide rectangle
fig_num = 2464;
vertices = [0 0; 1 0; 1 0.5; 0 0.5; 0 0]*10;
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num);


% assert that final vertices are within 5% error of having the same x
% position
final_vertices = new_vertices{end}(1,:);
expected_final_vertices = [5 2.5];
error_tolerance = 0.05;
vertical_error = sum(sum((final_vertices - expected_final_vertices).^2,1).^0.5,2);
% note this error is currently about 33% as of the commit this line was
% added
assert(vertical_error <= error_tolerance,['Expected and actual final points are: ',...
    '\n\t(%d,%d)  (expected) \n\t(%d,%d) (actual) \nyielding a vertical error of %d. Error tolerance was %d.\n'],...
    expected_final_vertices(1,1), expected_final_vertices(1,2),...
    final_vertices(1,1),final_vertices(1,2),vertical_error,error_tolerance);


%% Example case 4: tall rectangle
fig_num = 2465;
vertices = [0 0; 0.5 0; 0.5 1; 0 1; 0 0]*10;
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num);

% assert that final vertices are within 5% error of having the same x
% position
final_vertices = new_vertices{end}(1,:);
expected_final_vertices = [2.5 5];
error_tolerance = 0.05;
vertical_error = sum(sum((final_vertices - expected_final_vertices).^2,1).^0.5,2);
% note this error is currently about 33% as of the commit this line was
% added
assert(vertical_error <= error_tolerance,['Expected and actual final points are: ',...
    '\n\t(%d,%d)  (expected) \n\t(%d,%d) (actual) \nyielding a vertical error of %d. Error tolerance was %d.\n'],...
    expected_final_vertices(1,1), expected_final_vertices(1,2),...
    final_vertices(1,1),final_vertices(1,2),vertical_error,error_tolerance);


%% Example case 5: goofy polytope
fig_num = 47464;
vertices = [0 0; 10 0; 5 15; 4 17; 1 13; 0 5; 0 0];
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(vertices,fig_num);



%% Random polytope calculation
% Set up polytopes
close all;

seedGeneratorNames = 'haltonset';
seedGeneratorRanges = [1 100];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (-1));


edge_cut_step = 0.005;
bounding_box = [0,0, 1,1];
trim_polytopes = fcn_MapGen_polytopesDeleteByAABB(polytopes,bounding_box);

% Pick a random polytope
Npolys = length(trim_polytopes);
rand_poly = 1+floor(rand*Npolys);
shrinker = trim_polytopes(rand_poly);

% Do skeleton calculation
[new_vertices, new_projection_vectors, cut_distance] = fcn_MapGen_polytopeFindVertexSkeleton(shrinker.vertices,2727);


fig_num = 11;
figure(fig_num); 
clf;
axis equal;
hold on;

% for edge_cut = edge_cut_step:edge_cut_step:(shrinker.max_radius+edge_cut_step)
for cut = edge_cut_step:edge_cut_step:(shrinker.max_radius+edge_cut_step)
    
    % Find the shape that is less than or equal to the cut
    shape_index = find(cut_distance<=cut,1,'last');
    
    % Grab vertices to start from, cut to start from
    starting_vertices = new_vertices{shape_index};
    starting_cut = cut_distance(shape_index);
    
    % Calculate projection distance
    projection_distance = cut - starting_cut;
    
    % Determine final vertices
    final_vertices = starting_vertices + new_projection_vectors{shape_index}*projection_distance;
    
    % Plot results
    plot(final_vertices(:,1),final_vertices(:,2),'.-','Linewidth',2,'Markersize',20);
end

% script_test_fcn_MapGen_polytopeFindSelfIntersections
% Tests function: fcn_MapGen_polytopeFindSelfIntersections

% REVISION HISTORY:
% 2021_08_03
% -- first written by S. Brennan
% 2025_07_11 - S. Brennan, sbrennan@psu.edu
% -- updated script testing to standard form

%% Set up the workspace
close all

%% Code demos start here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                              ____   __    _____          _
%  |  __ \                            / __ \ / _|  / ____|        | |
%  | |  | | ___ _ __ ___   ___  ___  | |  | | |_  | |     ___   __| | ___
%  | |  | |/ _ \ '_ ` _ \ / _ \/ __| | |  | |  _| | |    / _ \ / _` |/ _ \
%  | |__| |  __/ | | | | | (_) \__ \ | |__| | |   | |___| (_) | (_| |  __/
%  |_____/ \___|_| |_| |_|\___/|___/  \____/|_|    \_____\___/ \__,_|\___|
%
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Demos%20Of%20Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 1

close all;
fprintf(1,'Figure: 1XXXXXX: DEMO cases\n');

%% DEMO case: self-intersection
fig_num = 10001;
titleString = sprintf('DEMO case: self-intersection');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

vertices = [0 0; 1 0; 0.5 1.5; 1 1; 0 1; 0 0];
verticesIncludingSelfIntersections = fcn_MapGen_polytopeFindSelfIntersections(...
    vertices, -1);

interiorPoint = [0.5 0.5];

% Call the function
[projectedPoints] = ...
    fcn_MapGen_polytopeProjectVerticesOntoWalls(...,
    interiorPoint,...
    verticesIncludingSelfIntersections,...
    verticesIncludingSelfIntersections(1:end-1,:),...
    verticesIncludingSelfIntersections(2:end,:),...
    (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(projectedPoints));

% Check variable sizes
Nvertices = length(verticesIncludingSelfIntersections(:,1));
assert(size(projectedPoints,1)==Nvertices);
assert(size(projectedPoints,2)==2);

% Check variable values
assert(isequal(round(projectedPoints,4),round(...
    [...
    0         0
    0         0
    1.0000         0
    0.7500    0.7500
    0.6667    1.0000
    0.6667    1.0000
    0.5000    1.0000
    0    1.0000
    ]...
    ,4)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test cases start here. These are very simple, usually trivial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  _______ ______  _____ _______ _____
% |__   __|  ____|/ ____|__   __/ ____|
%    | |  | |__  | (___    | | | (___
%    | |  |  __|  \___ \   | |  \___ \
%    | |  | |____ ____) |  | |  ____) |
%    |_|  |______|_____/   |_| |_____/
%
%
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=TESTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 2

close all;
fprintf(1,'Figure: 2XXXXXX: TEST mode cases\n');
% 
% %% TEST case: simple crossing at origin
% fig_num = 20001;
% titleString = sprintf('TEST case: simple crossing at origin');
% fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
% figure(fig_num); clf;


%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 8

close all;
fprintf(1,'Figure: 8XXXXXX: FAST mode cases\n');

%% Basic example - NO FIGURE
fig_num = 80001;
fprintf(1,'Figure: %.0f: FAST mode, empty fig_num\n',fig_num);
figure(fig_num); close(fig_num);

vertices = [0 0; 1 0; 0.5 1.5; 1 1; 0 1; 0 0];
verticesIncludingSelfIntersections = fcn_MapGen_polytopeFindSelfIntersections(...
    vertices, -1);

interiorPoint = [0.5 0.5];

% Call the function
[projectedPoints] = ...
    fcn_MapGen_polytopeProjectVerticesOntoWalls(...,
    interiorPoint,...
    verticesIncludingSelfIntersections,...
    verticesIncludingSelfIntersections(1:end-1,:),...
    verticesIncludingSelfIntersections(2:end,:),...
    ([]));

% Check variable types
assert(isnumeric(projectedPoints));

% Check variable sizes
Nvertices = length(verticesIncludingSelfIntersections(:,1));
assert(size(projectedPoints,1)==Nvertices);
assert(size(projectedPoints,2)==2);

% Check variable values
assert(isequal(round(projectedPoints,4),round(...
    [...
    0         0
    0         0
    1.0000         0
    0.7500    0.7500
    0.6667    1.0000
    0.6667    1.0000
    0.5000    1.0000
    0    1.0000
    ]...
    ,4)));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

vertices = [0 0; 1 0; 0.5 1.5; 1 1; 0 1; 0 0];
verticesIncludingSelfIntersections = fcn_MapGen_polytopeFindSelfIntersections(...
    vertices, -1);

interiorPoint = [0.5 0.5];

% Call the function
[projectedPoints] = ...
    fcn_MapGen_polytopeProjectVerticesOntoWalls(...,
    interiorPoint,...
    verticesIncludingSelfIntersections,...
    verticesIncludingSelfIntersections(1:end-1,:),...
    verticesIncludingSelfIntersections(2:end,:),...
    (-1));

% Check variable types
assert(isnumeric(projectedPoints));

% Check variable sizes
Nvertices = length(verticesIncludingSelfIntersections(:,1));
assert(size(projectedPoints,1)==Nvertices);
assert(size(projectedPoints,2)==2);

% Check variable values
assert(isequal(round(projectedPoints,4),round(...
    [...
    0         0
    0         0
    1.0000         0
    0.7500    0.7500
    0.6667    1.0000
    0.6667    1.0000
    0.5000    1.0000
    0    1.0000
    ]...
    ,4)));


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

vertices = [0 0; 1 0; 0.5 1.5; 1 1; 0 1; 0 0];
verticesIncludingSelfIntersections = fcn_MapGen_polytopeFindSelfIntersections(...
    vertices, -1);

interiorPoint = [0.5 0.5];

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [projectedPoints] = ...
        fcn_MapGen_polytopeProjectVerticesOntoWalls(...,
        interiorPoint,...
        verticesIncludingSelfIntersections,...
        verticesIncludingSelfIntersections(1:end-1,:),...
        verticesIncludingSelfIntersections(2:end,:),...
        ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [projectedPoints] = ...
        fcn_MapGen_polytopeProjectVerticesOntoWalls(...,
        interiorPoint,...
        verticesIncludingSelfIntersections,...
        verticesIncludingSelfIntersections(1:end-1,:),...
        verticesIncludingSelfIntersections(2:end,:),...
        (-1));
end
fast_method = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

% Plot results as bar chart
figure(373737);
clf;
hold on;

X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% BUG cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____  _    _  _____
% |  _ \| |  | |/ ____|
% | |_) | |  | | |  __    ___ __ _ ___  ___  ___
% |  _ <| |  | | | |_ |  / __/ _` / __|/ _ \/ __|
% | |_) | |__| | |__| | | (_| (_| \__ \  __/\__ \
% |____/ \____/ \_____|  \___\__,_|___/\___||___/
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=BUG%20cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All bug case figures start with the number 9

% close all;

%% BUG

%% Fail conditions
if 1==0
    %

end


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


% %% fcn_INTERNAL_loadExampleData
% function [seed_points, V, C] = fcn_INTERNAL_loadExampleData
%
%
% % pull halton set
% halton_points = haltonset(2);
% points_scrambled = scramble(halton_points,'RR2'); % scramble values
%
% % pick values from halton set
% Halton_range = [1801 1901];
% low_pt = Halton_range(1,1);
% high_pt = Halton_range(1,2);
% seed_points = points_scrambled(low_pt:high_pt,:);
% [V,C] = voronoin(seed_points);
% % V = V.*stretch;
% end % Ends fcn_INTERNAL_loadExampleData