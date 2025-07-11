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