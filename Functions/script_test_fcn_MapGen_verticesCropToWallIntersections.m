% script_test_fcn_MapGen_verticesCropToWallIntersections
% Tests: fcn_MapGen_verticesCropToWallIntersections

% 
% REVISION HISTORY:
% 
% 2021_07_15 by Sean Brennan
% -- first write of script
% -- remove dependence on test fixture
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

%% DEMO case: basic example with 2 intersections
fig_num = 10001;
titleString = sprintf('DEMO case: basic example with 2 intersections');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

walls = [0 0; 2 0; 1 3; 0 0];
inputVertices = [1 -1; 1 1; 2 2; 2 2.5; -1 2.5; -1 -1; 1 -1];

% Call the function
[croppedVertices, NwallsHit] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));
assert(isnumeric(NwallsHit));

% Check variable sizes
Nvertices = 5;
assert(isequal(size(croppedVertices),[Nvertices 2]));
assert(isequal(size(NwallsHit),[1 1]));

% Check variable values
assert(isequal(round(croppedVertices,4),round(...
    [...
    1.0000         0
    1.0000    1.0000
    1.5000    1.5000
    1.1667    2.5000
    0.8333    2.5000
    ]...
    ,4)));
assert(isequal(NwallsHit,3));

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

%% TEST case: Going from out to in
fig_num = 20001;
titleString = sprintf('TEST case: Going from out to in');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

walls = [0 0; 1 0; 1 1; 0 1; 0 0];
inputVertices = [-0.3 0.2; 0.3 0.2];

% Call the function
[croppedVertices, NwallsHit] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));
assert(isnumeric(NwallsHit));

% Check variable sizes
Nvertices = 2;
assert(isequal(size(croppedVertices),[Nvertices 2]));
assert(isequal(size(NwallsHit),[1 1]));

% Check variable values
assert(all(([0 0.2; 0.3 0.2]-eps*[1 1; 1 1])<croppedVertices,'all') && all([0 0.2; 0.3 0.2]+eps*[1 1; 1 1]>croppedVertices,'all'));
assert(isequal(NwallsHit,1))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: Going from in to out
fig_num = 20002;
titleString = sprintf('TEST case: Going from in to out');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

walls = [0 0; 1 0; 1 1; 0 1; 0 0];
inputVertices = [0.3 0.2; 1.3 0.2];

% Call the function
[croppedVertices, NwallsHit] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));
assert(isnumeric(NwallsHit));

% Check variable sizes
Nvertices = 2;
assert(isequal(size(croppedVertices),[Nvertices 2]));
assert(isequal(size(NwallsHit),[1 1]));

% Check variable values
assert(all(([0.3 0.2; 1 0.2]-eps*[1 1; 1 1])<croppedVertices,'all') && all([0.3 0.2; 1 0.2]+eps*[1 1; 1 1]>croppedVertices,'all'));
assert(isequal(NwallsHit,1))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: all inside, nothing to crop
fig_num = 20003;
titleString = sprintf('TEST case: all inside, nothing to crop');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

walls = [0 0; 1 0; 1 1; 0 1; 0 0];
inputVertices = [0.3 0.2; 0.4 0.2];

% Call the function
[croppedVertices, NwallsHit] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));
assert(isnumeric(NwallsHit));

% Check variable sizes
Nvertices = 2;
assert(isequal(size(croppedVertices),[Nvertices 2]));
assert(isequal(size(NwallsHit),[1 1]));

% Check variable values
assert(all(([0.3 0.2; 0.4 0.2]-eps*[1 1; 1 1])<croppedVertices,'all') && all([0.3 0.2; 0.4 0.2]+eps*[1 1; 1 1]>croppedVertices,'all'));
assert(isequal(NwallsHit,0))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: all outside, returns empty
fig_num = 20004;
titleString = sprintf('TEST case: all outside, returns empty');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

walls = [0 0; 1 0; 1 1; 0 1; 0 0];
inputVertices = [-0.3 0.2; -0.4 0.2];

% Call the function
[croppedVertices, NwallsHit] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));
assert(isnumeric(NwallsHit));

% Check variable sizes
assert(isempty(croppedVertices));
assert(isequal(size(NwallsHit),[1 1]));

% Check variable values
assert(isequal([],croppedVertices));
assert(isequal(NwallsHit,0))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: crossing over, crops both sides keeping inside
fig_num = 20005;
titleString = sprintf('TEST case: crossing over, crops both sides keeping inside');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

walls = [0 0; 1 0; 1 1; 0 1; 0 0];
inputVertices = [-0.3 0.2; 1.4 0.2];

% Call the function
[croppedVertices, NwallsHit] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));
assert(isnumeric(NwallsHit));

% Check variable sizes
Nvertices = 2;
assert(isequal(size(croppedVertices),[Nvertices 2]));
assert(isequal(size(NwallsHit),[1 1]));

% Check variable values
true_answer = [1 0.2; 0 0.2];
assert(all((true_answer-eps*[1 1; 1 1])<croppedVertices,'all') && all(true_answer+eps*[1 1; 1 1]>croppedVertices,'all'));
assert(isequal(NwallsHit,2))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: aligned with edge, across, crops both sides keeping inside
fig_num = 20006;
titleString = sprintf('TEST case: aligned with edge, across, crops both sides keeping inside');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

walls = [0 0; 1 0; 1 1; 0 1; 0 0];
inputVertices = [-0.3 0; 1.4 0];

% Call the function
[croppedVertices, NwallsHit] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));
assert(isnumeric(NwallsHit));

% Check variable sizes
Nvertices = 2;
assert(isequal(size(croppedVertices),[Nvertices 2]));
assert(isequal(size(NwallsHit),[1 1]));

% Check variable values
true_answer = [0 0; 1 0];
assert(all((true_answer-eps*[1 1; 1 1])<croppedVertices,'all') && all(true_answer+eps*[1 1; 1 1]>croppedVertices,'all'));
assert(isequal(NwallsHit,3))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: aligned with edge, inside to outside, crops outside but keeps edge
fig_num = 20007;
titleString = sprintf('TEST case: aligned with edge, inside to outside, crops outside but keeps edge');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

walls = [0 0; 1 0; 1 1; 0 1; 0 0];
inputVertices = [0.5 0; 1.4 0];

% Call the function
[croppedVertices, NwallsHit] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));
assert(isnumeric(NwallsHit));

% Check variable sizes
Nvertices = 2;
assert(isequal(size(croppedVertices),[Nvertices 2]));
assert(isequal(size(NwallsHit),[1 1]));

% Check variable values
true_answer = [0.5 0; 1 0];
assert(all((true_answer-eps*[1 1; 1 1])<croppedVertices,'all') && all(true_answer+eps*[1 1; 1 1]>croppedVertices,'all'));
assert(isequal(NwallsHit,2))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: aligned with edge, on corner to outside, keeps only corner
fig_num = 20008;
titleString = sprintf('TEST case: aligned with edge, on corner to outside, keeps only corner');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

walls = [0 0; 1 0; 1 1; 0 1; 0 0];
inputVertices = [1 0; 1.4 0];

% Call the function
[croppedVertices, NwallsHit] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));
assert(isnumeric(NwallsHit));

% Check variable sizes
Nvertices = 1;
assert(isequal(size(croppedVertices),[Nvertices 2]));
assert(isequal(size(NwallsHit),[1 1]));

% Check variable values
true_answer = [1 0];
assert(all((true_answer-eps*[1 1; 1 1])<croppedVertices,'all') && all(true_answer+eps*[1 1; 1 1]>croppedVertices,'all'));
assert(isequal(NwallsHit,2))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



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

walls = [0 0; 2 0; 1 3; 0 0];
inputVertices = [1 -1; 1 1; 2 2; 2 2.5; -1 2.5; -1 -1; 1 -1];

% Call the function
[croppedVertices, NwallsHit] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, ([]));

% Check variable types
assert(isnumeric(croppedVertices));
assert(isnumeric(NwallsHit));

% Check variable sizes
Nvertices = 5;
assert(isequal(size(croppedVertices),[Nvertices 2]));
assert(isequal(size(NwallsHit),[1 1]));

% Check variable values
assert(isequal(round(croppedVertices,4),round(...
    [...
    1.0000         0
    1.0000    1.0000
    1.5000    1.5000
    1.1667    2.5000
    0.8333    2.5000
    ]...
    ,4)));
assert(isequal(NwallsHit,3));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

walls = [0 0; 2 0; 1 3; 0 0];
inputVertices = [1 -1; 1 1; 2 2; 2 2.5; -1 2.5; -1 -1; 1 -1];

% Call the function
[croppedVertices, NwallsHit] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, (-1));

% Check variable types
assert(isnumeric(croppedVertices));
assert(isnumeric(NwallsHit));

% Check variable sizes
Nvertices = 5;
assert(isequal(size(croppedVertices),[Nvertices 2]));
assert(isequal(size(NwallsHit),[1 1]));

% Check variable values
assert(isequal(round(croppedVertices,4),round(...
    [...
    1.0000         0
    1.0000    1.0000
    1.5000    1.5000
    1.1667    2.5000
    0.8333    2.5000
    ]...
    ,4)));
assert(isequal(NwallsHit,3));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

walls = [0 0; 2 0; 1 3; 0 0];
inputVertices = [1 -1; 1 1; 2 2; 2 2.5; -1 2.5; -1 -1; 1 -1];

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [croppedVertices, NwallsHit] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [croppedVertices, NwallsHit] = fcn_MapGen_verticesCropToWallIntersections(inputVertices, walls, (-1));
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