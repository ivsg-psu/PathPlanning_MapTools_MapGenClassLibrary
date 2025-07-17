% script_test_fcn_MapGen_verticesCropToAABB
% Tests function: fcn_MapGen_verticesCropToAABB

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

%% TEST case: simple square
fig_num = 10001;
titleString = sprintf('DEMO case: simple square');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

verticies = [1 1; 1 -1; -1 -1; -1 1; 1 1]*0.5;
interiorPoint = [0.25 0.25];
AABB = [0 0 1 1];

% Call the function
[croppedVertices] = fcn_MapGen_verticesCropToAABB(verticies, interiorPoint, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));

% Check variable sizes
Nvertices = length(verticies(:,1));
assert(size(croppedVertices,1)<=Nvertices);
assert(size(croppedVertices,2)==2);

% Check variable values
assert(isequal(round(croppedVertices,4),[0,0;0.5000,0;0.5000,0.5000;0,0.5000;0,0]));

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

%% TEST case: simple square, interiorPoint on top of AABB
fig_num = 20001;
titleString = sprintf('TEST case: simple square, query point on top of AABB');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

verticies = [1 1; 1 -1; -1 -1; -1 1; 1 1]*0.5;
interiorPoint = [0 0];
AABB = [0 0 1 1];

% Call the function
[croppedVertices] = fcn_MapGen_verticesCropToAABB(verticies, interiorPoint, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));

% Check variable sizes
Nvertices = length(verticies(:,1));
assert(size(croppedVertices,1)<=Nvertices);
assert(size(croppedVertices,2)==2);

% Check variable values
assert(isequal(round(croppedVertices,4),[0,0;0.5000,0;0.5000,0.5000;0,0.5000;0,0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% TEST case: simple triangle, interiorPoint on AABB
fig_num = 20002;
titleString = sprintf('TEST case: simple triangle, interiorPoint on AABB');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

verticies = [1 -0.5; -1 -1; -0.5 1; 1 -0.5];
interiorPoint = [0 0];
AABB = [0 0 1 1];


% Call the function
[croppedVertices] = fcn_MapGen_verticesCropToAABB(verticies, interiorPoint, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));

% Check variable sizes
Nvertices = length(verticies(:,1));
assert(size(croppedVertices,1)<=Nvertices);
assert(size(croppedVertices,2)==2);

% Check variable values
assert(isequal(round(croppedVertices,4),[         0         0;     0.5000         0;          0    0.5000;          0         0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% TEST case: vertices enclosing region
fig_num = 20003;
titleString = sprintf('TEST case: vertices enclosing region');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

verticies = [1 1; 1 -1; -1 -1; -1 1; 1 1]*2;
interiorPoint = [0 0];
AABB = [0 0 1 1];

% Call the function
[croppedVertices] = fcn_MapGen_verticesCropToAABB(verticies, interiorPoint, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));

% Check variable sizes
Nvertices = length(verticies(:,1));
assert(size(croppedVertices,1)<=Nvertices);
assert(size(croppedVertices,2)==2);

% Check variable values
assert(isequal(round(croppedVertices,4),[     0     0;      1     0;      1     1;      0     1;      0     0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: vertices all within region
fig_num = 20004;
titleString = sprintf('TEST case: vertices all within region');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

verticies = [ 0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25; 0.25 0.25;];
interiorPoint = [0.5 0.5];
AABB = [0 0 1 1];

% Call the function
[croppedVertices] = fcn_MapGen_verticesCropToAABB(verticies, interiorPoint, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));

% Check variable sizes
Nvertices = length(verticies(:,1));
assert(size(croppedVertices,1)<=Nvertices);
assert(size(croppedVertices,2)==2);

% Check variable values
assert(isequal(round(croppedVertices,4),[    0.2500    0.2500;     0.7500    0.2500;     0.7500    0.7500;     0.2500    0.7500;    0.2500    0.2500]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% TEST case: infinity in the numbers (crops to AABB)
fig_num = 20005;
titleString = sprintf('TEST case: infinity in the numbers (crops to AABB)');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

verticies = [...
   -0.0457    0.0471
   -0.4075    0.0851
   -2.7148    0.1857
  -11.3670    0.2644
       Inf       Inf
    0.3841  -53.6209
    0.2255   -6.8725
    0.1203   -1.0609
    0.0352   -0.0163
   -0.0081    0.0368
   -0.0457    0.0471];

interiorPoint = [0 0];
AABB = [0 0 1 1];

% Call the function
[croppedVertices] = fcn_MapGen_verticesCropToAABB(verticies, interiorPoint, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));

% Check variable sizes
Nvertices = length(verticies(:,1));
assert(size(croppedVertices,1)<=Nvertices);
assert(size(croppedVertices,2)==2);

% Check variable values
assert(isequal(round(croppedVertices,4),[    ...
    0         0
    0.0219         0
    0    0.0269
    0         0
    ]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% TEST case: Points on edges
fig_num = 20006;
titleString = sprintf('TEST case: Points on edges');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

verticies = [
   0.951225451552411   0.038481963258910
   0.951225451552411                   0
   1.000000000000000   0.055571820180916
   0.972711461970130   0.055571820180916
   0.967138332311524   0.053930059059772
   0.951225451552411   0.038481963258910];
interiorPoint = [0.978515625000000   0.035665294924554];
AABB = [0 0 1 1];

% Call the function
[croppedVertices] = fcn_MapGen_verticesCropToAABB(verticies, interiorPoint, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));

% Check variable sizes
Nvertices = length(verticies(:,1));
assert(size(croppedVertices,1)<=Nvertices);
assert(size(croppedVertices,2)==2);

% Check variable values
assert(isequal(round(croppedVertices,4),[     0.9512         0;    1.0000    0.0556;    0.9727    0.0556;    0.9671    0.0539;    0.9512    0.0385;    0.9512         0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% TEST case: Points valid except infinity
fig_num = 20007;
titleString = sprintf('TEST case: Points valid except infinity');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

interiorPoint = [0.9404    0.0133];

verticies = [
    0.9512    0.0385
       Inf       Inf
    0.9275    0.0315
    0.9318    0.0417
    0.9512    0.0385];

AABB = [0 0 1 1];

% Call the function
[croppedVertices] = fcn_MapGen_verticesCropToAABB(verticies, interiorPoint, AABB, (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(croppedVertices));

% Check variable sizes
Nvertices = length(verticies(:,1));
assert(size(croppedVertices,1)>=Nvertices); % Infinite values add vertices!
assert(size(croppedVertices,2)==2);

% Check variable values
assert(isequal(round(croppedVertices,4),[ ...
    0.9275         0
    0.9512         0
    0.9512    0.0385
    0.9318    0.0417
    0.9275    0.0315
    0.9275         0
    ]));

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

verticies = [1 1; 1 -1; -1 -1; -1 1; 1 1]*0.5;
interiorPoint = [0.25 0.25];
AABB = [0 0 1 1];

% Call the function
[croppedVertices] = fcn_MapGen_verticesCropToAABB(verticies, interiorPoint, AABB, ([]));

% Check variable types
assert(isnumeric(croppedVertices));

% Check variable sizes
Nvertices = length(verticies(:,1));
assert(size(croppedVertices,1)<=Nvertices);
assert(size(croppedVertices,2)==2);

% Check variable values
assert(isequal(round(croppedVertices,4),[0,0;0.5000,0;0.5000,0.5000;0,0.5000;0,0]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

verticies = [1 1; 1 -1; -1 -1; -1 1; 1 1]*0.5;
interiorPoint = [0.25 0.25];
AABB = [0 0 1 1];

% Call the function
[croppedVertices] = fcn_MapGen_verticesCropToAABB(verticies, interiorPoint, AABB, (-1));

% Check variable types
assert(isnumeric(croppedVertices));

% Check variable sizes
Nvertices = length(verticies(:,1));
assert(size(croppedVertices,1)<=Nvertices);
assert(size(croppedVertices,2)==2);

% Check variable values
assert(isequal(round(croppedVertices,4),[0,0;0.5000,0;0.5000,0.5000;0,0.5000;0,0]));



% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

verticies = [1 1; 1 -1; -1 -1; -1 1; 1 1]*0.5;
interiorPoint = [0.25 0.25];
AABB = [0 0 1 1];

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [croppedVertices] = fcn_MapGen_verticesCropToAABB(verticies, interiorPoint, AABB, ([]));

end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [croppedVertices] = fcn_MapGen_verticesCropToAABB(verticies, interiorPoint, AABB, (-1));
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