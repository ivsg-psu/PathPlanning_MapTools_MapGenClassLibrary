% script_test_fcn_MapGen_polytopesIncreaseVertexCount
% Tests: fcn_MapGen_polytopesIncreaseVertexCount

% REVISION HISTORY:
%
% 2022_10_21 by Steve Harnett
% - first write of script
% 
% 2025_07_11 by Sean Brennan, sbrennan@psu.edu
% - updated script testing to standard form
% 
% 2025_11_20 by Sean Brennan, sbrennan@psu.edu
% - Updated rev history to be in Markdown format
% - Replaced fig_+num with figNum

% TO-DO:
% 
% 2025_11_20 by Sean Brennan, sbrennan@psu.edu
% - fill in to-do items here.


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


%% DEMO case: basic demo, one polytope set
figNum = 10001;
titleString = sprintf('DEMO case: basic demo, one polytope set');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Load test data
polytopes = fcn_INTERNAL_loadExampleData;

resolution = 0.05;

% Call the function
interpolatedPolytopes = fcn_MapGen_polytopesIncreaseVertexCount(polytopes,resolution/2, figNum);

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(interpolatedPolytopes));
assert(isfield(interpolatedPolytopes,'vertices'));
assert(isfield(interpolatedPolytopes,'xv'));
assert(isfield(interpolatedPolytopes,'yv'));
assert(isfield(interpolatedPolytopes,'distances'));
assert(isfield(interpolatedPolytopes,'mean'));
assert(isfield(interpolatedPolytopes,'area'));
assert(isfield(interpolatedPolytopes,'max_radius'));
assert(isfield(interpolatedPolytopes,'min_radius'));
assert(isfield(interpolatedPolytopes,'mean_radius'));
assert(isfield(interpolatedPolytopes,'radii'));
assert(isfield(interpolatedPolytopes,'cost'));
assert(isfield(interpolatedPolytopes,'parent_poly_id'));

% Check variable sizes
assert(isequal(length(polytopes),length(interpolatedPolytopes))); 

% Check variable values
% (can't - randomly generated!)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

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
% figNum = 20001;
% titleString = sprintf('TEST case: simple crossing at origin');
% fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
% figure(figNum); clf;


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
figNum = 80001;
fprintf(1,'Figure: %.0f: FAST mode, empty figNum\n',figNum);
figure(figNum); close(figNum);

% Load test data
polytopes = fcn_INTERNAL_loadExampleData;

resolution = 0.05;

% Call the function
interpolatedPolytopes = fcn_MapGen_polytopesIncreaseVertexCount(polytopes,resolution/2, ([]));

% Check variable types
assert(isstruct(interpolatedPolytopes));
assert(isfield(interpolatedPolytopes,'vertices'));
assert(isfield(interpolatedPolytopes,'xv'));
assert(isfield(interpolatedPolytopes,'yv'));
assert(isfield(interpolatedPolytopes,'distances'));
assert(isfield(interpolatedPolytopes,'mean'));
assert(isfield(interpolatedPolytopes,'area'));
assert(isfield(interpolatedPolytopes,'max_radius'));
assert(isfield(interpolatedPolytopes,'min_radius'));
assert(isfield(interpolatedPolytopes,'mean_radius'));
assert(isfield(interpolatedPolytopes,'radii'));
assert(isfield(interpolatedPolytopes,'cost'));
assert(isfield(interpolatedPolytopes,'parent_poly_id'));

% Check variable sizes
assert(isequal(length(polytopes),length(interpolatedPolytopes))); 

% Check variable values
% (can't - randomly generated!)

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: FAST mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

% Load test data
polytopes = fcn_INTERNAL_loadExampleData;

resolution = 0.05;

% Call the function
interpolatedPolytopes = fcn_MapGen_polytopesIncreaseVertexCount(polytopes,resolution/2, (-1));

% Check variable types
assert(isstruct(interpolatedPolytopes));
assert(isfield(interpolatedPolytopes,'vertices'));
assert(isfield(interpolatedPolytopes,'xv'));
assert(isfield(interpolatedPolytopes,'yv'));
assert(isfield(interpolatedPolytopes,'distances'));
assert(isfield(interpolatedPolytopes,'mean'));
assert(isfield(interpolatedPolytopes,'area'));
assert(isfield(interpolatedPolytopes,'max_radius'));
assert(isfield(interpolatedPolytopes,'min_radius'));
assert(isfield(interpolatedPolytopes,'mean_radius'));
assert(isfield(interpolatedPolytopes,'radii'));
assert(isfield(interpolatedPolytopes,'cost'));
assert(isfield(interpolatedPolytopes,'parent_poly_id'));

% Check variable sizes
assert(isequal(length(polytopes),length(interpolatedPolytopes))); 

% Check variable values
% (can't - randomly generated!)

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',figNum);
figure(figNum);
close(figNum);

% Load test data
polytopes = fcn_INTERNAL_loadExampleData;

resolution = 0.05;

Niterations = 10;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    interpolatedPolytopes = fcn_MapGen_polytopesIncreaseVertexCount(polytopes,resolution/2, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    interpolatedPolytopes = fcn_MapGen_polytopesIncreaseVertexCount(polytopes,resolution/2, (-1));
end
fast_method = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));

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
assert(~any(figHandles==figNum));


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


%% fcn_INTERNAL_loadExampleData
function polytopes = fcn_INTERNAL_loadExampleData
% axis_box = [0 1 0 1];

seedGeneratorNames = 'haltonset';
seedGeneratorRanges = [1 100];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[polytopes] = fcn_MapGen_generatePolysFromSeedGeneratorNames(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (-1));

end % Ends fcn_INTERNAL_loadExampleData