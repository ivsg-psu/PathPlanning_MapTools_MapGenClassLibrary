% script_test_fcn_MapGen_polytopesShrinkEvenly
% Tests: fcn_MapGen_polytopesShrinkEvenly

% REVISION HISTORY:
%
% 2025_07_29 by Sean Brennan, sbrennan@psu.edu
% - first write of script
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

%% DEMO case: shrinkage of a single polytope
figNum = 10001;
titleString = sprintf('DEMO case: shrinkage of a single polytope');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);
unshrunkPolytopeVerticesOnly.vertices = [
    1.0000    0.5217
    1.0000    0.5242
    0.9300    0.6329
    0.8472    0.6479
    0.8921    0.5627
    1.0000    0.5217
];
unshrunkPolytopes = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));

% Set parameters
cutDistance = 0.01; % Set the shrinkage distance

% Call the function
shrunkPolytopes = fcn_MapGen_polytopesShrinkEvenly(unshrunkPolytopes, cutDistance, figNum);

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(shrunkPolytopes));
assert(isfield(shrunkPolytopes,'vertices'));
assert(isfield(shrunkPolytopes,'xv'));
assert(isfield(shrunkPolytopes,'yv'));
assert(isfield(shrunkPolytopes,'distances'));
assert(isfield(shrunkPolytopes,'mean'));
assert(isfield(shrunkPolytopes,'area'));
assert(isfield(shrunkPolytopes,'max_radius'));
assert(isfield(shrunkPolytopes,'min_radius'));
assert(isfield(shrunkPolytopes,'mean_radius'));
assert(isfield(shrunkPolytopes,'radii'));
assert(isfield(shrunkPolytopes,'cost'));
assert(isfield(shrunkPolytopes,'parent_poly_id'));

% Check variable sizes
assert(isequal(length(unshrunkPolytopes),length(shrunkPolytopes))); 

% Check variable values
assert(shrunkPolytopes.area<unshrunkPolytopes.area); 
assert(all(shrunkPolytopes.max_radius<unshrunkPolytopes.max_radius)); 
assert(all(shrunkPolytopes.mean_radius<unshrunkPolytopes.mean_radius)); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% DEMO case: shrinkage of a polytope field
figNum = 10002;
titleString = sprintf('DEMO case: shrinkage of a polytope field');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

map_name = "HST 30 450 SQT 0 1 0 1 SMV 0.02 0.005 1e-6 1234";
plot_flag = 1; disp_name = [1, 0.05 -0.05, 0.5 0.5 0.5, 10];
line_style = '-'; line_width = 2; color = [0 0 1];
axis_limits = [0 1 -0.1 1]; axis_style = 'square';
fill_info = [1 1 0 1 0.5];
figNum = 7;

[unshrunkPolytopes,~] =fcn_MapGen_generatePolysFromName(...
    map_name,...
    plot_flag,...
    disp_name,...
    -1,...
    line_style,...
    line_width,....
    color,...
    axis_limits,...
    axis_style,...
    fill_info);

% Set shrinkage parameters
cutDistance = 0.01; % Set the shrinkage distance

% Call the function
shrunkPolytopes = fcn_MapGen_polytopesShrinkEvenly(unshrunkPolytopes, cutDistance, figNum);

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(shrunkPolytopes));
assert(isfield(shrunkPolytopes,'vertices'));
assert(isfield(shrunkPolytopes,'xv'));
assert(isfield(shrunkPolytopes,'yv'));
assert(isfield(shrunkPolytopes,'distances'));
assert(isfield(shrunkPolytopes,'mean'));
assert(isfield(shrunkPolytopes,'area'));
assert(isfield(shrunkPolytopes,'max_radius'));
assert(isfield(shrunkPolytopes,'min_radius'));
assert(isfield(shrunkPolytopes,'mean_radius'));
assert(isfield(shrunkPolytopes,'radii'));
assert(isfield(shrunkPolytopes,'cost'));
assert(isfield(shrunkPolytopes,'parent_poly_id'));

% Check variable sizes
assert(isequal(length(unshrunkPolytopes),length(shrunkPolytopes))); 

% Check variable values
% Too many values

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

unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);

unshrunkPolytopeVerticesOnly.vertices = [
    1.0000    0.5217
    1.0000    0.5242
    0.9300    0.6329
    0.8472    0.6479
    0.8921    0.5627
    1.0000    0.5217
];
unshrunkPolytopes = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));

% Set parameters

cutDistance = 0.01; % Set the shrinkage distance

% Call the function
shrunkPolytopes = fcn_MapGen_polytopesShrinkEvenly(unshrunkPolytopes, cutDistance, []);

% Check variable types
assert(isstruct(shrunkPolytopes));
assert(isfield(shrunkPolytopes,'vertices'));
assert(isfield(shrunkPolytopes,'xv'));
assert(isfield(shrunkPolytopes,'yv'));
assert(isfield(shrunkPolytopes,'distances'));
assert(isfield(shrunkPolytopes,'mean'));
assert(isfield(shrunkPolytopes,'area'));
assert(isfield(shrunkPolytopes,'max_radius'));
assert(isfield(shrunkPolytopes,'min_radius'));
assert(isfield(shrunkPolytopes,'mean_radius'));
assert(isfield(shrunkPolytopes,'radii'));
assert(isfield(shrunkPolytopes,'cost'));
assert(isfield(shrunkPolytopes,'parent_poly_id'));

% Check variable sizes
assert(isequal(length(unshrunkPolytopes),length(shrunkPolytopes))); 

% Check variable values
assert(shrunkPolytopes.area<unshrunkPolytopes.area); 
assert(all(shrunkPolytopes.max_radius<unshrunkPolytopes.max_radius)); 
assert(all(shrunkPolytopes.mean_radius<unshrunkPolytopes.mean_radius));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: FAST mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);

unshrunkPolytopeVerticesOnly.vertices = [
    1.0000    0.5217
    1.0000    0.5242
    0.9300    0.6329
    0.8472    0.6479
    0.8921    0.5627
    1.0000    0.5217
];
unshrunkPolytopes = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));


% Set parameters

cutDistance = 0.01; % Set the shrinkage distance

% Call the function
shrunkPolytopes = fcn_MapGen_polytopesShrinkEvenly(unshrunkPolytopes, cutDistance, -1);

% Check variable types
assert(isstruct(shrunkPolytopes));
assert(isfield(shrunkPolytopes,'vertices'));
assert(isfield(shrunkPolytopes,'xv'));
assert(isfield(shrunkPolytopes,'yv'));
assert(isfield(shrunkPolytopes,'distances'));
assert(isfield(shrunkPolytopes,'mean'));
assert(isfield(shrunkPolytopes,'area'));
assert(isfield(shrunkPolytopes,'max_radius'));
assert(isfield(shrunkPolytopes,'min_radius'));
assert(isfield(shrunkPolytopes,'mean_radius'));
assert(isfield(shrunkPolytopes,'radii'));
assert(isfield(shrunkPolytopes,'cost'));
assert(isfield(shrunkPolytopes,'parent_poly_id'));

% Check variable sizes
assert(isequal(length(unshrunkPolytopes),length(shrunkPolytopes))); 

% Check variable values
assert(shrunkPolytopes.area<unshrunkPolytopes.area); 
assert(all(shrunkPolytopes.max_radius<unshrunkPolytopes.max_radius)); 
assert(all(shrunkPolytopes.mean_radius<unshrunkPolytopes.mean_radius));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',figNum);
figure(figNum);
close(figNum);

unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);

unshrunkPolytopeVerticesOnly.vertices = [
    1.0000    0.5217
    1.0000    0.5242
    0.9300    0.6329
    0.8472    0.6479
    0.8921    0.5627
    1.0000    0.5217
];
unshrunkPolytopes = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));


% Set parameters

cutDistance = 0.01; % Set the shrinkage distance

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    shrunkPolytopes = fcn_MapGen_polytopesShrinkEvenly(unshrunkPolytopes, cutDistance, []);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    shrunkPolytopes = fcn_MapGen_polytopesShrinkEvenly(unshrunkPolytopes, cutDistance, []);
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