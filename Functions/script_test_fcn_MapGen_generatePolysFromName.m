% script_test_fcn_MapGen_generatePolysFromName
% Tests function: fcn_MapGen_generatePolysFromName

% REVISION HISTORY:
% 
% 2021_06_06 by Sean Brennan, sbrennan@psu.edu
% - first written by S. Brennan.
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

%% DEMO case: basic demo
figNum = 10001;
titleString = sprintf('DEMO case: basic demo');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

map_name = "HST 1 100 SQT 0 1 0 1 SMV 0.01 0.001 1e-6 1111";
plot_flag = 1; 
disp_name = 0; 

line_style = 'r-';
line_width = 2;

% Call the function
[polytopes, h_fig] = fcn_MapGen_generatePolysFromName(map_name, plot_flag, disp_name, (figNum), (line_style), (line_width));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(polytopes));
assert(isfield(polytopes,'vertices'));
assert(isfield(polytopes,'xv'));
assert(isfield(polytopes,'yv'));
assert(isfield(polytopes,'distances'));
assert(isfield(polytopes,'mean'));
assert(isfield(polytopes,'area'));
assert(isfield(polytopes,'max_radius'));
assert(isfield(polytopes,'min_radius'));
assert(isfield(polytopes,'mean_radius'));
assert(isfield(polytopes,'radii'));
assert(isfield(polytopes,'cost'));
assert(isfield(polytopes,'parent_poly_id'));
assert(ishandle(h_fig));

% Check variable sizes
Npolys = 100;
assert(isequal(Npolys,length(polytopes))); 
assert(isequal(size(h_fig),[1 1]));

% Check variable values
assert(isequal(h_fig.Number,figNum));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% DEMO case: advanced demo
figNum = 10002;
titleString = sprintf('DEMO case: advanced demo');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

map_name = "HST 30 450 SQT 0 1 0 1 SMV 0.02 0.005 1e-6 1234";
plot_flag = 1; 
disp_name = [1, 0.05 -0.05, 0.5 0.5 0.5, 10];

line_style = '-'; 
line_width = 2; 
color = [0 0 1];
axis_limits = [0 1 -0.1 1]; 
axis_style = 'square';
fill_info = [1 1 0 1 0.5];

% Call the function
[polytopes, h_fig] = fcn_MapGen_generatePolysFromName(map_name, plot_flag, disp_name,...
    (figNum), (line_style), (line_width), (color), (axis_limits), (axis_style), (fill_info));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(polytopes));
assert(isfield(polytopes,'vertices'));
assert(isfield(polytopes,'xv'));
assert(isfield(polytopes,'yv'));
assert(isfield(polytopes,'distances'));
assert(isfield(polytopes,'mean'));
assert(isfield(polytopes,'area'));
assert(isfield(polytopes,'max_radius'));
assert(isfield(polytopes,'min_radius'));
assert(isfield(polytopes,'mean_radius'));
assert(isfield(polytopes,'radii'));
assert(isfield(polytopes,'cost'));
assert(isfield(polytopes,'parent_poly_id'));
assert(ishandle(h_fig));

% Check variable sizes
Npolys = 421;
assert(isequal(Npolys,length(polytopes))); 
assert(isequal(size(h_fig),[1 1]));

% Check variable values
assert(isequal(h_fig.Number,figNum));

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

map_name = "HST 1 100 SQT 0 1 0 1 SMV 0.01 0.001 1e-6 1111";
plot_flag = 1; 
disp_name = 0; 

line_style = 'r-';
line_width = 2;

% Call the function
[polytopes, h_fig] = fcn_MapGen_generatePolysFromName(map_name, plot_flag, disp_name, ([]), (line_style), (line_width));

% Check variable types
assert(isstruct(polytopes));
assert(isfield(polytopes,'vertices'));
assert(isfield(polytopes,'xv'));
assert(isfield(polytopes,'yv'));
assert(isfield(polytopes,'distances'));
assert(isfield(polytopes,'mean'));
assert(isfield(polytopes,'area'));
assert(isfield(polytopes,'max_radius'));
assert(isfield(polytopes,'min_radius'));
assert(isfield(polytopes,'mean_radius'));
assert(isfield(polytopes,'radii'));
assert(isfield(polytopes,'cost'));
assert(isfield(polytopes,'parent_poly_id'));
assert(isempty(h_fig));

% Check variable sizes
Npolys = 100;
assert(isequal(Npolys,length(polytopes))); 
% assert(isempty((size(h_fig),[1 1]));

% Check variable values
% assert(isequal(h_fig.Number,figNum));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: FAST mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

map_name = "HST 1 100 SQT 0 1 0 1 SMV 0.01 0.001 1e-6 1111";
plot_flag = 1; 
disp_name = 0; 

line_style = 'r-';
line_width = 2;

% Call the function
[polytopes, h_fig] = fcn_MapGen_generatePolysFromName(map_name, plot_flag, disp_name, (-1), (line_style), (line_width));

% Check variable types
assert(isstruct(polytopes));
assert(isfield(polytopes,'vertices'));
assert(isfield(polytopes,'xv'));
assert(isfield(polytopes,'yv'));
assert(isfield(polytopes,'distances'));
assert(isfield(polytopes,'mean'));
assert(isfield(polytopes,'area'));
assert(isfield(polytopes,'max_radius'));
assert(isfield(polytopes,'min_radius'));
assert(isfield(polytopes,'mean_radius'));
assert(isfield(polytopes,'radii'));
assert(isfield(polytopes,'cost'));
assert(isfield(polytopes,'parent_poly_id'));
assert(isempty(h_fig));

% Check variable sizes
Npolys = 100;
assert(isequal(Npolys,length(polytopes))); 
% assert(isempty((size(h_fig),[1 1]));

% Check variable values
% assert(isequal(h_fig.Number,figNum));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',figNum);
figure(figNum);
close(figNum);

map_name = "HST 1 100 SQT 0 1 0 1 SMV 0.01 0.001 1e-6 1111";
plot_flag = 1; 
disp_name = 0; 

line_style = 'r-';
line_width = 2;

Niterations = 10;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [polytopes, h_fig] = fcn_MapGen_generatePolysFromName(map_name, plot_flag, disp_name, ([]), (line_style), (line_width));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [polytopes, h_fig] = fcn_MapGen_generatePolysFromName(map_name, plot_flag, disp_name, (-1), (line_style), (line_width));
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