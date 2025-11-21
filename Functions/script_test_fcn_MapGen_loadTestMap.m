% script_test_fcn_MapGen_loadTestMap
% tests fcn_MapGen_loadTestMap

% REVISION HISTORY:
%
% As: fcn_BoundedAStar_loadTestMap
% 
% 2025_08_14 - K. Hayes, kxh1031@psu.edu
% - initial write of script
% 
% 2025_08_22 - K. Hayes
% - added test cases to check all maps
%
% As: fcn_MapGen_loadTestMap
% 
% 2025_10_31 by Sean Brennan, sbrennan@psu.edu
% - moved code into MapGen library and renamed it accordingly
% - cleaned up header docstrings to clearly name cases
% - changed map_idx to mapIndex
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

%% DEMO case: load a pre-generated polytope map (mapIndex = 9, Halton set)
figNum = 10001;
titleString = sprintf('DEMO case: load a pre-generated polytope map (mapIndex = 9, Halton set)');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

mapIndex = 9;
add_boundary = 0;

[polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = fcn_MapGen_loadTestMap(mapIndex, (add_boundary), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(starts));
assert(isstruct(polytopes));
assert(isnumeric(finishes));
assert(isnumeric(resolution_scale));
assert(isnumeric(length_cost_weights));
assert(isnumeric(navigated_portions));

% Check variable sizes
Npoly = 51;
assert(isequal(Npoly,length(polytopes))); 

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

%% TEST case: Load mapIndex = 1, generic canyon map
figNum = 20001;
titleString = sprintf('TEST case:  Load mapIndex = 1, generic canyon map');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

mapIndex = 1;
add_boundary = 0;

[polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = fcn_MapGen_loadTestMap(mapIndex, (add_boundary), (figNum));


sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(starts));
assert(isstruct(polytopes));
assert(isnumeric(finishes));
assert(isnumeric(resolution_scale));
assert(isnumeric(length_cost_weights));
assert(isnumeric(navigated_portions));

% Check variable sizes
Npoly = 43;
assert(isequal(Npoly,length(polytopes))); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: Load mapIndex = 2, the lower triangular flood plain
figNum = 20002;
titleString = sprintf('TEST case:  Load mapIndex = 2, the lower triangular flood plain');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

mapIndex = 2;
add_boundary = 0;

[polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = fcn_MapGen_loadTestMap(mapIndex, (add_boundary), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(starts));
assert(isstruct(polytopes));
assert(isnumeric(finishes));
assert(isnumeric(resolution_scale));
assert(isnumeric(length_cost_weights));
assert(isnumeric(navigated_portions));

% Check variable sizes
Npoly = 3;
assert(isequal(Npoly,length(polytopes))); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: Load mapIndex = 3, the mustafar mining rig map (the comb)
figNum = 20003;
titleString = sprintf('TEST case:  Load mapIndex = 3, the mustafar mining rig map (the comb)');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

mapIndex = 3;
add_boundary = 0;

[polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = fcn_MapGen_loadTestMap(mapIndex, (add_boundary), (figNum));


sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(starts));
assert(isstruct(polytopes));
assert(isnumeric(finishes));
assert(isnumeric(resolution_scale));
assert(isnumeric(length_cost_weights));
assert(isnumeric(navigated_portions));

% Check variable sizes
Npoly = 8;
assert(isequal(Npoly,length(polytopes))); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: Load mapIndex = 4, long river valleys
figNum = 20004;
titleString = sprintf('TEST case:  Load mapIndex = 4, long river valleys');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

mapIndex = 4;
add_boundary = 0;

[polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = fcn_MapGen_loadTestMap(mapIndex, (add_boundary), (figNum));


sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(starts));
assert(isstruct(polytopes));
assert(isnumeric(finishes));
assert(isnumeric(resolution_scale));
assert(isnumeric(length_cost_weights));
assert(isnumeric(navigated_portions));

% Check variable sizes
Npoly = 4;
assert(isequal(Npoly,length(polytopes))); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: Load mapIndex = 5, bridge map, good for random edge deletion case
figNum = 20005;
titleString = sprintf('TEST case:  Load mapIndex = 5, bridge map, good for random edge deletion case');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

mapIndex = 5;
add_boundary = 0;

[polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = fcn_MapGen_loadTestMap(mapIndex, (add_boundary), (figNum));


sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(starts));
assert(isstruct(polytopes));
assert(isnumeric(finishes));
assert(isnumeric(resolution_scale));
assert(isnumeric(length_cost_weights));
assert(isnumeric(navigated_portions));

% Check variable sizes
Npoly = 13;
assert(isequal(Npoly,length(polytopes))); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: Load mapIndex = 6, large map, good for dilation case, nearly fully tiled
figNum = 20006;
titleString = sprintf('TEST case:  Load mapIndex = 6, large map, good for dilation case, nearly fully tiled');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

mapIndex = 6;
add_boundary = 0;

[polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = fcn_MapGen_loadTestMap(mapIndex, (add_boundary), (figNum));


sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(starts));
assert(isstruct(polytopes));
assert(isnumeric(finishes));
assert(isnumeric(resolution_scale));
assert(isnumeric(length_cost_weights));
assert(isnumeric(navigated_portions));

% Check variable sizes
Npoly = 30;
assert(isequal(Npoly,length(polytopes))); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: Load mapIndex = 7, generic polytope map
figNum = 20007;
titleString = sprintf('TEST case:  Load mapIndex = 7, generic polytope map');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

mapIndex = 7;
add_boundary = 0;

[polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = fcn_MapGen_loadTestMap(mapIndex, (add_boundary), (figNum));


sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(starts));
assert(isstruct(polytopes));
assert(isnumeric(finishes));
assert(isnumeric(resolution_scale));
assert(isnumeric(length_cost_weights));
assert(isnumeric(navigated_portions));

% Check variable sizes
Npoly = 51;
assert(isequal(Npoly,length(polytopes))); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: Load mapIndex = 8, Josh's polytope map from 24 April 2024
figNum = 20008;
titleString = sprintf('TEST case:  Load mapIndex = 8, Josh''s polytope map from 24 April 2024');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

mapIndex = 8;
add_boundary = 0;

[polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = fcn_MapGen_loadTestMap(mapIndex, (add_boundary), (figNum));


sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(starts));
assert(isstruct(polytopes));
assert(isnumeric(finishes));
assert(isnumeric(resolution_scale));
assert(isnumeric(length_cost_weights));
assert(isnumeric(navigated_portions));

% Check variable sizes
Npoly = 25;
assert(isequal(Npoly,length(polytopes))); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%% TEST case: Load mapIndex = 9, Halton set
figNum = 20009;
titleString = sprintf('TEST case:  Load mapIndex = 9, Halton set');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

mapIndex = 9;
add_boundary = 0;

[polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = fcn_MapGen_loadTestMap(mapIndex, (add_boundary), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(starts));
assert(isstruct(polytopes));
assert(isnumeric(finishes));
assert(isnumeric(resolution_scale));
assert(isnumeric(length_cost_weights));
assert(isnumeric(navigated_portions));

% Check variable sizes
Npoly = 51;
assert(isequal(Npoly,length(polytopes))); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


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

mapIndex = 9;
add_boundary = 0;

[polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = ...
    fcn_MapGen_loadTestMap(mapIndex, (add_boundary), ([]));

% Check variable types
assert(isnumeric(starts));
assert(isstruct(polytopes));
assert(isnumeric(finishes));
assert(isnumeric(resolution_scale));
assert(isnumeric(length_cost_weights));
assert(isnumeric(navigated_portions));

% Check variable sizes
Npoly = 51;
assert(isequal(Npoly,length(polytopes))); 

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: FAST mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

mapIndex = 9;
add_boundary = 0;

[polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = ...
    fcn_MapGen_loadTestMap(mapIndex, (add_boundary), (-1));

% Check variable types
assert(isnumeric(starts));
assert(isstruct(polytopes));
assert(isnumeric(finishes));
assert(isnumeric(resolution_scale));
assert(isnumeric(length_cost_weights));
assert(isnumeric(navigated_portions));

% Check variable sizes
Npoly = 51;
assert(isequal(Npoly,length(polytopes))); 

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',figNum);
figure(figNum);
close(figNum);

mapIndex = 9;
add_boundary = 0;

Niterations = 50;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = ...
        fcn_MapGen_loadTestMap(mapIndex, (add_boundary), ([]));

end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [polytopes, starts, finishes, resolution_scale, length_cost_weights, navigated_portions] = ...
        fcn_MapGen_loadTestMap(mapIndex, (add_boundary), (-1));

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
        %% Fails because start_definition is not correct type
        
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