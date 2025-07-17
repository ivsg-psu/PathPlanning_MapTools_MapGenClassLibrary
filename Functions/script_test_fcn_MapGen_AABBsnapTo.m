% script_test_fcn_MapGen_AABBsnapTo
% Tests: fcn_MapGen_AABBsnapTo

% REVISION HISTORY:
% 2021_07_14 by Sean Brennan
% -- first write of script
% 2025_04_24 by Sean Brennan
% -- fixed calls to match revised function
% 2025_07_11 - S. Brennan, sbrennan@psu.edu
% -- updated script testing to standard form

close all;

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

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.3 0.2];
snapType = 0;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0.1667         0],round(snapPoint,4)))
assert(isequal(1,wallNumber))

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

%% Snap Type 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                     _______                    ___
%  / ____|                   |__   __|                  / _ \
% | (___  _ __   __ _ _ __      | |_   _ _ __   ___    | | | |
%  \___ \| '_ \ / _` | '_ \     | | | | | '_ \ / _ \   | | | |
%  ____) | | | | (_| | |_) |    | | |_| | |_) |  __/   | |_| |
% |_____/|_| |_|\__,_| .__/     |_|\__, | .__/ \___|    \___/
%                    | |            __/ | |
%                    |_|           |___/|_|

% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Snap%20Type%20%20%200
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All snap type 0 start with number 2

close all;
fprintf(1,'Figure: 2XXXXXX: TEST mode cases with Snap Type 0\n');

%% TEST case: inside, close to bottom
fig_num = 20001;
titleString = sprintf('TEST case: inside, close to bottom');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.3 0.2];
snapType = 0;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0.1667         0],round(snapPoint,4)))
assert(isequal(1,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% TEST case: inside, close to right
fig_num = 20002;
titleString = sprintf('TEST case: inside, close to right');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.8 0.4];
snapType = 0;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([1.0000    0.3333],round(snapPoint,4)))
assert(isequal(2,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: inside, close to top
fig_num = 20003;
titleString = sprintf('TEST case: inside, close to top');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.6 0.9];
snapType = 0;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0.6250    1.0000],round(snapPoint,4)))
assert(isequal(3,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: inside, close to left
fig_num = 20004;
titleString = sprintf('TEST case: inside, close to left');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.2 0.7];
snapType = 0;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0    0.8333],round(snapPoint,4)))
assert(isequal(4,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: outside, close to bottom
fig_num = 20005;
titleString = sprintf('TEST case: outside, close to bottom');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.3 -0.2];
snapType = 0;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0.3000   -0.2000],round(snapPoint,4)))
assert(isequal(1,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% check snap type 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                     _______                   __
%  / ____|                   |__   __|                 /_ |
% | (___  _ __   __ _ _ __      | |_   _ _ __   ___     | |
%  \___ \| '_ \ / _` | '_ \     | | | | | '_ \ / _ \    | |
%  ____) | | | | (_| | |_) |    | | |_| | |_) |  __/    | |
% |_____/|_| |_|\__,_| .__/     |_|\__, | .__/ \___|    |_|
%                    | |            __/ | |
%                    |_|           |___/|_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Snap%20Type%20%20%201
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All Snap Type 1 tests start with number 3
close all;
fprintf(1,'Figure: 3XXXXXX: TEST mode cases with Snap Type 1\n');

%% TEST case: Snap type 1, inside, close to bottom
fig_num = 30001;
titleString = sprintf('TEST case: Snap type 1, inside, close to bottom');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.3 0.2];
snapType = 1;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0.3         0],round(snapPoint,4)))
assert(isequal(1,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: Snap type 1, inside, close to right
fig_num = 30002;
titleString = sprintf('TEST case: Snap type 1, inside, close to right');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.8 0.4];
snapType = 1;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([1.0000    0.4],round(snapPoint,4)))
assert(isequal(2,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: Snap type 1, inside, close to top
fig_num = 30003;
titleString = sprintf('TEST case: Snap type 1, inside, close to top');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.6 0.9];
snapType = 1;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0.6    1.0000],round(snapPoint,4)))
assert(isequal(3,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% TEST case: Snap type 1, inside, close to left
fig_num = 30004;
titleString = sprintf('TEST case: Snap type 1, inside, close to left');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.2 0.7];
snapType = 1;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0    0.7],round(snapPoint,4)))
assert(isequal(4,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: Snap type 1, inside, close to bottom
fig_num = 30005;
titleString = sprintf('TEST case: Snap type 1, inside, close to bottom');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.3 -0.2];
snapType = 1;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0.3000   -0.2000],round(snapPoint,4)))
assert(isequal(1,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: Snap type 1, outside, close to right
fig_num = 30006;
titleString = sprintf('TEST case: Snap type 1, outside, close to right');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [1.2 0.4];
snapType = 1;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([1.2000    0.4000],round(snapPoint,4)))
assert(isequal(2,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: Snap type 1, outside, close to top
fig_num = 30006;
titleString = sprintf('TEST case: Snap type 1, outside, close to top');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.6 1.1];
snapType = 1;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0.6000    1.1000],round(snapPoint,4)))
assert(isequal(3,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: Snap type 1, outside, close to top
fig_num = 30007;
titleString = sprintf('TEST case: Snap type 1, outside, close to top');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.6 1.1];
snapType = 1;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0.6000    1.1000],round(snapPoint,4)))
assert(isequal(3,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: Snap type 1, outside, close to left
fig_num = 30008;
titleString = sprintf('TEST case: Snap type 1, outside, close to left');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [-0.2 0.7];
snapType = 1;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([-0.2000    0.7000],round(snapPoint,4)))
assert(isequal(4,wallNumber))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% check snap type 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                     _______                   ___
%  / ____|                   |__   __|                 |__ \
% | (___  _ __   __ _ _ __      | |_   _ _ __   ___       ) |
%  \___ \| '_ \ / _` | '_ \     | | | | | '_ \ / _ \     / /
%  ____) | | | | (_| | |_) |    | | |_| | |_) |  __/    / /_
% |_____/|_| |_|\__,_| .__/     |_|\__, | .__/ \___|   |____|
%                    | |            __/ | |
%                    |_|           |___/|_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Snap%20Type%20%20%202
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All snap type 2's start with the number 4
close all;
fprintf(1,'Figure: 4XXXXXX: TEST mode cases with Snap Type 2\n');

%% TEST case: Snap type 2, user-defined vector
fig_num = 40001;
titleString = sprintf('TEST case: Snap type 2, user-defined vector');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.2 0.7; 0.2 0.9];
snapType = 2;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (fig_num) );

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0.2000    1.0000],round(snapPoint,4)))
assert(isequal(3,wallNumber))

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

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.3 0.2];
snapType = 0;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), ([]) );

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0.1667         0],round(snapPoint,4)))
assert(isequal(1,wallNumber))

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.3 0.2];
snapType = 0;

% Call the function
[snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (-1) );

% Check variable types
assert(isnumeric(snapPoint));
assert(isnumeric(wallNumber));

% Check variable sizes
assert(isequal(size(snapPoint),[1 2]));
assert(isequal(size(wallNumber),[1 1]));

% Check variable values
assert(isequal([0.1667         0],round(snapPoint,4)))
assert(isequal(1,wallNumber))

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

axisAlignedBoundingBox = [0 0 1 1];
testPoint = [0.3 0.2];
snapType = 0;

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), ([]) );
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [snapPoint, wallNumber] = fcn_MapGen_AABBsnapTo(axisAlignedBoundingBox, testPoint, (snapType), (-1) );
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