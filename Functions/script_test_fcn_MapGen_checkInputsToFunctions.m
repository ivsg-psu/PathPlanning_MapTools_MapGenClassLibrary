% script_test_fcn_MapGen_checkInputsToFunctions.m
% Tests fcn_MapGen_checkInputsToFunctions

% Revision history:
%      2021_06_06:
%      -- first write of the code copying functionality from fcn_FastestTraversal_checkInputsToFunctions

%% column_of_numbers
%             _                               __                       _
%            | |                             / _|                     | |
%    ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___
%   / __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
%  | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
%   \___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
%                                  ______   ______
%                                 |______| |______|
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs

%% Test the column_of_numbers type (success)
column_of_numbers_test = 4;
fcn_MapGen_checkInputsToFunctions(column_of_numbers_test, 'column_of_numbers');

column_of_numbers_test = [4; 3; 2];
fcn_MapGen_checkInputsToFunctions(column_of_numbers_test, 'column_of_numbers');

column_of_numbers_test = [4; 3; 2];
fcn_MapGen_checkInputsToFunctions(column_of_numbers_test, 'column_of_numbers',3);

%% 2column_of_numbers
% 
%   ___           _                               __                       _                   
%  |__ \         | |                             / _|                     | |                  
%     ) |___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___ 
%    / // __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
%   / /| (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
%  |____\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
%                                      ______   ______                                         
%                                     |______| |______|                                        
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs

%% Test the column_of_numbers type (success)
Twocolumn_of_numbers_test = [4 2];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers');

Twocolumn_of_numbers_test = [4 1; 3 0; 2 5];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers');

Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',3);

% Minimum length is 2 or greater
Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[2 3]);

% Maximum length is 5 or less
Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[5 4]);

% Maximum length is 3
Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[3 3]);


%% 2or3column_of_numbers
% 
%   ___           ____            _                               __                       _                   
%  |__ \         |___ \          | |                             / _|                     | |                  
%     ) |___  _ __ __) | ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___ 
%    / // _ \| '__|__ < / __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
%   / /| (_) | |  ___) | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
%  |____\___/|_| |____/ \___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
%                                                      ______   ______                                         
%                                                     |______| |______|                                        
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs


%% Test the column_of_numbers type (success)
% Test 1 by 2
TwoOrThreeColumn_of_numbers_test = [4 2];
fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

% Test 1 by 3
TwoOrThreeColumn_of_numbers_test = [4 2 1];
fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

% Test multiple points - 2 columns
TwoOrThreeColumn_of_numbers_test = [4 1; 3 0; 2 5];
fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

% Test multiple points - 3 columns
TwoOrThreeColumn_of_numbers_test = [4 1 3; 3 0 5; 2 5 7];
fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

% Test specified length - 2 columns
TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',3);

% Test specified length - 3 columns
TwoOrThreeColumn_of_numbers_test = [4 1 5; 3 9 5; 2 7 5];
fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',3);

% Minimum length is 2 or greater - 2 columns
TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[2 3]);

% Minimum length is 2 or greater - 3 columns
TwoOrThreeColumn_of_numbers_test = [4 1 5; 3 9 5; 2 7 5];
fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[2 3]);

% Maximum length is 5 or less - 2 column
TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[5 4]);

% Maximum length is 5 or less - 3 column
TwoOrThreeColumn_of_numbers_test = [4 1 4; 3 9 4; 2 7 4];
fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[5 4]);

% Length MUST be 3 - 2 column
TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[3 3]);

% Length MUST be 3 - 3 column
TwoOrThreeColumn_of_numbers_test = [4 1 3; 3 9 3; 2 7 3];
fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[3 3]);

%% 2column_of_integers
% 
% 
%   ___           _                               __ _       _                           
%  |__ \         | |                             / _(_)     | |                          
%     ) |___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ _ __ | |_ ___  __ _  ___ _ __ ___ 
%    / // __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| | '_ \| __/ _ \/ _` |/ _ \ '__/ __|
%   / /| (_| (_) | | |_| | | | | | | | | || (_) | | | | | | | ||  __/ (_| |  __/ |  \__ \
%  |____\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_|_| |_|\__\___|\__, |\___|_|  |___/
%                                      ______   ______                __/ |              
%                                     |______| |______|              |___/               
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs

%% Test the 2column_of_integers type (success)
Twocolumn_of_integers_test = [4 2];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers');

Twocolumn_of_integers_test = [4 1; 3 0; 2 5];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers');

Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',3);

% Minimum length is 2 or greater
Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',[2 3]);

% Maximum length is 5 or less
Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',[5 4]);

% Maximum length is 3
Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',[3 3]);


%% Fail conditions
if 1==0
    %% column_of_numbers
    %             _                               __                       _
    %            | |                             / _|                     | |
    %    ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___
    %   / __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
    %  | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
    %   \___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
    %                                  ______   ______
    %                                 |______| |______|
    %
    %
    
    
    %% Test the column_of_numbers type (FAILURE because 1 x 2)
    column_of_numbers_test = [4 1];
    fcn_MapGen_checkInputsToFunctions(column_of_numbers_test, 'column_of_numbers');
    
    %% Test the column_of_numbers type (FAILURE because 3 long, not 2)
    column_of_numbers_test = [4; 3; 2];
    fcn_MapGen_checkInputsToFunctions(column_of_numbers_test, 'column_of_numbers',2);
       
    %% Test the column_of_numbers type (FAILURE because NaN)
    column_of_numbers_test = [4; nan; 2];
    fcn_MapGen_checkInputsToFunctions(column_of_numbers_test, 'column_of_numbers',3);
    
    %% 2column_of_numbers
    %
    %   ___           _                               __                       _
    %  |__ \         | |                             / _|                     | |
    %     ) |___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___
    %    / // __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
    %   / /| (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
    %  |____\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
    %                                      ______   ______
    %                                     |______| |______|
    %
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
    
    %% Test the column_of_numbers type (FAILURE because 1 x 1)
    Twocolumn_of_numbers_test = [4];
    fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers');
    
    %% Test the column_of_numbers type (FAILURE because 3 long, not 2)
    Twocolumn_of_numbers_test = [4 1; 3 1; 2 2];
    fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',2);
       
    %% Test the column_of_numbers type (FAILURE because NaN)
    Twocolumn_of_numbers_test = [4 1; nan 1; 2 0];
    fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',3);
    
    
    %% Minimum length is 4 or greater
    Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[4 5]);
    
    %% Maximum length is 2 or less
    Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[2 1]);
    
    %% Maximum length is 2
    Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[2 2]);
    
    %% 2or3column_of_numbers
    %
    %   ___           ____            _                               __                       _
    %  |__ \         |___ \          | |                             / _|                     | |
    %     ) |___  _ __ __) | ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___
    %    / // _ \| '__|__ < / __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
    %   / /| (_) | |  ___) | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
    %  |____\___/|_| |____/ \___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
    %                                                      ______   ______
    %                                                     |______| |______|
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
    
    

    %% Test the column_of_numbers type (FAILURE because 1 x 1)
    TwoOrThreeColumn_of_numbers_test = [4];
    fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

    %% Test the column_of_numbers type (FAILURE because 1 x 4)
    TwoOrThreeColumn_of_numbers_test = [4 1 1 1];
    fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

    %% Test the column_of_numbers type (FAILURE because 3 long, not 2)
    TwoOrThreeColumn_of_numbers_test = [4 1; 3 1; 2 2];
    fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',2);
       
    %% Test the column_of_numbers type (FAILURE because NaN)
    TwoOrThreeColumn_of_numbers_test = [4 1; nan 1; 2 0];
    fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',3);
    
    
    %% Minimum length is 4 or greater
    TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[4 5]);
    
    %% Maximum length is 2 or less
    TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[2 1]);
    
    %% Maximum length is 2
    TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_MapGen_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[2 2]);
    
    %% 2column_of_integers
    %
    %
    %   ___           _                               __ _       _
    %  |__ \         | |                             / _(_)     | |
    %     ) |___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ _ __ | |_ ___  __ _  ___ _ __ ___
    %    / // __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| | '_ \| __/ _ \/ _` |/ _ \ '__/ __|
    %   / /| (_| (_) | | |_| | | | | | | | | || (_) | | | | | | | ||  __/ (_| |  __/ |  \__ \
    %  |____\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_|_| |_|\__\___|\__, |\___|_|  |___/
    %                                      ______   ______                __/ |
    %                                     |______| |______|              |___/
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
    %% Test the 2column_of_integers type (FAILURE because not integers)
    Twocolumn_of_integers_test = [4 3.2];
    fcn_MapGen_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers');
    
    %% Test the 2column_of_integers type (FAILURE because 1 x 1)
    Twocolumn_of_integers_test = [4];
    fcn_MapGen_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers');
    
    %% Test the 2column_of_integers type (FAILURE because 3 long, not 2)
    Twocolumn_of_integers_test = [4 1; 3 1; 2 2];
    fcn_MapGen_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',2);
       
    %% Test the 2column_of_integers type (FAILURE because NaN)
    Twocolumn_of_integers_test = [4 1; nan 1; 2 0];
    fcn_MapGen_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',3);
    
    
    %% Minimum length is 4 or greater
    Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
    fcn_MapGen_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',[4 5]);
    
    %% Maximum length is 2 or less
    Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
    fcn_MapGen_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',[2 1]);
    
    %% Maximum length is 2
    Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
    fcn_MapGen_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',[2 2]);

end