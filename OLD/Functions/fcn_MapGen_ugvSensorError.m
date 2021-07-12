function [ ...
DX_err, ...
DY_err, ...
DZ_err ...
] = ...
fcn_MapGen_ugvSensorError( ...
Scanning_Results, ...
Position_Uncertainty, ...
Angular_Uncertainty, ...
Laser_Uncertainty, ...
varargin...
)
% fcn_MapGen_ugvSensorError
% calculates error in sensed locations from a UGV perspective
% 
%    Detailed description is as follows:    - scanning_results is a Nx3 
% matrix, N being the number of received    lidar points or sensor scans   
%      -- entry 1: R (range or distance of laser point)        -- entry 2: 
% beta (scanning angle of respective laser point)        -- entry 3: kappa 
% (heading of vehicle wrt Euclidian line)         - position_uncertainty 
% is a vector of constants        -- entry 1: dx0 (uncertainty of sensor 
% location in x direction)        -- entry 2: dy0 (uncertainty of sensor 
% location in y direction)        -- entry 3: dz0 (uncertainty of sensor 
% location in z direction)         -angular_uncertainty is a vector of 
% constants        -- entry 1: domega (uncertainty of sensor pointing 
% angle around the        x axis)        -- entry 1: dpsi (uncertainty of 
% sensor pointing angle around the        y axis)        -- entry 1: 
% dkappa (uncertainty of sensor pointing angle around the        z axis)   
%       -laser_uncertainty is a vector of constants        -- entry 1: 
% dbeta (uncertainty of scanning angle around x axis)        -- entry 2: 
% dR (uncertainty of scanning range value)     Enter in units of meters 
% and degrees.     This function is designed to output a 'bubble' or 
% 'shadow' of uncertainty that pertains to how well a sensor on a UGV 
% perceives the objects within its scanning range. The inputs to this 
% function are primarily error values specified by the sensor itself, as 
% well as the distance and scanning angles to the object or point (e.g. 
% for lidar, consider one reflection from an emitted pulse at a given 
% scanning angle). Error values will be defined on a global (map) 
% coordinate system, not the local (sensor/vehicle) system.     - x,y,z : 
% local coordinate system of sensor, beam origin [meters]    - x : 
% direction of forward travel    - y : direction of sideways travel (i.e. 
% left, right)    - z : altitude or distance above ground - X,Y,Z : global 
% coordinate system of ground or targeted object, nadir [meters]    - X : 
% intended path of travel (parallel to Euclidian distance line)    - Y : 
% deviations from intended path occur along Y (for flat ground)    - Z : 
% global reference to ground for altitude/height - omega : roll, rotation 
% around x axis [degrees] - psi : pitch, rotation around y axis [degrees] 
% - kappa : yaw, rotation around z axis, angle from Euclidian line to 
% actual    line of travel, should include the heading uncertainty 
% [degrees] - beta : scan angle from x0,y0,z0 to X0,Y0,Z0 [degrees] - R : 
% range or distance from perceived object [meters]     Assumptions:    - 
% flat terrain    - rake scanning across y axis (or y-z plane but not 
% considering    elevation or angle), path of laser is parallel to y-x 
% plane (flat)    - assumes that scanning frequency (scanning paths 
% completed per second)    is much greater than vehicle velocity 
% (essentially, the vehicle moves a    negligible distance between each 
% return of the laser within a scan)    - ignoring other error sources: 
% time offset (scanner to clock),    calibration offset or misalignment 
% between sensors, possible errors in    the transformation 
% (post-processing)in the local coordinate system,     number, 
% distribution, and distance of GPS reference stations, quality    of the 
% GPS/INS postprocessing, correction of the relative errors    through 
% block adjustment of the scan strips, etc. *Total Error is actually 
% magnitude error, calculated by taking the square root of the sum of the 
% square of individual error values. 
% 
% FORMAT:
% 
%    [ ...
%    DX_err, ...
%    DY_err, ...
%    DZ_err ...
%    ] = ...
%    fcn_MapGen_ugvSensorError( ...
%    Scanning_Results, ...
%    Position_Uncertainty, ...
%    Angular_Uncertainty, ...
%    Laser_Uncertainty, ...
%    (fig_num) ...
%    )
% 
% INPUTS:
% 
%     Scanning_Results: a Nx3 matrix, N being the number of received lidar 
%     points or sensor scans
% 
%     Position_Uncertainty: a 1x3 vector of constants defining uncertainty 
%     in the x, y, and z directions
% 
%     Angular_Uncertainty: a 1x3 vector of constants defining uncertainty 
%     in the x, y, and z pointing angle (in degrees)
% 
%     Laser_Uncertainty: a 1x2 vector of constants defining uncertainty in 
%     the LIDAR (see long description)
% 
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a 
%     figure to be drawn showing results.
% 
% 
% OUTPUTS:
% 
%     DX_err: Error in the x position
% 
%     DY_err: Error in the y position
% 
%     DZ_err: Error in the z position
% 
% 
% DEPENDENCIES:
% 
%     fcn_MapGen_checkInputsToFunctions
% 
% 
% EXAMPLES:
% 
% See the script: script_test_fcn_MapGen_ugvSensorError
% for a full test suite.
% 
% This function was written on 2021_07_07 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

% 
% REVISION HISTORY:
% 
% 2021_07_07 by Sean Brennan
% -- first write of function

% 
% TO DO:
% 
% -- fill in to-do items here.

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments 
flag_do_plot = 0;      % Set equal to 1 for plotting 
flag_do_debug = 0;     % Set equal to 1 for debugging 

if flag_do_debug
    fig_for_debug = 961;
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end 

%% check input arguments?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if 1 == flag_check_inputs

    % Are there the right number of inputs?
    if nargin < 4 || nargin > 5
        error('Incorrect number of input arguments')
    end

    % Check the Scanning_Results input, make sure it is 'positive_3column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        Scanning_Results, 'positive_3column_of_numbers');
 
    % Check the Position_Uncertainty input, make sure it is 'positive_3column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        Position_Uncertainty, 'positive_3column_of_numbers',[1 1 ]);
 
    % Check the Angular_Uncertainty input, make sure it is 'positive_3column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        Angular_Uncertainty, 'positive_3column_of_numbers',[1 1 ]);
 
    % Check the Laser_Uncertainty input, make sure it is 'positive_2column_of_numbers' type
    fcn_MapGen_checkInputsToFunctions(...
        Laser_Uncertainty, 'positive_2column_of_numbers',[1 1 ]);
 
end

% Does user want to show the plots?
if  5== nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
        flag_do_plot = 1;
    end
end

%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%See: http://patorjk.com/software/taag/#p=display&f=Big&t=Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§





% Revisions:
%     Nick Carder 2/15/21
%     v3  - for specific use in scripts, intended for handling vectors for use
%     in plotting
%         - equations based on simple 2D trigonometry
%         - beginning to explicitly apply assumptions for 2D flat terrain

%% Initialization and Variable Assignment

%Re-assigning variable names
R =     Scanning_Results{1};
beta =  Scanning_Results{2};
kappa = Scanning_Results{3};

dx0 =   Position_Uncertainty{1};
dy0 =   Position_Uncertainty{2};
dz0 =   Position_Uncertainty{3};

domega =    Angular_Uncertainty{1};
dpsi =      Angular_Uncertainty{2};
dkappa =    Angular_Uncertainty{3};

dbeta = Laser_Uncertainty{1};
dR =    Laser_Uncertainty{2};

clear Scanning_Results Position_Uncertainty Angular_Uncertainty Laser_Uncertainty

%Initialization
% DX=[];   DY=[];   DZ=[];

%% For 3D Cases (Future Work)
%{

%% Influence of Roll (domega) 
% negligible for 2D, flat terrain

dx=0;

% dy=h*(sin(beta+domega)-sin(beta))/cos(beta);
% appoximated to:
dy = D*domega;

DX_omega = -dy*sind(kappa);
DX = [DX , DX_omega]; 

DY_omega = dy*cosd(kappa);
DY = [DY , DY_omega];

% DZ_omega = h*(1-cos(beta+domega)/cos(beta))];
% approximated to:
DZ_omega =  D*domega*tand(beta);
DZ = [DZ , DZ_omega];

%% Influence of Pitch (dpsi)
% negligible for 2D, flat terrain

dx=(-D*sind(dpsi));
% approximates to:
% dx = h*dpsi;

dy=0;

DX_psi = dx*cosd(kappa);
DX = [DX , DX_psi];

DY_psi = dx*sind(kappa);
DY = [DY , DY_psi];

DZ_psi = D*(1-cosd(dpsi));
DZ = [DZ , DZ_psi];

%}
%% Influence of Yaw (dkappa)

DX_kappa = R.*cosd(beta+dbeta+kappa+dkappa)-R.*cosd(beta+dbeta+kappa);
% DX = [DX , DX_kappa];
    
DY_kappa = R.*sind(beta+dbeta+kappa+dkappa)-R.*sind(beta+dbeta+kappa);
% DY = [DY , DY_kappa];

% DZ_kappa = 0;
% DZ = [DZ , DZ_kappa];


%% Influence of Scan Mirror Angle (dbeta)

DX_beta = R.*cosd(beta+dbeta+kappa+dkappa)-R.*cosd(beta+kappa+dkappa);
% DX = [DX , DX_beta];

DY_beta = R.*sind(beta+dbeta+kappa+dkappa)-R.*sind(beta+kappa+dkappa);
% DY = [DY , DY_beta];

% DZ_beta = h*(1 - cosd(beta + dbeta)/cosd(beta));
% approximates to:
% DZ_beta = D*dbeta*tand(beta);
% DZ = [DZ , DZ_beta];

%% Influence of Ranging (dR)

DX_R = dR*cosd(beta+dbeta+kappa+dkappa);
% DX = [DX , DX_R];

DY_R = dR*sind(beta+dbeta+kappa+dkappa);
% DY = [DY , DY_R];

% DZ_R = -dR*cosd(abs(beta));
% DZ = [DZ , DZ_R];

%% Influence of Laser Beam Origin Position (dx0, dy0, dz0)

DX_x0 = dx0*ones(size(R));
DX_y0 = zeros(size(R));
% DX_z0 = zeros(size(dx0));
% DX = [DX , DX_x0, DX_y0];

DY_x0 = zeros(size(R));
DY_y0 = dy0*ones(size(R));
% DY_z0 = zeros(size(dx0));
% DY = [DY , DY_x0, DY_y0];

% DZ_x0 = 0;
% DZ_y0 = 0;
% DZ_z0 = dz0;
% DZ = [DZ , DZ_x0, DZ_y0, DZ_z0];

%% Total Error Values

% 3D case
% DX_err=sqrt(sum(DX.^2,2));
% DY_err=sqrt(sum(DY.^2,2));
% DZ_err=sqrt(sum(DZ.^2));

% no D_omega or D_psi or Dz0 for 2D
DX_err= sqrt(DX_kappa.^2 + DX_beta.^2 + DX_R.^2 + DX_x0.^2 + DX_y0.^2);
DY_err= sqrt(DY_kappa.^2 + DY_beta.^2 + DY_R.^2 + DY_x0.^2 + DY_y0.^2);

% for displaying individual error values
% DX_table=table(h, beta, kappa, DX_omega, DX_psi, DX_kappa, DX_beta, DX_R, DX_x0, DX_y0, DX_z0);
% DY_table=table(h, beta, kappa, DY_omega, DY_psi, DY_kappa, DY_beta, DY_R, DY_x0, DY_y0, DY_z0);
% DZ_table=table(h, beta, kappa, DZ_omega, DZ_psi, DZ_kappa, DZ_beta, DZ_R, DZ_x0, DZ_y0, DZ_z0);

% 

%§
%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if flag_do_plot
    % Nothing to plot here
end % Ends the flag_do_plot if statement    

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends the function

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


end
