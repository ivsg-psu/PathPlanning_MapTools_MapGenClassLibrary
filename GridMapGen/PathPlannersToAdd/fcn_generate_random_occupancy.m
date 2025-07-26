function occupancy_map = fcn_generate_random_occupancy(n,m,varargin)
% FCN_GENERATE_RANDOM_OCCUPANCY  generates a random occupancy map
%
% OCCUPANCY_MAP = FCN_GENERATE_RANDOM_OCCUPANCY(N,M) returns an N-by-M
% matrix containing 1's and 0's, with 1 being occupied and 0 clear. Default
% values are used where 20% of the map is occupied, dilation is level 3.
%
% OCCUPANCY_MAP = FCN_GENERATE_RANDOM_OCCUPANCY(N,M,R) allows the user to
% specify the ratio of occupied elements to total elements. For example, if
% R = 0.5, then 50% of the elements will be occupied. The default value is
% 20%.
%
% OCCUPANCY_MAP = FCN_GENERATE_RANDOM_OCCUPANCY(N,M,R,D) allows the user to
% specify the dilation level (integer), which in turn affects the
% smoothness. Dilation levels of 0 are the same as the random plot. A value
% of 1, 2, or 3 usually gives good maps, though larger values can be used.
% The default is 2.
% 
% OCCUPANCY_MAP = FCN_GENERATE_RANDOM_OCCUPANCY(N,M,R,D,SEED_MAP) allows
% the user to specify a seed map as an N-by-M matrix. This is especially
% useful to compare the same map data at different resolutions or levels of
% smoothing. Note: the size of the seed map overwrites the values of N and
% M in the first two arguments.
% 
% Examples:
%      
%      % BASIC example
%      occupancy_map = fcn_generate_random_occupancy(100,100);
%      image(occupancy_map+1);
%      colormap([1 1 1;0 0 0])
%
%      % ADVANCED example
%      % Give a specified seed matrix, an occupancy ratio of 40%, and iterate over various smoothness
%      seed_map = rand(100,100);
%      r = 0.4;
%      total_maps = seed_map<r+1;
%      levels = 20;
%      for D = 1:levels
%          occupancy_map = fcn_generate_random_occupancy(100,100,r,D,seed_map);
%          total_maps = total_maps + occupancy_map;
%          maps{D} = occupancy_map; %#ok<AGROW>
%      end
%      figure;
%      for D = 1:levels
%          shade = (levels-D)*0.75/levels;
%          image(maps{D}+1); colormap([1 1 1; shade shade shade]); hold on;
%          percent_occupied = fcn_generate_occupancy_stats(maps{D});
%          title(sprintf('Percent occupied at dialation %d, r of %.2f: %.4f\n',D,r,percent_occupied));
%          pause(1);
%      end
% 
% This function was written on 2008_10_18 by Dr. Sean Brennan
% Questions or comments? sbrennan@psu.edu 
%


%% Check input arguments and, if needed, generate initial occupancy map
if nargin<2 || nargin>5
    error('Incorrect number of arguments');
end

% check if user has specified ratio of occupied
if nargin>=3
    r = varargin{1};
else
    r = 0.2; % Only 20% should be occupied
end

% Check if user has specified dilation factor
if nargin>=4
    dilation = varargin{2};
else
    dilation = 2;  % Dilation = 2 gives fairly smooth results
end

% Check if user has specified the initial random map
if nargin==5
    rand_map = varargin{3};
    [n,m] = size(rand_map); %#ok<NASGU>
else
    % Generate inital occupancy map
    rand_map = rand(n,m);
end

%% Initialize occupancy map (for debugging)
%occupancy_map = (rand_map<r);   %#ok<NASGU>

%% Display map (for debugging)
% figure('Name','Original Random','NumberTitle','off');
% set(gcf,'Position',[16 149 293 231]);
% image(occupancy_map+1);
% colormap([1 1 1;0 0 0])
% percent_occupied = fcn_generate_occupancy_stats(occupancy_map);
% title(sprintf('Percent occupied at dialation %d, r of %.2f: %.4f\n',1,r,percent_occupied));

%% Dilate on original random image with bad dialation function
rand_map_d1 = fcn_dilate_by_n_bad(rand_map,dilation);

%% Iterate to find the right value of r that preserves ratio
low = min(min(rand_map_d1));
high = max(max(rand_map_d1));

levels = low:(high-low)/100:high;
percentages = zeros(length(levels),1);
for i=1:length(levels)
    temp_occupancy = rand_map_d1<levels(i);
    percentages(i) = fcn_generate_occupancy_stats(temp_occupancy);
end
indices = find(percentages>r);
r_d = levels(indices(1));


%% Save final result
occupancy_map = rand_map_d1<r_d;

%% 





function dilated = fcn_dilate_by_n_bad(matrix,level)
% FCN_DIALATE_BY_N  dilates an image by n pixels
% DIALATED = FCN_DILATE_BY_N(matrix,n), given an N-by-M matrix,
% returns an N-by-M matrix that has been dialated by n pixels. NOTE: this
% is NOT true dilation. For true dilation, see the image processing
% toolbox.
% 
% Examples:
%       
%      test_occupancy = rand(100,100)<0.2; % 20 percent occupied
%      dilated = fcn_dilate_by_n(test_occupancy,1)>=1;
%      % Display dialated on top of original
%      image(2*test_occupancy+dilated+1);
%      colormap([1 1 1;0.5 0.5 0.5;0 0 0])
%
% This function was written on 2008_10_18 by Dr. Sean Brennan
% Questions or comments? sbrennan@psu.edu 
%


%% Check input arguments
if nargin~=2
    error('Two arguments for function call. First is number of rows, second is number of columns.');
end

[n,m] = size(matrix);
dilated = matrix;

%% Generate inital occupancy map
multiplier1 = diag(ones(n,1),0) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
multiplier2 = diag(ones(m,1),0) + diag(ones(m-1,1),1) + diag(ones(m-1,1),-1);
for i=1:level
    dilated = multiplier1*dilated*multiplier2;
    dilated(:,1) = max(max(dilated));
    dilated(:,end) = max(max(dilated));
    dilated(1,:) = max(max(dilated));
    dilated(end,:) = max(max(dilated));
end

%% 




function percent_occupied = fcn_generate_occupancy_stats(matrix)
% FCN_GENERATE_OCCUPANCY_STATS  calculates percentage occupancy
% PERCENT_OCCUPIED = FCN_GENERATE_OCCUPANCY_STATS(matrix), given an N-by-M matrix,
% returns percentage of pixels greater than or equal to 1 (occupied).
% 
% Examples:
%       
%      occupancy_map = fcn_generate_random_occupancy(100,100);
%      percent_occupied = fcn_generate_occupancy_stats(occupancy_map)
%
% This function was written on 2008_10_18 by Dr. Sean Brennan
% Questions or comments? sbrennan@psu.edu 
%


%% Check input arguments
if nargin~=1
    error('One argument needed for function call, an occupancy matrix');
end

[n,m] = size(matrix);

%% Calculate occupancy statistics
%percent_occupied = sum(sum(matrix>=1))/(n*m);
percent_occupied = sum(sum(matrix))/(n*m);


