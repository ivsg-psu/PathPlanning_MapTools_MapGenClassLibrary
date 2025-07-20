function occupancy_map = fcn_GridMapGen_generateRandomOccupancyMap_OLD(n,m,varargin)
% fcn_GridMapGen_generateRandomOccupancyMap_OLD  generates a random occupancy map
%
% OCCUPANCY_MAP = fcn_GridMapGen_generateRandomOccupancyMap_OLD(N,M) returns an N-by-M
% matrix containing 1's and 0's, with 1 being occupied and 0 clear. Default
% values are used where 20% of the map is occupied, dilation is level 3.
%
% OCCUPANCY_MAP = fcn_GridMapGen_generateRandomOccupancyMap_OLD(N,M,R) allows the user to
% specify the ratio of occupied elements to total elements. For example, if
% R = 0.5, then 50% of the elements will be occupied. The default value is
% 20%.
%
% OCCUPANCY_MAP = fcn_GridMapGen_generateRandomOccupancyMap_OLD(N,M,R,D) allows the user to
% specify the dilation level (integer), which in turn affects the
% smoothness. Dilation levels of 0 are the same as the random plot. A value
% of 1, 2, or 3 usually gives good maps, though larger values can be used.
% The default is 2.
% 
% OCCUPANCY_MAP = fcn_GridMapGen_generateRandomOccupancyMap_OLD(N,M,R,D,SEED_MAP) allows
% the user to specify a seed map as an N-by-M matrix. This is especially
% useful to compare the same map data at different resolutions or levels of
% smoothing. Note: the size of the seed map overwrites the values of N and
% M in the first two arguments.
% 
% Examples:
%      
%      % BASIC example
%      occupancy_map = fcn_GridMapGen_generateRandomOccupancyMap_OLD(100,100);
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
%          occupancy_map = fcn_GridMapGen_generateRandomOccupancyMap_OLD(100,100,r,D,seed_map);
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
% percent_occupied = fcn_GridMapGen_dilateOccupancyStats(occupancy_map);
% title(sprintf('Percent occupied at dialation %d, r of %.2f: %.4f\n',1,r,percent_occupied));

%% Dilate on original random image with bad dialation function
rand_map_d1 = fcn_GridMapGen_dilateByN(rand_map,dilation,18181);

%% Iterate to find the right value of r that preserves ratio
low = min(min(rand_map_d1));
high = max(max(rand_map_d1));

levels = low:(high-low)/100:high;
percentages = zeros(length(levels),1);
for i=1:length(levels)
    temp_occupancy = rand_map_d1<levels(i);
    percentages(i) = fcn_GridMapGen_dilateOccupancyStats(temp_occupancy);
end
indices = find(percentages>r);
r_d = levels(indices(1));


%% Save final result
occupancy_map = rand_map_d1<r_d;

%% 
end



