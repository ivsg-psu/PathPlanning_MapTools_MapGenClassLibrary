% script_test_fcn_MapGen_nameToMap
% Tests function: fcn_MapGen_nameToMap

% REVISION HISTORY:
% 2021_06_06
% -- first written by S. Brennan.

close all;

%% Basic Example:
map_name = "HST 1 100 SQT 0 1 0 1 SMV 0.01 0.001 1e-6 1111";
plot_flag = 1; disp_name = 0; fig_num = []; line_style = 'r-';
line_width = 2;
[polytopes,fig]=...
    fcn_MapGen_nameToMap(map_name,plot_flag,disp_name,fig_num,line_style,line_width); %#ok<*ASGLU>

assert(true);

%% Advanced Example
map_name = "HST 30 450 SQT 0 1 0 1 SMV 0.02 0.005 1e-6 1234";
plot_flag = 1; disp_name = [1, 0.05 -0.05, 0.5 0.5 0.5, 10];
fig_num = 999; line_style = '-'; line_width = 2; color = [0 0 1];
axis_limits = [0 1 -0.1 1]; axis_style = 'square';
fill_info = [1 1 0 1 0.5];
[polytopes,fig]=fcn_MapGen_nameToMap(map_name,plot_flag,disp_name,fig_num,line_style,line_width,color,axis_limits,axis_style,fill_info);
