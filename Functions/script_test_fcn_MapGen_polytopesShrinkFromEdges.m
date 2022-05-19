% script_fcn_MapGen_polytopeShrinkFromEdges
% Tests function: fcn_MapGen_polytopeShrinkFromEdges

% REVISION HISTORY:
% 2022_01_17
% -- first written by S. Harnett using
% script_test_fcn_MapGen_polytopesShrinkToRadius as a template

%% Set up variables
close all;clear all;
polytopes = fcn_MapGen_haltonVoronoiTiling([1 100]);
fig_num = 1;
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);
pre_shrink_stats = fcn_MapGen_polytopesStatistics(trim_polytopes);
R_bar_initial = pre_shrink_stats.average_max_radius;
G_bar_all = [];
G_perim_all = [];
des_gap_size_all = [];
rho_all = [];
%% loop through gap sizes from rather small to quite large
for i = linspace(0.001,0.08,20)
    fig_num = 2;
    des_gap_size = i;
    shrunk_polytopes1=...
        fcn_MapGen_polytopesShrinkFromEdges(...
        trim_polytopes,des_gap_size);
    field_stats = fcn_MapGen_polytopesStatistics(shrunk_polytopes1);
    G_bar = field_stats.average_gap_size_G_bar;
    G_perim = field_stats.perimeter_gap_size;
    rho = field_stats.linear_density_mean;
    G_bar_all = [G_bar_all; G_bar];
    G_perim_all = [G_perim_all; G_perim];
    des_gap_size_all = [des_gap_size_all; des_gap_size];
    rho_all = [rho_all; rho];
end

% Defaults for this blog post
width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

% The new defaults will not take effect if there are any open figures. To
% use them, we close all figures, and then repeat the first example.
close all;

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultAxesFontSize',fsz);
set(0,'defaultLegendFontSize',fsz);
set(0,'defaultAxesLineWidth',alw);
% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
set(0,'defaultAxesBox','on')

figure
hold on
box on
grid on
plot(des_gap_size_all,(G_perim_all-des_gap_size_all))
plot(des_gap_size_all,(G_bar_all-des_gap_size_all))
legend('perimeter estimate',...
'density estimate')
xlabel('desired or commanded gap size')
ylabel('gap size estimate error')
% print('C:\Users\sjh6473\github\sjharnett\figures\exported\after_gvsets\g_bar_err','-dpng','-r300');
% savefig('C:\Users\sjh6473\github\sjharnett\figures\figs\after_gvsets\g_bar_err')

figure
hold on
box on
grid on
plot(des_gap_size_all./R_bar_initial*100,G_bar_all)
plot(des_gap_size_all./R_bar_initial*100,G_perim_all)
legend('density estimate',...
'perimeter estimate')
xlabel('desired or commanded gap size as percent of initial average max radius')
ylabel('gap size estimate')
% print('C:\Users\sjh6473\github\sjharnett\figures\exported\after_gvsets\runocc_err','-dpng','-r300');
% savefig('C:\Users\sjh6473\github\sjharnett\figures\figs\after_gvsets\runocc_err')
