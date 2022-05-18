% script_fcn_MapGen_polytopeShrinkFromEdges
% Tests function: fcn_MapGen_polytopeShrinkFromEdges

% REVISION HISTORY:
% 2022_01_17
% -- first written by S. Harnett using
% script_test_fcn_MapGen_polytopesShrinkToRadius as a template

%% Set up variables
close all;clear all;
polytopes = fcn_MapGen_haltonVoronoiTiling([1 100]);
%% big shrink
fig_num = 1;
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);
pre_shrink_stats = fcn_MapGen_polytopesStatistics(trim_polytopes);
R_bar_initial = pre_shrink_stats.average_max_radius;
%% Basic example of uniform shrinking
G_bar_all = [];
G_perim_all = [];
des_gap_size_all = [];
r_occ_meas_all = [];
rho_all = [];
A_unocc_est_density_all = [];
A_unocc_est_perim_all = [];
A_unocc_est_perim_improved_all = [];
r_unocc_meas_all = [];
A_unocc_est_parallelogram_all = [];
A_unocc_est_avg_parallelogram_all = [];
A_unocc_est_parallelograms_and_kites_avg_all = [];
A_unocc_est_parallelograms_and_kites_all = [];
for i = linspace(0.001,0.08,20)
    fig_num = 2;
    des_gap_size = i;
    shrunk_polytopes1=...
        fcn_MapGen_polytopesShrinkFromEdges(...
        trim_polytopes,des_gap_size);
    field_stats = fcn_MapGen_polytopesStatistics(shrunk_polytopes1);
    r_occ_meas = field_stats.occupancy_ratio; % calculated occupancy ratio
    r_unocc_meas = field_stats.unoccupancy_ratio;
    G_bar = field_stats.average_gap_size_G_bar;
    G_perim = field_stats.perimeter_gap_size;
    rho = field_stats.linear_density_mean;
    P_tot = field_stats.total_perimeter;
    N_vert = field_stats.NtotalVertices;
    A_unocc_est_density = des_gap_size^2*rho; % theoretial occupancy ratio from gap size
    A_unocc_est_perim = 1/2*des_gap_size*P_tot;
    A_unocc_est_perim_improved = A_unocc_est_perim + N_vert/3*des_gap_size^2*sqrt(3)/4;
    G_bar_all = [G_bar_all; G_bar];
    G_perim_all = [G_perim_all; G_perim];
    des_gap_size_all = [des_gap_size_all; des_gap_size];
    r_occ_meas_all = [r_occ_meas_all; r_occ_meas];
    rho_all = [rho_all; rho];
    r_unocc_meas_all = [r_unocc_meas_all; r_unocc_meas];
    A_unocc_est_density_all = [A_unocc_est_density_all; A_unocc_est_density];
    A_unocc_est_perim_all = [A_unocc_est_perim_all; A_unocc_est_perim];
    A_unocc_est_perim_improved_all = [A_unocc_est_perim_improved_all; A_unocc_est_perim_improved];
    % TODO make estimate subtracting one parallelogram from each vertex
    % TODO make the same code subtracting one parellelogram per angle for the average angle size
    % note this is interior angles
    angles = field_stats.angle_column_no_nan;
    parallelogram_areas = des_gap_size/2*des_gap_size/2*sin(angles);
    total_parallelogram_area = sum(parallelogram_areas);
    average_angle = field_stats.average_vertex_angle;
    total_avg_parallelogram_area = des_gap_size/2*des_gap_size/2*sin(average_angle)*length(parallelogram_areas);
    A_unocc_est_avg_parallelogram = A_unocc_est_perim + total_avg_parallelogram_area;
    A_unocc_est_parallelogram = A_unocc_est_perim + total_parallelogram_area;
    A_unocc_est_avg_parallelogram_all = [A_unocc_est_avg_parallelogram_all; A_unocc_est_avg_parallelogram];
    A_unocc_est_parallelogram_all = [A_unocc_est_parallelogram_all; A_unocc_est_parallelogram];
    % TODO if the angle is over 90, 180-the angle forms a roughly isoceles triangle,
    % but it actually forms a kite
    % assert(isequal(round(G_bar,4),round(des_gap_size,4)));
    % assert(isequal(round(r_unocc_meas,4),round(estimated_unoccupancy_ratio,4)));
    angles_acute_logical = angles <= pi/2;
    angles_acute = nonzeros(angles_acute_logical.*angles)';
    avg_acute = mean(angles_acute);
    angles_obtuse_logical = angles > pi/2;
    angles_obtuse = nonzeros(angles_obtuse_logical.*angles)';
    avg_obtuse = mean(angles_obtuse);
    parallelogram_areas_acute = des_gap_size/2*des_gap_size/2*sin(angles_acute);
    total_parallelogram_area_acute_avg = length(angles_acute)*des_gap_size/2*des_gap_size/2*sin(avg_acute);
    total_parallelogram_area_acute = sum(parallelogram_areas_acute);
    % kite area calculation
    a = des_gap_size/2*cos((pi-angles_obtuse)/2);
    b = des_gap_size/2*sin((pi-angles_obtuse)/2);
    c = b/tan(angles_obtuse/2);
    a_avg = des_gap_size/2*cos((pi-avg_obtuse)/2);
    b_avg = des_gap_size/2*sin((pi-avg_obtuse)/2);
    c_avg = b_avg/tan(avg_obtuse/2);
    kite_areas_obtuse = (a+c).*(2*b)/2;
    kite_areas_obtuse_avg = (a_avg+c_avg).*(2*b_avg)/2;
    total_kite_areas = sum(kite_areas_obtuse);
    total_kite_areas_avg = length(angles_obtuse)*kite_areas_obtuse_avg;
    A_unocc_est_parallelograms_and_kites = A_unocc_est_perim + total_parallelogram_area_acute + total_kite_areas;
    A_unocc_est_parallelograms_and_kites_all = [A_unocc_est_parallelograms_and_kites_all; A_unocc_est_parallelograms_and_kites];
    A_unocc_est_parallelograms_and_kites_avg = A_unocc_est_perim + total_parallelogram_area_acute_avg + total_kite_areas_avg;
    A_unocc_est_parallelograms_and_kites_avg_all = [A_unocc_est_parallelograms_and_kites_avg_all; A_unocc_est_parallelograms_and_kites_avg];
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
ylim([-100,140])
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_density_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_perim_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_perim_improved_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_parallelogram_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_avg_parallelogram_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_parallelograms_and_kites_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_parallelograms_and_kites_avg_all - r_unocc_meas_all)./r_unocc_meas_all*100)
legend('density estimate',...
'perimeter estimate',...
'perimeter estimate with triangles',...
'perimeter estimate with parallelograms',...
'perimeter estimate with average parallelograms',...
'perimeter estimate with parllelograms and kites',...
'perimeter estimate with average parallelograms and kites')
xlabel('desired or commanded gap size as percent of initial average max radius')
ylabel('unoccpancy ratio estimate percent error')
% print('C:\Users\sjh6473\github\sjharnett\figures\exported\after_gvsets\runocc_err','-dpng','-r300');
% savefig('C:\Users\sjh6473\github\sjharnett\figures\figs\after_gvsets\runocc_err')

figure
hold on
box on
grid on
plot(des_gap_size_all./R_bar_initial*100,r_unocc_meas_all,'k-');
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_density_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_perim_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_perim_improved_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_parallelogram_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_avg_parallelogram_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_parallelograms_and_kites_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_parallelograms_and_kites_avg_all);
legend('measured unoccupancy ratio',...
'density estimate',...
'perimeter estimate',...
'perimeter estimate with triangles',...
'perimeter estimate with parallelograms',...
'perimeter estimate with average parallelograms',...
'perimeter estimate with parllelograms and kites',...
'perimeter estimate with average parallelograms and kites')
xlabel('desired or commanded gap size as percent of initial average max radius')
ylabel('unoccpancy ratio')
