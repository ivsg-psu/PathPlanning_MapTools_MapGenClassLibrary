% Scratch Paper Examples for Polytope Error Bubble Generation
% Nick Carder
% 3/26/21

%% Generating Error Bubbles and Plotting

% create polytopes
stretch = [200, 200]; % stretch in the x and y directions
Halton_range = [1 1000]; % range of Halton points to use to generate the tiling
polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,stretch);

% Plot the polytopes
fig_num = 22;
line_width = 2;
axis_limits = [0 200 0 200];
axis_stype = 'square';
fcn_MapGen_plotPolytopes(polytopes,fig_num,'b',line_width,axis_limits,axis_stype);


% remove the edge polytopes that extend past the high and low points
fig_num = 23;
bounding_box = [0 0; 200 200];
trimmed_polytopes = ...
    fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);

%shrink polytopes to create space
scale = ones(1,size(trim_polytopes,2));
des_rad = 1; sigma_radius = 0.5; min_rad = 0.25;
shrunk_polytopes2=fcn_polytope_editing_shrink_to_average_max_radius_with_variance(trim_polytopes,des_rad,sigma_radius,min_rad);








%generate error bubbles
[err] = err_ugv_bubble_v3(shrunk_polytopes2, 0, 5);

%verify
figure('name','UGV Positioning Bubbles')
ax.ugv1=gca;
hold on
for ii=1:length(shrunk_polytopes2)
    plot(shrunk_polytopes2(ii).vertices(:,1),shrunk_polytopes2(ii).vertices(:,2),'k')
    plot(ax.ugv1,err(ii).circ_x(err(ii).bubble),err(ii).circ_y(err(ii).bubble),'r')
end
hold off
xlabel('X Distance [m]')
ylabel('Y Distance [m]')
title('UGV Positioning Bubbles')
legend('Real Object','Perceived Object')
grid on
axis equal

%% Testing work with contour plots

%generating R and beta values for a grid of points
x=linspace(0,200,200);
y=linspace(100,-100,200);
[X,Y] = meshgrid(x,y);

R=sqrt(X.^2+Y.^2);  % perceived distance, nearly constant for flat object
beta=rad2deg(atan(Y./X));
kappa=zeros(size(R));

[err.x, err.y] = err_ugv_v3 ({R, beta, kappa}, ...
    {0.08, 0.08, 0.08}, {0.03, 0.03, 0.4}, {0.02, -0.05});

err.R = sqrt(err.x.^2 + err.y.^2);

figure('name','Error in X')
contour(x,y,err.x,'ShowText','on')
h=colorbar;
title('Error in X Dimension [m]')
xlabel('X [m]')
ylabel('Y [m]')

figure('name','Error in Y')
contour(x,y,err.y,'ShowText','on')
h=colorbar;
title('Error in Y Dimension [m]')
xlabel('X [m]')
ylabel('Y [m]')

figure('name','Total Error (Bubble Radius) [m]')
contour(x,y,err.R,'ShowText','on')
h=colorbar;
title('Total Error (Bubble Radius) [m]')
xlabel('X [m]')
ylabel('Y [m]')