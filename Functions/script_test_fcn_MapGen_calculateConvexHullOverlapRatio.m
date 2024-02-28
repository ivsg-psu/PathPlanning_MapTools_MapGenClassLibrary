% script_test_fcn_MapGen_calculateConvexHullOverlapRatio
% Tests: fcn_MapGen_calculateConvexHullOverlapRatio

%
% REVISION HISTORY:
%
% 2024_02_28 by Steve Harnett
% -- first write of script
%%%%%%%%%%%%%%§

close all; clear all; clc;

%% basic polytope case
polytopes(1).vertices = [0 0; 10 0; 10 1; 0 1; 0 0];
polytopes(2).vertices = polytopes(1).vertices+[8,0];
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes);
fig_num = 1;
[ ...
convex_hull_overlap_ratio,...
A_overlap,...
A_occupied...
] = ...
fcn_MapGen_calculateConvexHullOverlapRatio( ...
polytopes, ...
fig_num...
)

%% flood plain map case
load(strcat(pwd,'\..\..\PathPlanning_GridFreePathPlanners_BoundedAStar\Test_Fixtures\flood_plains\flood_plain_4.mat'));
shrunk_polytopes = flood_plain_4;

%% convert from LLA to QGS84
centre_co_avg_alt = 351.7392;
new_polytopes = [];
for i = 1:length(shrunk_polytopes)
    poly = shrunk_polytopes(i);
    lats = poly.vertices(:,2);
    longs = poly.vertices(:,1);
    alts = centre_co_avg_alt*ones(size(lats));
    wgs_verts = [];
    for j = 1:length(lats)
        xyz = INTERNAL_WGSLLA2xyz(lats(j),longs(j),alts(j));
        xyz = xyz/1000;
        wgs_verts(j,:) = [xyz(1),xyz(2)];
    end
    new_polytopes(i).vertices = wgs_verts;
end
shrunk_polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(new_polytopes);

fig_num = fig_num + 1;
[ ...
convex_hull_overlap_ratio,...
A_overlap,...
A_occupied...
] = ...
fcn_MapGen_calculateConvexHullOverlapRatio( ...
shrunk_polytopes, ...
fig_num...
)

xlabel('x [km]')
ylabel('y [km]')
xlim([min([shrunk_polytopes.xv]) max([shrunk_polytopes.xv])]);
ylim([min([shrunk_polytopes.yv]) max([shrunk_polytopes.yv])]);

function xyz = INTERNAL_WGSLLA2xyz(wlat, wlon, walt)
    %Function xyz = wgslla2xyz(lat, lon, alt) returns the
    %equivalent WGS84 XYZ coordinates (in meters) for a
    %given geodetic latitude "lat" (degrees), longitude "lon"
    %(degrees), and altitude above the WGS84 ellipsoid
    %in meters.  Note: N latitude is positive, S latitude
    %is negative, E longitude is positive, W longitude is
    %negative.
    %
    %Ref: Decker, B. L., World Geodetic System 1984,
    %Defense Mapping Agency Aerospace Center.

    A_EARTH = 6378137;
    flattening = 1/298.257223563;
    NAV_E2 = (2-flattening)*flattening; % also e^2
    deg2rad = pi/180;

    slat = sin(wlat*deg2rad);
    clat = cos(wlat*deg2rad);
    r_n = A_EARTH/sqrt(1 - NAV_E2*slat*slat);
    xyz = [ (r_n + walt)*clat*cos(wlon*deg2rad);
            (r_n + walt)*clat*sin(wlon*deg2rad);
            (r_n*(1 - NAV_E2) + walt)*slat ];

    if ((wlat < -90.0) | (wlat > +90.0) | (wlon < -180.0) | (wlon > +360.0))
        error('WGS lat or WGS lon out of range');
    end
end
