% script_test_fcn_MapGen_cropPolytopeToRange

%% Test case 1: simple square

verticies = [1 1; 1 -1; -1 -1; -1 1]*0.5;
interior_point = [0 0];
[cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies,interior_point);

figure(1);
clf;
hold on;
plot(verticies(:,1),verticies(:,2),'r-');
plot(cropped_vertices(:,1),cropped_vertices(:,2),'g-');


%% Test case 2: simple triangle
verticies = [1 -0.5; -1 -1; -0.5 1];
interior_point = [0 0];
[cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies,interior_point);

figure(1);
clf;
hold on;
plot(verticies(:,1),verticies(:,2),'r-');
plot(cropped_vertices(:,1),cropped_vertices(:,2),'g-');

%% Test case 3: vertices enclosing region
verticies = [1 1; 1 -1; -1 -1; -1 1]*2;
interior_point = [0 0];
[cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies,interior_point);

figure(1);
clf;
hold on;
plot(verticies(:,1),verticies(:,2),'r-');
plot(cropped_vertices(:,1),cropped_vertices(:,2),'g-');

%% Test case 4: vertices all within region
verticies = [ 0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
interior_point = [0.5 0.5];
[cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies,interior_point);

figure(1);
clf;
hold on;
plot(verticies(:,1),verticies(:,2),'r-');
plot(cropped_vertices(:,1),cropped_vertices(:,2),'g-');

%% Test case 5: infinity in the numbers
verticies = [...
   -0.0457    0.0471
   -0.4075    0.0851
   -2.7148    0.1857
  -11.3670    0.2644
       Inf       Inf
    0.3841  -53.6209
    0.2255   -6.8725
    0.1203   -1.0609
    0.0352   -0.0163
   -0.0081    0.0368];

interior_point = [0 0];
[cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies,interior_point);

figure(1);
clf;
hold on;
plot(verticies(:,1),verticies(:,2),'r-');
plot(cropped_vertices(:,1),cropped_vertices(:,2),'g-');

%% Points on edges
verticies = [
   0.951225451552411   0.038481963258910
   0.951225451552411                   0
   1.000000000000000   0.055571820180916
   0.972711461970130   0.055571820180916
   0.967138332311524   0.053930059059772];
interior_point = [0.978515625000000   0.035665294924554];
[cropped_vertices] = fcn_MapGen_cropPolytopeToRange(verticies,interior_point);

figure(1);
clf;
hold on;
plot(verticies(:,1),verticies(:,2),'r-');
plot(cropped_vertices(:,1),cropped_vertices(:,2),'g-');
