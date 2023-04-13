close all; clear all; clc
n = 25;
m = 25;
sz = [n m]; % size of board
surface_grid = rand(n,m);
cells = [];
figure; hold on
for ind = 1:1:n
    yvec = [(ind-1)/n (ind-1)/n ind/n ind/n];
    for jind = 1:1:m
        xvec = [(jind-1)/m jind/m jind/m (jind-1)/m];
        face_alpha = surface_grid(ind,jind);
        fill(xvec,yvec,'black','FaceAlpha',face_alpha)
        cells(ind*jind).xv = xvec;
        cells(ind*jind).yv = yvec;
        cells(ind*jind).cost = face_alpha;
    end
end

% make a grid Thrusday
% convert some grid cells to polytopes Friday-Monday
% add these two a set of other polytopes Monday
% subtract polytopes from background and triangulate Monday-Tuesday