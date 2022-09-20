p1 = polyshape([1 0 0 1],[1 1 0 0]);
p2 = polyshape([1.5 0.5 0.5 1.5],[1.5 1.5 0.5 0.5]);
figure
hold on
plot(p1)
plot(p2)
p3 = intersect(p1,p2)
plot(p3)
p1_new = subtract(p1,p3)
p2_new = subtract(p2,p3)
figure
hold on
plot(p1_new)
plot(p2_new)
plot(p3)
p1_new = simplify(p1_new)
p2_new = simplify(p2_new)
p3 = simplify(p3)