// Gmsh project created on Thu Aug 23 16:25:03 2018
h = 5e-3;
//+
Point(1) = {0, 0, 0, h};
//+
Point(2) = {0.5, 0, 0, h};
//+
Point(3) = {0.5, 0.5, 0, h};
//+
Point(4) = {1, 0.5, 0, h};
//+
Point(5) = {1, 1, 0, h};
//+
Point(6) = {0, 1, 0, h};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Line Loop(1) = {6, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Surface("Domain") = {1};
