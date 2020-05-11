function [area] = area(A,B,C)
% Input: A,B,C are the three sommets of a triangle.
% This function calculates the area of the triangle whose sommets are A,B
% and C.

T1 = A(1)*(B(2)-C(2));
T2 = B(1)*(C(2)-A(2));
T3 = C(1)*(A(2)-B(2));

area = abs(1/2*(T1+T2+T3));