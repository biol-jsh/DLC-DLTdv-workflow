function area = triangle_area(p1, p2, p3)
%TRIANGLE_AREA Compute the area of a triangle in 3D space.
%
%   INPUTS:
%   p1 - [x1, y1, z1], coordinates of the first vertex of the triangle
%   p2 - [x2, y2, z2], coordinates of the second vertex of the triangle
%   p3 - [x3, y3, z3], coordinates of the third vertex of the triangle
%
%   OUTPUT:
%   area - the area of the triangle

% Remove singleton dimensions from the input arrays
p1 = squeeze(p1);
p2 = squeeze(p2);
p3 = squeeze(p3);

% Compute vectors that define the triangle
v1 = p2 - p1;
v2 = p3 - p1;

% Compute the cross product of these vectors
cross_product = cross(v1, v2);

% Compute the area of the triangle
area = 0.5 * norm(cross_product);
end
