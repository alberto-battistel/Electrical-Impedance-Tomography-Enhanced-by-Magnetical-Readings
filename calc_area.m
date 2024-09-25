function [areas] = calc_area(nodes, elems)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

areas = zeros(length(elems),1);

for ii = 1:length(elems)
    % Get the vertices of the i-th triangle
    x1 = nodes(elems(ii,1), 1);
    y1 = nodes(elems(ii,1), 2);
    x2 = nodes(elems(ii,2), 1);
    y2 = nodes(elems(ii,2), 2);
    x3 = nodes(elems(ii,3), 1);
    y3 = nodes(elems(ii,3), 2);


    areas(ii) = 0.5 * abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
end
end