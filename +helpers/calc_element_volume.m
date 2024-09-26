function vol = calc_element_volume(elements, nodes)
len = size(elements,1);
vol = zeros(len,1);
for ii = 1:len
    element = elements(ii,:);
    % Extract the node coordinates
    p1 = nodes(element(1), :);
    p2 = nodes(element(2), :);
    p3 = nodes(element(3), :);
    p4 = nodes(element(4), :);
    
    % Volume of a tetrahedron = (1/6) * |det([p2-p1; p3-p1; p4-p1])|
    vol(ii,1) = abs(det([p2 - p1; p3 - p1; p4 - p1])) / 6;
end
end
