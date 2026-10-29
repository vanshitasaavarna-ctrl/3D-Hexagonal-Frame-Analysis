% =========================================================================
% CVL7411 ADVANCED STRUCTURAL ANALYSIS - PROGRAMMING ASSIGNMENT GROUP (D)
% 3D SINGLE-STOREY HEXAGONAL SPACE FRAME WITH RIGID SLAB
% =========================================================================

function [coords, members] = generate_geometry(L, H)
coords = zeros(12, 3);

% Base Nodes 1-6 (Z = 0)
for i = 1:6
    theta = (i - 1) * (pi / 3);
    coords(i, :) = [L * cos(theta), L * sin(theta), 0];
end

% Roof Nodes 7-12 (Z = H)
for i = 1:6
    theta = (i - 1) * (pi / 3);
    coords(i + 6, :) = [L * cos(theta), L * sin(theta), H];
end

members = zeros(12, 2);

% Members 1-6: Vertical Columns
for i = 1:6
    members(i, :) = [i, i + 6];
end

% Members 7-12: Perimeter Roof Ring Beams
for i = 1:6
    node_start = i + 6;
    node_end = mod(i, 6) + 7;
    members(i + 6, :) = [node_start, node_end];
end
end