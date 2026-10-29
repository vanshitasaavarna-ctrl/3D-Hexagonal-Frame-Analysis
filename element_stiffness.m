% =========================================================================
% CVL7411 ADVANCED STRUCTURAL ANALYSIS - PROGRAMMING ASSIGNMENT GROUP (D)
% 3D SINGLE-STOREY HEXAGONAL SPACE FRAME WITH RIGID SLAB
% =========================================================================

function [K_global, T, k_loc] = element_stiffness(node1_coords, node2_coords, E, G, A, Iy, Iz, J)
dx = node2_coords(1) - node1_coords(1);
dy = node2_coords(2) - node1_coords(2);
dz = node2_coords(3) - node1_coords(3);
L_e = sqrt(dx^2 + dy^2 + dz^2);

% Local Element Stiffness Matrix (12 x 12)
k_loc = zeros(12, 12);
EA_L = E * A / L_e;
k_loc([1, 7], [1, 7]) = [EA_L, -EA_L; -EA_L, EA_L];

GJ_L = G * J / L_e;
k_loc([4, 10], [4, 10]) = [GJ_L, -GJ_L; -GJ_L, GJ_L];

EIz_12 = 12 * E * Iz / L_e^3; EIz_6 = 6 * E * Iz / L_e^2;
EIz_4  = 4 * E * Iz / L_e;   EIz_2 = 2 * E * Iz / L_e;
k_loc([2, 6, 8, 12], [2, 6, 8, 12]) = [ ...
    EIz_12,  EIz_6, -EIz_12,  EIz_6; ...
    EIz_6,   EIz_4,  -EIz_6,  EIz_2; ...
    -EIz_12, -EIz_6,  EIz_12, -EIz_6; ...
    EIz_6,   EIz_2,  -EIz_6,  EIz_4];

EIy_12 = 12 * E * Iy / L_e^3; EIy_6 = 6 * E * Iy / L_e^2;
EIy_4  = 4 * E * Iy / L_e;   EIy_2 = 2 * E * Iy / L_e;
k_loc([3, 5, 9, 11], [3, 5, 9, 11]) = [ ...
    EIy_12, -EIy_6, -EIy_12, -EIy_6; ...
    -EIy_6,   EIy_4,   EIy_6,   EIy_2; ...
    -EIy_12,  EIy_6,  EIy_12,  EIy_6; ...
    -EIy_6,   EIy_2,   EIy_6,   EIy_4];

% 3D Transformation Matrix T (12 x 12)
cx = dx / L_e; cy = dy / L_e; cz = dz / L_e;

if abs(cx) < 1e-5 && abs(cy) < 1e-5
    if cz > 0
        R_sub = [0 0 1; 0 1 0; -1 0 0];
    else
        R_sub = [0 0 -1; 0 1 0; 1 0 0];
    end
else
    D = sqrt(cx^2 + cy^2);
    R_sub = [cx, cy, cz; ...
        -cy/D, cx/D, 0; ...
        -cx*cz/D, -cy*cz/D, D];
end

T = blkdiag(R_sub, R_sub, R_sub, R_sub);

% Global Transformation
K_global = T' * k_loc * T;
end