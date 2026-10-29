% =========================================================================
% CVL7411 ADVANCED STRUCTURAL ANALYSIS - PROGRAMMING ASSIGNMENT GROUP (D)
% 3D SINGLE-STOREY HEXAGONAL SPACE FRAME WITH RIGID SLAB
% =========================================================================

function [coords, U_nodes, R_nodes, F_members] = analyze_space_frame(...
    L, H, Px, E, G, bc, hc, bb, hb)

% Columns - Cross Sectional Properties
A_c  = bc * hc;
Iy_c = (bc * hc^3) / 12;
Iz_c = (hc * bc^3) / 12;
a_c  = max(bc, hc); b_c = min(bc, hc);
J_c  = a_c * (b_c^3) * (1/3 - 0.21 * (b_c/a_c) * (1 - (b_c^4)/(12*a_c^4)));

% Beams with Rigid Slab Diaphragm Constraint 
rigid_factor = 1e6;
A_b  = (bb * hb) * rigid_factor; 
Iy_b = (bb * hb^3) / 12;
Iz_b = ((hb * bb^3) / 12) * rigid_factor; 
a_b  = max(bb, hb); b_b = min(bb, hb);
J_b  = a_b * (b_b^3) * (1/3 - 0.21 * (b_b/a_b) * (1 - (b_b^4)/(12*a_b^4)));

% Geometry and its connectivity is specified
[coords, members] = generate_geometry(L, H);

% Assembling Global Stiffness Matrix 
K_global = zeros(72, 72); % 12 x 6 = 72 DOFs
for m = 1:12
    n1 = members(m, 1); n2 = members(m, 2);
    c1 = coords(n1, :); c2 = coords(n2, :);

    if m <= 6
        [Ke, ~, ~] = element_stiffness(c1, c2, E, G, A_c, Iy_c, Iz_c, J_c);
    else
        [Ke, ~, ~] = element_stiffness(c1, c2, E, G, A_b, Iy_b, Iz_b, J_b);
    end
    dofs = [(n1-1)*6+(1:6), (n2-1)*6+(1:6)];
    K_global(dofs, dofs) = K_global(dofs, dofs) + Ke;
end

% Applying Boundary Conditions & Load Vector
fixed_dofs = 1:36;  % Base nodes 1-6 fixed (36 DOFs)
free_dofs  = 37:72; % Roof nodes 7-12 free (36 DOFs)

F_global = zeros(72, 1);
for node = 7:12
    dof_x = (node - 1) * 6 + 1;
    F_global(dof_x) = Px; % Lateral load Px in +X direction at roof level
end

% Solving Free DOFs via Custom Gaussian Elimination 
K_ff = K_global(free_dofs, free_dofs);
F_f  = F_global(free_dofs);
U_f  = custom_gauss_solver(K_ff, F_f);

U_sys = zeros(72, 1);
U_sys(free_dofs) = U_f;

% Computation of Base Support Reactions
K_rf = K_global(fixed_dofs, free_dofs);
R_fixed = K_rf * U_f;

U_nodes = reshape(U_sys, 6, 12)'; % [12 x 6] Displacements
R_nodes = reshape(R_fixed, 6, 6)'; % [6 x 6] Reactions at Nodes 1-6

% Computation of Member Internal End Actions 
F_members = zeros(12, 12);
for m = 1:12
    n1 = members(m, 1); n2 = members(m, 2);
    dofs = [(n1-1)*6+(1:6), (n2-1)*6+(1:6)];
    c1 = coords(n1, :); c2 = coords(n2, :);

    if m <= 6
        [~, T, k_loc] = element_stiffness(c1, c2, E, G, A_c, Iy_c, Iz_c, J_c);
    else
        [~, T, k_loc] = element_stiffness(c1, c2, E, G, A_b, Iy_b, Iz_b, J_b);
    end

    u_global = U_sys(dofs);
    u_local  = T * u_global;
    f_local  = k_loc * u_local; 
    F_members(m, :) = f_local';
end
end