clear all; close all; clc
% Symbolic Solver
syms delta P1 P2

%% Inputs
% Keep units consistent
Inputs = load('reticulate.mat');
% P      = Inputs.model.P;
% M      = Inputs.model.M;
% N      = Inputs.model.N;
% E      = Inputs.model.E; % Does not include 1st column (element number)
% BC     = Inputs.model.BC;
% BC_NH  = Inputs.model.BC_NH;
% F      = Inputs.model.F;
Yield  = 50e3;
y      = 1;


% Property Matrix = property number, material number, area, inertia
P = [1, 1, 10, 200;
     2, 1, 10, 100];

 
% Material Matrix = material number, modulus of elasticity
M = [1, 30e6]; 


% Node Matrix
% 1st column = node number, 2nd = x-component, 3rd = y-component
N = [1, -60, 0;
     2, -60, 120;
     3, 60, 120;
     4, 60, 0];


% Element Connectivity Matrix
% 1st column = element number, 2nd = node 1, 3rd = node 2, 4th = property
E = [1, 1, 2, 1;
     2, 2, 3, 2;
     3, 3, 4, 1];


% Boundary Conditions
% 1st column = node number, 2nd = x, 3rd = y, 4th = rotation (1 = fixed)
BC = [1, 1, 1, 1;
      4, 1, 1, 1];
  
  
% Displacement Input 
% 1st column = node number, 2nd = x, 3rd = y, 4th = rotation
BC_NH=[];


% Force Input
% 1st column = node number, 2nd = x, 3rd = y, 4th = rotation, 5th = value
F = [2, 1, 0, 0, 10000;
     3, 0, 0, 1, 5000];

        
%% FEA Solver
Number_of_Nodes = size(N,1);
Number_of_Elements = size(E,1);

DOF = 3; % Degrees of freedom
DOFS = Number_of_Nodes*DOF; % Total degrees of freedom

Connectivity = cell(1,Number_of_Elements);
K_Global = zeros(Number_of_Nodes*DOF); % add sym() if symbolic

count_i = 0;
count_j = 0;
for i = 1:Number_of_Elements
    % Find Elemental Nodes and Plot
    Node_1 = E(i,2);
    Node_2 = E(i,3);
%   Node_1 = E(i,1);
%   Node_2 = E(i,2);
    N1 = N(Node_1,2:3); % x and y of first node
    N2 = N(Node_2,2:3); % x and y of second node
    plot([N1(1);N2(1)],[N1(2);N2(2)],'ko--'); hold on
   
    % Find Element Length and Angle
    V = [N2(1)-N1(1),N2(2)-N1(2)];
    L = norm(V);
    D = atand(V(2)/V(1));
    S = sind(D);
    C = cosd(D);
    
    % Find Properties
    propID      =   E(i,end);
    propMat     =   P(P(:,1) == propID,2);
    A           =   P(P(:,1) == propID,3);
    I           =   P(P(:,1) == propID,4);
    Elastic     =   M(M(:,1) == propMat,2);
    
    % Local Stiffness Matrix (Global Coordinates)
    K_Element_Global = (Elastic/L).*...
        [A*C^2+12*I*S^2/L^2, (A-12*I/L^2)*C*S, -6*I*S/L, ...
        -(A*C^2+12*I*S^2/L^2), -(A-12*I/L^2)*C*S, -6*I*S/L;
        (A-12*I/L^2)*C*S, A*S^2+12*I*C^2/L^2, 6*I*C/L, ...
        -(A-12*I/L^2)*C*S, -(A*S^2+12*I*C^2/L^2), 6*I*C/L;
        -6*I*S/L, 6*I*C/L, 4*I, 6*I*S/L, -6*I*C/L, 2*I;
        -(A*C^2+12*I*S^2/L^2), -(A-12*I/L^2)*C*S, 6*I*S/L, ...
        A*C^2+12*I*S^2/L^2, (A-12*I/L^2)*C*S, 6*I*S/L;
        -(A-12*I/L^2)*C*S, -(A*S^2+12*I*C^2/L^2), -6*I*C/L, ...
        (A-12*I/L^2)*C*S, A*S^2+12*I*C^2/L^2, -6*I*C/L;
        -6*I*S/L, 6*I*C/L, 2*I, 6*I*S/L, -6*I*C/L, 4*I];
    
    % Cell Array stores row/column numbers
    Connectivity{i} = [Node_1*DOF-2,...
                       Node_1*DOF-1,...
                       Node_1*DOF,...
                       Node_2*DOF-2,...
                       Node_2*DOF-1,...
                       Node_2*DOF];             
    
    % Construct Global Stiffnes Matrix
    for ii=Connectivity{i}
        count_i=count_i+1;
        for jj=Connectivity{i}
             count_j=count_j+1;
             K_Global(ii,jj) = K_Global(ii,jj)+...
                 K_Element_Global(count_i,count_j);
        end
        count_j=0;
    end
    count_i=0;
end


% Force Vector
Force = zeros(DOFS,1); % add sym() if symbolic
for i=1:size(F,1)
    node_f = F(i,1);
    if F(i,2) == 1 % x-component
        Force(node_f*DOF-DOF+1) = F(i,end);
    end
    if F(i,3) == 1 % y-component
        Force(node_f*DOF-DOF+2) = F(i,end);
    end
    if F(i,4) == 1 % moment
        Force(node_f*DOF-DOF+3) = F(i,end);
    end
end
% Non-homogeneous BC Vector
Force_NH = zeros(DOFS,1); % add sym() if symbolic
for i=1:size(BC_NH,1)
    node_d = BC_NH(i,1);
    if BC_NH(i,2) ~= 0 % x-component
        Force_NH(node_d*DOF-DOF+1) = BC_NH(i,2);
    end
    if BC_NH(i,3) ~= 0 % y-component
        Force_NH(node_d*DOF-DOF+2) = BC_NH(i,3);
    end
    if BC_NH(i,4) ~= 0 % rotation
        Force_NH(node_d*DOF-DOF+3) = BC_NH(i,4);
    end
end
% Update Force Vector
for i=(find(Force_NH)).'
    Force = Force-K_Global(:,i).*Force_NH(i);
end
 
 
% Row Vector (if non-zero then eliminate DOF)
BC_eliminate = zeros(DOFS,1);
counter = 0;
for i = (BC(:,1)).'
    counter = counter+1;
    BC_eliminate(i*DOF) = BC(counter,4);
    BC_eliminate(i*DOF-1) = BC(counter,3);
    BC_eliminate(i*DOF-2) = BC(counter,2);
end
% Find where F_nh is non-zero and eliminate DOF
Force_nh_binary = Force_NH;
mask = Force_nh_binary~=0;
Force_nh_binary(mask) = 1;
BC_eliminate = BC_eliminate+Force_nh_binary;
BC_eliminate(all(~K_Global,2),:) = 1;


% Reduce Stiffness and Force Matrices
rows_eliminate = find(BC_eliminate)'; % Get index of rows to eliminate
rows_keep = find(~BC_eliminate)'; % Get index of rows to keep
K_Global_Reduced = K_Global;
K_Global_Reduced(rows_eliminate,:) = [];
K_Global_Reduced(:,rows_eliminate) = [];
Force_Reduced = Force;
Force_Reduced(rows_eliminate) = [];


% Find Displacements
Displacement_Reduced = K_Global_Reduced\Force_Reduced;
Displacement(rows_keep) = Displacement_Reduced;
Displacement(find(Force_NH)) = nonzeros(Force_NH);
Displacement(rows_eliminate) = 0;
Displacement = Displacement.';


% Solve for Reaction Forces
Force_Reaction = K_Global*Displacement;


% Element Forces and Stresses
Element_Forces = zeros(2*DOF,Number_of_Elements);
Stresses = zeros(2*DOF,Number_of_Elements);
for i = 1:Number_of_Elements
    % Find Elemental Nodes
    Node_1 = E(i,2);
    Node_2 = E(i,3);
%   Node_1 = E(i,1);
%   Node_2 = E(i,2);

    N1 = N(Node_1,2:3); % x and y of first node
    N2 = N(Node_2,2:3); % x and y of second node
   
    % Find Element Length and Angle
    V = [N2(1)-N1(1),N2(2)-N1(2)];
    L = norm(V);
    D = atand(V(2)/V(1));
    S = sind(D);
    C = cosd(D);
    
    % Find Properties
    propID      =   E(i,end);
    propMat     =   P(P(:,1) == propID,2);
    A           =   P(P(:,1) == propID,3);
    I           =   P(P(:,1) == propID,4);
    Elastic     =   M(M(:,1) == propMat,2);
    
    % Get Forces and Stresses
    C1 = A*Elastic/L;
    C2 = Elastic*I/L^3;
    T = [C,S,0,0,0,0;
         -S,C,0,0,0,0;
         0,0,1,0,0,0;
         0,0,0,C,S,0;
         0,0,0,-S,C,0;
         0,0,0,0,0,1];
    K_Element_Local = [C1, 0, 0, -C1, 0, 0;
                       0, 12*C2, 6*C2*L, 0, -12*C2, 6*C2*L;
                       0, 6*C2*L, 4*C2*L^2, 0, -6*C2*L, 2*C2*L^2;
                       -C1, 0, 0, C1, 0, 0;
                       0, -12*C2, -6*C2*L, 0, 12*C2, -6*C2*L;
                       0, 6*C2*L, 2*C2*L^2, 0, -6*C2*L, 4*C2*L^2];
    u = [Displacement(Node_1*DOF-2);
         Displacement(Node_1*DOF-1);
         Displacement(Node_1*DOF);
         Displacement(Node_2*DOF-2);
         Displacement(Node_2*DOF-1);
         Displacement(Node_2*DOF)];
    
    Element_Forces(:,i) = K_Element_Local*T*u;
    % Normal Stresses
    Stresses(1:2,i) = Element_Forces(1:2,i)/A;
    Stresses(4:5,i) = Element_Forces(4:5,i)/A;
    % Bending Stresses
    Stresses(3,i) = Element_Forces(3,i)*y/I;
    Stresses(6,i) = Element_Forces(6,i)*y/I;
end


% Maxiumums and Factor of Safety
[max_stress_values, stress_indices] = max(abs(Stresses));
[max_stress_value, stress_index] = max(max_stress_values);
Max_Stress = Stresses(stress_indices(stress_index),stress_index);
Max_Displacement =  max(abs(Displacement));
Safety_Factor = min(Yield/abs(Max_Stress));


% Plot New Displacement
Magnification = 100;
for i = 1:Number_of_Elements
    Node_1 = E(i,2);
    Node_2 = E(i,3);
%   Node_1 = E(i,1);
%   Node_2 = E(i,2);
    N1_x_new = N(Node_1,2)+Magnification*Displacement(Node_1*DOF-2);
    N1_y_new = N(Node_1,3)+Magnification*Displacement(Node_1*DOF-1);
    N2_x_new = N(Node_2,2)+Magnification*Displacement(Node_2*DOF-2);
    N2_y_new = N(Node_2,3)+Magnification*Displacement(Node_2*DOF-1);
    plot([N1_x_new;N2_x_new],[N1_y_new;N2_y_new],'ro-');hold on
    
    % Add Node and Element Numbers
    midline_x = N1_x_new+(N2_x_new-N1_x_new)/2;
    midline_y = N1_y_new+(N2_y_new-N1_y_new)/2;
    text(midline_x,midline_y,num2str(i),'Fontsize',12); hold on
    text(N1_x_new,N1_y_new,num2str(Node_1),'Fontsize',12,...
        'VerticalAlignment','Bottom');hold on
    text(N2_x_new,N2_y_new,num2str(Node_2),'Fontsize',12,...
        'VerticalAlignment','Bottom');hold on
end
grid on; grid minor
xlabel('X-Direction','Fontsize',12,'Interp','Latex')
ylabel('Y-Direction','Fontsize',12,'Interp','Latex')
title(['Undeformed \& Deformed Plot (Mag  x',...
    num2str(Magnification),')'],'Interp','Latex','Fontsize',14)