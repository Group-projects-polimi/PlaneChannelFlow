% Physical Parameters
rho = 1;  % Density kg/m^3
cp = 10;  % Specific heat capacity J/(K*kg)
k = 0.12;  % Thermal conductivity W/(m*K)
nu = 1e-2;  % Kinematic viscosity m^2/s
Gamma = k/(rho*cp);  % Heat diffusion coefficient m^2/s
L = 10;  % Channel length m
H = 1;  % Channel height m
Pe = 16.5;  % Peclet number
u_mean = Pe*Gamma/(2*H*rho);  % Average speed m/s
Re = 2*H*u_mean/nu;  % Reynolds number

% Mesh Parameters
nx = 50;  % Coarse mesh x number
ny = 5;  % Coarse mesh y number
dx = L/nx;
dy = H/ny;
N = nx*ny; % DOF number

% Boundary Temperature Parameters
T_in = 50;  % Inlet temperature °C
T_wall = 100;  % Wall temperature °C

% Velocity Field Function
u_x = @(y) 6*u_mean*(y/H)*(1-(y/H));

% Initialize Coefficient Matrix A and Right-hand Side Vector b
A = zeros(N,N);
b = zeros(N,1);

% Fill in the Equations for the Internal CV
for j = 2:(ny-1)
    for i = 2:(nx-1)
        n = (j-1)*nx+i;
        % Calculate the Advection-Diffusion Coefficient
        y_e = (j-1)*dy+dy/2;  % y_eastern
        ux_e = u_x(y_e);  % Mean velocity
        F_e = rho*ux_e*dy;  % Eastbound convection flux kg/s
        F_w = F_e;  % Incompressible flow
        % Diffusion coefficient
        D_e = Gamma*dy/dx;  % Eastern
        D_w = Gamma*dy/dx;  % Western
        D_n = Gamma*dx/dy;  % Northern
        D_s = Gamma*dx/dy;  % Southern
        % Coefficients of Algebraic Equations
        a_E = D_e+max(-F_e,0);
        a_W = D_w+max(F_w,0);
        a_N = D_n;
        a_S = D_s;
        a_P = a_E+a_W+a_N+a_S+(F_e-F_w);
        % Fill matrix A (The internal CV has no source term, and b(n) remains zero.)
        A(n,n) = a_P;
        A(n,n+1) = -a_E;
        A(n,n-1) = -a_W;
        A(n,n+nx) = -a_N;
        A(n,n-nx) = -a_S;
    end
end

% Fill boundary conditions at the inlet (Dirichlet T=T_in)
for j = 1:ny
    n = (j-1)*nx+1;
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_in;
end

% Fill the boundary conditions at the outlet (Neumann ∂T/∂x=0)
for j = 1:ny
    n = (j-1)*nx+nx;
    A(n,:) = 0;
    A(n,n) = 1;
    A(n,n-1) = -1;
    b(n) = 0;
end

% Fill boundary conditions at the walls including 4 corner CVs (Dirichlet T=T_wall)
% Lower wall surface
for i = 1:nx
    n = (1-1)*nx+i; % j=1
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_wall;
end
% Upper wall surface
for i = 1:nx
    n = (ny-1)*nx+i;
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_wall;
end

% Solve the linear equation system
T = A\b;

T_2D = reshape(T, nx, ny)';
x_coords = linspace(0, L, nx);   
y_coords = linspace(0, H, ny);    
[X, Y] = meshgrid(x_coords, y_coords);

 
%%% 12.
% Verify inlet temperature (plotting the temperature curve at x=0)
figure('Name', 'Inlet Temperature Validation');
plot(y_coords, T_2D(:, 1), 'b-o'); 
ylim([T_in-10, T_wall+10]);             
xlabel('y (m)');
ylabel('Temperature (°C)');
title('Temperature Distribution at the Inlet Boundary');
grid on;

% Verify wall surface temperature (plotting the temperature curve at y=0 and y=H)
figure('Name', 'Lower Wall Temperature Verification');
plot(x_coords, T_2D(1, :), 'r-o'); 
hold on
plot(x_coords, T_2D(end, :), 'b'); 
hold off
ylim([T_wall-1, T_wall+1]);        
xlabel('x (m)');
ylabel('Temperature (°C)');
legend({'Lower wall (y=0)', 'Upper wall (y=H)'}, 'Location', 'best');
title('Temperature distribution on the lower wall surface');
grid on;


%%% 13. 
% Visualization and Validation of Temperature Fields
figure('Name', 'Coarse-mesh temperature field');
pcolor(X, Y, T_2D);   
shading interp;         
colorbar;               
xlabel('x (m)');       
ylabel('y (m)');        
title(['Channel temperature field（nx=' num2str(nx) ', ny=' num2str(ny) '）']);
axis equal tight;

%visualizartion of the outlet temperature
T_o=T_2D(:,nx);
figure('Name', 'Visualization of the outlet temperature');
plot(y_coords, T_o, 'r-o');        
xlabel('y (m)');
ylabel('Temperature (°C)');
title('Temperature distribution on the outlet');
grid on;

%visualization of the centerline temperature including entry legth
T_c=T_2D(ceil(ny/2),:);
fully_dev_idx = find(T_c >= 0.9 * T_wall, 1, 'first');
x_e=x_coords(fully_dev_idx);

figure('Name', 'Visualization of the centerline temperature');
plot(x_coords, T_c, 'r-o'); 
xline(x_e, '--k', 'Fully developed region', ...
    'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'right');
xlabel('x (m)');
ylabel('Temperature (°C)');
title('Temperature distribution on the centerline');
grid on;

%velocity weighted mean temperature
T_weighted=zeros(nx,1);
dy = y_coords(2) - y_coords(1);
for i = 1:nx%implemetned using a general forumlation for calculating the integral via trapezoidal rule
    for j = 2:ny
        y=y_coords(j);
        y_prev=y_coords(j-1);
        dy = y_coords(j) - y_coords(j-1);
        T_weighted(i) = T_weighted(i) + (u_x(y_prev)*T_2D(j-1,i) + u_x(y)*T_2D(j,i))*dy/2;   
    end
    T_weighted(i) = T_weighted(i)/(u_mean*H);

end

figure('Name', 'Visualization of the velocity weighted mean temperature');
plot(x_coords, T_weighted, 'r-o'); 
xlabel('x (m)');
ylabel('Temperature (°C)');
title('Horizontal distribution of temperature weighted along y by u_x');
grid on;