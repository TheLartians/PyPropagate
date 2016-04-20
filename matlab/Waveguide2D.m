%% Units
keV = 0.016021766208;
hbarc = 0.00003161526;
um = 10;

%% Simulation box
Nx = 500;
Ny = 500;
Nz = 1000;

sx = 0.2*um;
sy = 0.2*um;
sz = 800*um;

dx = sx/(Nx-1);
dy = sy/(Ny-1);
dz = sz/(Nz-1);

xindex = @(x) floor((x+sx/2)/dx+1);
yindex = @(y) floor((y+sy/2)/dy+1);
zindex = @(z) floor((z+sz/2)/dz+1);

%% Parameters
E = 12 * keV; 
r_waveguide = 0.024 * um;

n_Ge = 1 - 6.006E-06 + 6.32E-07j;

n = ones(Nx,Ny) * n_Ge;
n(xindex(-r_waveguide):xindex(r_waveguide),yindex(-r_waveguide):yindex(r_waveguide)) = 1;

%% Define wave equation
k = E/hbarc;
A = 1j/(2*k);
F = 1j*k/2*(n.^2 - 1);

%% Boundary Conditions 
u = ones(Nx,Ny);
u_boundary = @(zi) exp(F(1,1) * zi * dz);

%% Run simulation
ra = A * dz/dx^2/2;
rc = A * dz/dy^2/2;
rf = F * dz/2/2;

result = zeros(Nx,Nz);
result(:,1) = u(:,yindex(0));

for i=2:Nz
    up = u;
    u(1,:) = u_boundary(i-1.5);
    u(Nx,:) = u_boundary(i-1.5);
    u(:,1) = u_boundary(i-1.5);
    u(:,Ny) = u_boundary(i-1.5);
    
    u = FiniteDifference2DStep1(ra,rc,rf,rf,up,u);
    
    up = u;
    u(1,:) = u_boundary(i-1);
    u(Nx,:) = u_boundary(i-1);
    u(:,1) = u_boundary(i-1);
    u(:,Ny) = u_boundary(i-1);

    u = FiniteDifference2DStep2(ra,rc,rf,rf,up,u);
    
    result(:,i) = u(:,yindex(0));
    fprintf('\rStep %i of %i     ',i,Nz);
end

%% Plot the result
image(abs(result).^2,'CDataMapping','scaled')
colorbar

