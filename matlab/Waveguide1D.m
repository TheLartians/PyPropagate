%% Units
keV = 0.016021766208;
hbarc = 0.00003161526;
um = 10;

%% Simulation box
Nx = 800;
Nz = 1000;

sx = 0.2*um;
sz = 800*um;

dx = sx/(Nx-1);
dz = sz/(Nz-1);

xindex = @(x) floor((x+sx/2)/dx+1);
zindex = @(z) floor((z+sz/2)/dz+1);

%% Parameters
E = 12 * keV; 
r_waveguide = 0.024 * um;

n_Ge = 1-6.006E-06+6.32E-07j;

n = ones(Nx,1) * n_Ge;
n(xindex(-r_waveguide):xindex(r_waveguide)) = 1;

%% Define wave equation
k = E/hbarc;
A = 1j/(2*k);
F = 1j*k/2*(n.^2 - 1);

%% Boundary Conditions 
u = ones(Nx,1);
u_boundary = @(zi) exp(F(1,1) * zi * dz);

%% Run simulation
ra = A * dz/dx^2;
rf = F * dz/2;

result = zeros(Nx,Nz);
result(:,1) = u;

for i=2:Nz
    up = u;
    u([1,Nx]) = u_boundary(i-1);
    u = FiniteDifference1D(ra,rf,rf,up,u);
    result(:,i) = u;
end

ma%% Plot the result

image(abs(result).^2,'CDataMapping','scaled')
colorbar

