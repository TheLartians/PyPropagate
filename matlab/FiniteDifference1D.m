function u = FiniteDifference1D(ra,rfp,rf,up,u)

% Perform a 1D Finite Difference step
% This solves the PDE: du/dz u = A du/dx + F u 
%
% Arguments:
% ra = A Dz/Dx^2, where Dz is the step size in z direction and Dx is the
% step size in x direction
% rfp = F Dz/2, where F is evaluated at z = z_i
% rf = F Dz/2, where F is evaluated at z = z_i+1 
% up = u at z = z_i
% u = u at z = z_i+1 (only used for bounday condition)
%
% Returns:
% u at z = z_i+1

n = length(u);

A = -ra/2 * ones(n-2,1);
B = 1+ra-rf(2:n-1);
R = (up(3:n)+up(1:n-2))*ra/2 + up(2:n-1).*(1+rfp(2:n-1)-ra);

R(1) = R(1) + u(1) * ra/2.;
R(n-2) = R(n-2) + u(n) * ra/2;

u(2:n-1) = Tridiagonal(A,B,A,R);

