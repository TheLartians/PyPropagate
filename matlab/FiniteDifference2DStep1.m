function u = FiniteDifference2DStep1(ra,rc,rfp,rf,up,u)
% Perform the first 2D Finite Difference step
% This solves the PDE: du/dz u = A d^2u/dx^2 + C d^2u/dy^2 + F u 
%
% Arguments:
% ra = A Dz/Dx^2, where Dz is the half step size in z direction and Dx is 
% the step size in x direction
% rc = C Dz/Dy^2, where Dy is the step size in y direction
% rfp = F Dz/2, where F is evaluated at z = z_i
% rf = F Dz/2, where F is evaluated at z = z_i+1
% up = u at z = z_i
% u = u at z = z_i+1 (only used for bounday condition)
%
% Returns:
% u at z = z_i+1

nx = size(u,1);
ny = size(u,2);

A = -rc * ones(ny-2,1);

for i=2:nx-1
    B = 1+2*rc-rf(i,2:ny-1);
    R = (up(i+1,2:ny-1)+up(i-1,2:ny-1))*ra + up(i,2:ny-1).*(1+rfp(i,2:ny-1)-2*ra);

    R(1) = R(1) + u(i,1) * rc;
    R(ny-2) = R(ny-2) + u(i,ny) * rc;

    u(i,2:ny-1) = Tridiagonal(A,B,A,R);
end


end
