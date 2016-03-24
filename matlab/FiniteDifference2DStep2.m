function u = FiniteDifference2DStep2(ra,rc,rfp,rf,up,u)
% Perform the second 2D Finite Difference step
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

u = transpose(FiniteDifference2DStep1(rc,ra,transpose(rfp),transpose(rf),transpose(up),transpose(u)));

end
