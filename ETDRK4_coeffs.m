function [eL,eL2,Q,fu,fab,fc] = ETDRK4_coeffs(nu,h,ksq,dt,Nx,Ny)

% calculation of the EDTRK4 coefficients

 L  = - nu*(ksq.^h);
 
eL  = exp(L*dt);
eL2 = exp(L*dt/2);

M=32;

r  = 1*exp(2i*pi*((1:M)/M));
fu = zeros(Ny,Nx); fab = fu; fc = fu; Q = fu;

for j = 1:M
     z = r(j) + L*dt;
     Q =   Q + dt*( exp(z/2)-1 )./z;
    fu =  fu + dt*( -4 -z +exp(z).*(4 -3*z +z.^2) )./z.^3;
   fab = fab + dt*( +2 +z +exp(z).*(-2 +z) )./z.^3;
    fc =  fc + dt*( -4 -3*z -z.^2 +exp(z).*(4 -z) )./z.^3;
end

fu = (fu/M); fab = (fab/M); fc = (fc/M); Q = (Q/M);

% eL  = repmat(eL ,[1 1 Nz]);
% eL2 = repmat(eL2,[1 1 Nz]);
% fu  = repmat(fu ,[1 1 Nz]);
% fab = repmat(fab,[1 1 Nz]);
% fc  = repmat(fc ,[1 1 Nz]);
% Q   = repmat(Q  ,[1 1 Nz]);