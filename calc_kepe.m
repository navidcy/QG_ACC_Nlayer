function [ke,pe,KEb]=calc_kepe(psih,Ub)

global ksq Hj Nz H Nx Ny greduced f0

ke=zeros(Nz,1);pe=zeros(Nz-1,1);

for in=1:Nz
    ke_n = .5*(ksq.*abs(psih(:,:,in)).^2);
    ke(in) = Hj(in)/H * sum(ke_n(:))/(Nx*Ny)^2;
end


for in=1:Nz-1
    pe_n05 = .5*f0^2./greduced(in)*(abs(psih(:,:,in+1)-psih(:,:,in)).^2);
    pe(in) = sum(pe_n05(:))/(Nx*Ny)^2;
end

KEb = .5*Ub.^2;
