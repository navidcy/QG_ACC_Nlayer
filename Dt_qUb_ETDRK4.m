function [nlin_q,nlin_U] = Dt_qUb_ETDRK4(qh,Ub)

global FILTERsmooth
global KX ksq rek Uj etah topoflag taub Hj H rhoj rho forceqh Nz invS betaj
[Ny,Nx,Nz]=size(qh);

psih = tprod(invS,[1 2 3 -1],qh,[1 2 -1]);

% qhtopo = qh;

% if topoflag==1
%     qhtopo(:,:,Nz) = qhtopo(:,:,Nz) + etah;
% end

% Jpsiqh = Jac_psiq(psih,qhtopo);
% Jpsiqh = Jac_psiqUbUj(psih,qhtopo,Ub,Uj);

Jpsiqh = Jac_psiqUbUj(psih,qh,Ub,Uj);


nlin_q=0*qh;

for in=1:Nz
    nlin_q(:,:,in) = - Jpsiqh(:,:,in) - betaj(in)*1i*KX.*psih(:,:,in) ;%- (Ub+Uj(in))*1i*KX.*qhtopo(:,:,in);
end
psih_N=psih(:,:,Nz);

% wind stress at top layer
nlin_q(:,:,1) = nlin_q(:,:,1) + forceqh;

% Ekman drag at bottom layer
nlin_q(:,:,Nz) = nlin_q(:,:,Nz) - rek*(-ksq.*psih_N); 


formstress = real(mean(conj(psih_N(:)).*(1i*KX(:).*etah(:)))) / (Nx*Ny);
nlin_U = taub/(rhoj(1)*H) - Hj(Nz)/H * rek*Ub - Hj(Nz)/H * formstress;

% for in=1:Nz
%    nlin_q(:,:,in)=nlin_q(:,:,in).*FILTERsmooth;
% end


end