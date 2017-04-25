function [dqdt,dUbdt,psih,formstress] = Dt_qUb(qh,Ub,invS,betaj)

global KX nu4 ksq rek Uj etah topoflag taub Hj H rhoj forceqh Nz
[Ny,Nx,Nz]=size(qh);

% tic
% psih=PVinverse(qh,invS);
% toc

qht = qh;
if topoflag==1
    qht(:,:,Nz) = qht(:,:,Nz) - etah;
end

psih = tprod(invS,[1 2 3 -1],qht,[1 2 -1]);

% norm(psih(:)-psih2(:))/norm(psih(:))

Jpsiqh = Jac_psiq(psih,qh);

dqdt=0*qh;

for in=1:Nz
    dqdt(:,:,in) = -Jpsiqh(:,:,in) -betaj(in)*1i*KX.*psih(:,:,in) - (Ub+Uj(in))*1i*KX.*qh(:,:,in);
end

% wind stress at top layer
dqdt(:,:,1) = dqdt(:,:,1) + forceqh;


psih_N=psih(:,:,Nz);

% Ekman drag at bottom layer
dqdt(:,:,Nz) = dqdt(:,:,Nz) - rek*(-ksq.*psih_N); 


formstress = real(mean(conj(psih_N(:)).*(1i*KX(:).*etah(:)))) / (Nx*Ny);

dUbdt = taub/(rhoj(1)*Hj(1)) - rek*Ub - Hj(Nz)/H * formstress;