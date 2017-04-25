function [dqdt,psih] = Dt_q(qh,invS,FILTERsmooth,betaj)

global KX nu4 ksq rek Uj etah topoflag
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
    dqdt(:,:,in) = -Jpsiqh(:,:,in) -betaj(in)*1i*KX.*psih(:,:,in) - Uj(in)*1i*KX.*qh(:,:,in);
end

dqdt(:,:,3) = dqdt(:,:,3) - rek*(-ksq.*psih(:,:,3)); 