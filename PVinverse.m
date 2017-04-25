function psih = PVinverse(qh,invS)

[Ny,Nx,Nz]=size(qh);
psih=0*qh;
for il=1:Ny
    for ik=1:Nx
        qtemp = squeeze(qh(il,ik,:));
        invs=squeeze(invS(il,ik,:,:));
        psitemp = invs*qtemp;
        psih(il,ik,:)=psitemp;
    end
end
