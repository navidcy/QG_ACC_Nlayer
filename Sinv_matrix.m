function invS = Sinv_matrix(Fm,Fp,ksq)

% This function calculates S^{-1} that is needed to invert the QGPV to get
% the streamfunction at each layer. The matrix S an (Nz x Nz) matrix that
% relates the QGPV vector q with elements the (kx,ky) spectral components
% of q at each layer with the streamfunction vector \psi with elements the
% (kx,ky) spectral components of \psi at each layer:
%
% \hat{q}_j = S_{jk}\hat{\psi}_k


Nz=length(Fm)+1;
[Ny,Nx]=size(ksq);

[a,b]=size(Fm);
if b>a
    Fm=Fm.';
end
[a,b]=size(Fp);
if b>a
    Fp=Fp.';
end

if length(Fm)~=Nz-1
    error('Fm');
end
if length(Fp)~=Nz-1
    error('Fp');
end

invS = zeros(Ny,Nx,Nz,Nz);
I = eye(Nz);
for il=1:Ny
    for ik=1:Nx
        k2 = ksq(il,ik);
        if k2==0
            k2=1;
        end
        F = +diag(Fp,1) + diag(Fm,-1) - diag([Fp;0]+[0;Fm],0);
        S = - k2*I + F ;
        
        invS(il,ik,:,:)=I/S;
    end
end

invS(1,1,:,:)=0;