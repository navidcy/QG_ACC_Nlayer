function Jpsiqh = Jac_psiqUbUj(psih,qh,Ub,Uj)

global KX KY Nz eta topoflag;

%   u  = -real(ifft2(1i*KY.*psih));
%   v  = +real(ifft2(1i*KX.*psih));
% dqdx =  real(ifft2(1i*KX.*qh));
% dqdy =  real(ifft2(1i*KY.*qh));
% 
% % UU=0*psih;
% % UU(:,:,1)=Uj(1);UU(:,:,2)=Uj(2);UU(:,:,3)=Uj(3);
%  
% Jpsiq = u.*dqdx + v.*dqdy;
% 
% Jpsiqh = fft2(Jpsiq);


  u  = -real(ifft2(1i*KY.*psih));
  v  = +real(ifft2(1i*KX.*psih));
  
  q  =  real(ifft2(qh));

  qtopo = q;
if topoflag==1
    qtopo(:,:,Nz) = qtopo(:,:,Nz) + eta;
end

% Jpsiqh = 1i*KX.*fft2(u.*q) + 1i*KY.*fft2(v.*q);

Jpsiqh=0*qh;

for in=1:Nz
    Jpsiqh(:,:,in) = 1i*KX.*fft2((Ub+Uj(in)+u(:,:,in)).*qtopo(:,:,in)) + 1i*KY.*fft2(v(:,:,in).*qtopo(:,:,in));
end

% dqdx =  real(ifft2(1i*KX.*qh));
% dqdy =  real(ifft2(1i*KY.*qh));
% 
% % UU=0*psih;
% % UU(:,:,1)=Uj(1);UU(:,:,2)=Uj(2);UU(:,:,3)=Uj(3);
%  
% Jpsiq = u.*dqdx + v.*dqdy;
% 
% Jpsiqh = fft2(Jpsiq);

