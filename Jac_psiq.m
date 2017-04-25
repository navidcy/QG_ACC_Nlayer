function Jpsiqh = Jac_psiq(psih,qh)

global KX KY;

  u  = -real(ifft2(1i*KY.*psih));
  v  = +real(ifft2(1i*KX.*psih));
dqdx =  real(ifft2(1i*KX.*qh));
dqdy =  real(ifft2(1i*KY.*qh));

% UU=0*psih;
% UU(:,:,1)=Uj(1);UU(:,:,2)=Uj(2);UU(:,:,3)=Uj(3);
 
Jpsiq = u.*dqdx + v.*dqdy;

Jpsiqh = fft2(Jpsiq);


%   u  = -real(ifft2(1i*KY.*psih));
%   v  = +real(ifft2(1i*KX.*psih));  
%   q  =  real(ifft2(qh));

% Jpsiqh = 1i*KX.*fft2(u.*q) + 1i*KY.*fft2(v.*q);
% for in=1:3
%     u(:,:,in) = -real(ifft2(1i*KY.*psih(:,:,in)));
%     v(:,:,in) =  real(ifft2(1i*KX.*psih(:,:,in)));
%     q(:,:,in) =  real(ifft2(qh(:,:,in)));
% end
% for in=1:3
%     Jpsiqh(:,:,in) = 1i*KX.*fft2(u(:,:,in).*q(:,:,in)) + 1i*KY.*fft2(v(:,:,in).*q(:,:,in));
% end
% 
% dqdx =  real(ifft2(1i*KX.*qh));
% dqdy =  real(ifft2(1i*KY.*qh));
% 
% % UU=0*psih;
% % UU(:,:,1)=Uj(1);UU(:,:,2)=Uj(2);UU(:,:,3)=Uj(3);
%  
% Jpsiq = u.*dqdx + v.*dqdy;
% 
% Jpsiqh = fft2(Jpsiq);

