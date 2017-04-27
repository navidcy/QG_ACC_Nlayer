clearvars -global;
clearvars;

global KX KY Uj ksq rek etah Hj H Nx Ny Nz topoflag greduced f0 taub forceqh rhoj invS betaj nu h beta
global FILTERsmooth eta rho

cfllimit=0.8;
figNo=10;

nofilter=1;

Lx=2*pi;Ly=2*pi;

Lx=1e6; %m
Ly=1e6; %m
Ld = 15000; %m

topoflag=1;
figure(334);clf;

nu=1e9; h=2;

rek = 1e-7;

Nx=32; Ny=Nx; Nz=2;
dx=Lx/Nx;dy=Ly/Ny;
x=0:dx:Lx-dx;x=x+dx/2;
y=0:dy:Ly-dy;y=y+dy/2;
[X,Y]=meshgrid(x,y);


Hj=[500;1750;1750]; %m
Hj=[500;1750]; %m
% Hj=4000;
% Hj=[500;600;600;600;1200;1200];
H = sum(Hj(:));
rhoj=[1025;1025.275;1025.640]; %kg*m^(-3)
rhoj=[1025;1025.275]; %kg*m^(-3)
% rhoj=1035;
% rhoj=[1025;1025.275;1025.5;1026;1026.5;1027]; %kg*m^(-3)
% rhoj=[1025;1028;1035]; %kg*m^(-3)
rho = sum(Hj.*rhoj)/H;
Uj = [0.05;0.025;0.001]*1;
Uj=0;

Uj=[zeros(Nz,1)];

% Hj = [500:500:5000]';
% % rhoj = [1025:-.5:1020.5];
% rhoj = [1025:-.5:1020.5]';
% Uj = [.05:-.005:.005]';

beta=1*1.2130692965249345e-11; % m^(-1)s^(-1)
f0 = 0.0001236812857687059; % s^(-1)
g = 9.81; % m*s^(-2)

hrms=1500; % m
etarms = hrms*f0/(Hj(Nz));



if length(Hj)~=Nz,error('Hj');end
if length(rhoj)~=Nz,error('rhoj');end


greduced = g*(rhoj(2:end)-rhoj(1:end-1))./rhoj(1:end-1); %definition to match PYQG

% greduced = g*(rhoj(2:end)-rhoj(1:end-1))./rhoj(2:end); %CORRECT DEFINITION
% format long
% greduced
% format short



Fm = f0^2./( greduced .* Hj(2:end  ) ); % m^(-2)
Fp = f0^2./( greduced .* Hj(1:end-1) ); % m^(-2)

if Nz>1
    betaj=zeros(Nz,1);
    for in=1:Nz
        if in==1
            betaj(in) = beta - Fp(in)*(Uj(in+1)-Uj(in));
        elseif in==Nz
            betaj(in) = beta - Fm(in-1)*(Uj(in-1)-Uj(in));
        else
            betaj(in) = beta - Fm(in-1)*(Uj(in-1)-Uj(in)) - Fp(in)*(Uj(in+1)-Uj(in));
        end
    end
else
    betaj=beta;
end

%         
% betaj2 = [beta                       - Fp(1)*(Uj(2)-Uj(1));
%          beta - Fm(1)*(Uj(1)-Uj(2)) - Fp(2)*(Uj(3)-Uj(2));
%          beta - Fm(2)*(Uj(2)-Uj(3))                       ];




kx=2*pi/Lx*[0:Nx/2-1 -Nx/2:-1];
ky=2*pi/Ly*[0:Ny/2-1 -Ny/2:-1];

[KX,KY]=meshgrid(kx,ky);
ksq=KX.^2+KY.^2;


eta = 2*etarms*cos(5*2*pi/Lx*X).*cos(5*2*pi/Ly*Y);
eta = sqrt(2)*etarms*cos(8*2*pi/Lx*X);
etah = fft2(eta);
etah = (randn(Ny,Nx)+1i*randn(Ny,Nx)).*(ksq>(2*pi/Lx*5)^2&ksq<(2*pi/Lx*7)^2);

eta = real(ifft2(etah));
% eta = sqrt(2)*etarms*cos(5*2*pi/Lx*X);

% eta = 2*etarms*cos(5*2*pi/Lx*X).*cos(5*2*pi/Ly*Y);
eta = etarms*eta/sqrt(mean(eta(:).^2));
etah=fft2(eta);

% figure(1);clf;
% pcolor2(X,Y,eta)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HIGH-WAVENUMBER FILTER
% a is chosen so that the energy at the largest nondim
% wavenumber K*dx be zero whithin machine double precision
s=4;
Kmax = (2*pi/Lx)*Ny/2; Kmax_s=Kmax*dy;
kcut = 2/3*Kmax;kcut_s=kcut*dy;
a = -log(1e-15)/(Kmax_s-kcut_s)^s * dy^s;
K=sqrt(KX.^2+KY.^2);

FILTERsmooth2 = 1*ones(Ny,Nx).*abs(K<=(2*pi/Lx)*Ny/3) + exp(-a*(K-kcut).^s).*abs(K>(2*pi/Lx)*Ny/3);

cphi=.65*pi;
filterfac=23.6;
wv = sqrt((KX*dx).^2+(KY*dy).^2);
FILTERsmooth = exp(-filterfac*(wv-cphi).^4);
FILTERsmooth(wv<=cphi)=1;


% FILTER(KX==0&KY==0)=0;
if nofilter==1
    FILTER = ones(Ny,Nx);
else
    FILTER = FILTERsmooth;
end;


% FILTERsmooth  = ones(size(FILTERsmooth));
% FILTERsmooth2 = ones(size(FILTERsmooth2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% figure(1)
% subplot(131)
% pcolor2(fftshift(KX),fftshift(KY),fftshift(FILTERsmooth))
% axis square;
% subplot(132)
% pcolor2(fftshift(KX),fftshift(KY),fftshift(FILTERsmooth2))
% axis square;
% subplot(133)
% pcolor2(fftshift(KX),fftshift(KY),fftshift(abs(FILTERsmooth2-FILTERsmooth)))
% axis square;

% F = +diag(Fp,1) + diag(Fm,-1) - diag([Fp;0]+[0;Fm],0);
% num2str(F,'%1.5e')

invS = Sinv_matrix(Fm,Fp,ksq);

xi=randn(Ny,Nx,Nz)+1i*randn(Ny,Nx,Nz);
xi(1,1,:)=0*xi(1,1,:);
q = real(ifft2(1e-7*xi.*FILTERsmooth));
q = q-mean(q(:));
q=0*q;
  
% for in=1:Nz
%     q(:,:,in) = 1e-9*( cos(2*pi/Lx*3*X).*cos(2*pi/Ly*2*Y) + 0*.7*cos(2*pi/Lx*6*X).*cos(2*pi/Ly*6*Y) - 1*.4*cos(2*pi/Lx*3*X).*cos(2*pi/Ly*4*Y) );
% %     q(:,:,in) = 1e-7*( cos(2*pi/Lx*1*X).*cos(2*pi/Ly*1*Y));
% end
qimatlab=q;


Ub=0;

qh = fft2(q);
qih=qh;
% tic
% psih=PVinverse(qh,invS);
% 
% toc
% tic
% psih2 = tprod(invS,[1 2 3 -1],qh,[1 2 -1]);
% toc


% Jpsiqh = Jac_psiq(psih,qh);
% 
% Jpsiq = real(ifft2(Jpsiqh));
% 
% figure(32)
% for in=1:Nz
%     subplot(1,Nz,in);
%     pcolor2(X,Y,Jpsiq(:,:,in));shading interp;
% end

dt=1000;
% dt=1500;

 tau = 1*.25*sin(2*pi/Ly*(Y-Ly/2)); % N*m^(-2)
taub = 13*.05; % N*m^(-2)

forceqh = -1i*KY.*fft2(tau) / (rhoj(1)*Hj(1));
forceqh=1*forceqh;

Nt=2000000

KE   = zeros(Nz,Nt);
KEUb = zeros(1,Nt);
if Nz>1
    PE   = zeros(Nz-1,Nt);
end

Ubt = zeros(1,Nt);
formstress = zeros(1,Nt);


tic
[eL,eL2,Q,fu,fab,fc] = ETDRK4_coeffs(nu,h,ksq,dt,Nx,Ny);
toc


Ub=0;
% 
% % 
% % [nlin_q,nlin_U] = Dt_qUb_ETDRK4(qh,0);
% % Jh=nlin_q;
% % Jh(11,11,1)
% 
% 
% 
% 
psih = tprod(invS,[1 2 3 -1],qh,[1 2 -1]);
% 
% % psih = PVinverse(qh,invS);
% 
jacobmatlab=real(ifft2(Jac_psiqUbUj(psih,qh,0*Ub,0*Uj)));
% 
% 
% u = real(ifft2(-1i*KY.*psih));
% v = real(ifft2(+1i*KX.*psih));
% 
% JAC = real(ifft2(1i*KX.*fft2(u.*q)) + ifft2(1i*KY.*fft2(v.*q)));
% JAC = u.*real(ifft2(1i*KX.*qh)) + v.*real(ifft2(1i*KY.*qh));
% 
% 
% % figure(1);
% % pcolor2(JAC(:,:,1));colorbar;
% 
% 
% 
% Jh=Jac_psiq(psih,qh);
% for in=1:Nz
%     jacobmatlab(:,:,in)=real(ifft2(Jh(:,:,in)));
% end
% 
% JAC = real(ifft2(1i*KX.*fft2(u.*q)) + ifft2(1i*KY.*fft2(v.*q)));
% JAC = u.*real(ifft2(1i*KX.*qh)) + v.*real(ifft2(1i*KY.*qh));
% 
% jacobmatlab = JAC;
% jacobmatlab=real(ifft2(Jac_psiqUbUj(psih,qh,0*Ub,0*Uj)));
% 
% % for in=1:Nz
% %     q(:,:,in)=real(ifft2(qh(:,:,in)));
% %     psi(:,:,in)=real(ifft2(psih(:,:,in)));
% %     JAC(:,:,in)=real(ifft2(Jh(:,:,in)));
% % end
% % 
% % 
% 
% 
% % figure(1);
% % pcolor2(JAC(:,:,1));colorbar;
% 
% 
% save('/Users/navid/Desktop/jacobmatlab.mat','jacobmatlab');
% 

psiih = tprod(invS,[1 2 3 -1],qih,[1 2 -1]);

ui = real(ifft2(-1i*KY.*psiih));
vi = real(ifft2(+1i*KX.*psiih));

for it=1:Nt

    [qhnew ,Ubnew ] = time_step_ETDRK4(qh,Ub,eL,eL2,Q,fu,fab,fc,dt);
%     [qhnew,Ubnew] = time_step_RK4(qh,Ub,dt);
%     [qhnew,Ubnew] = time_step_Euler(qh,Ub,dt);
    Ub = Ubnew;
    
    
%     qhnew(1:4,1:4,1)
%     qhnew2(1:4,1:4,1)
%     qhnew3(1:4,1:4,1)
    
%     qh = qhnew.*FILTERsmooth;
    qh = qhnew;

    psih = tprod(invS,[1 2 3 -1],qh,[1 2 -1]);
    psih_N = psih(:,:,Nz);
    fs = real(mean(conj(psih_N(:)).*(1i*KX(:).*etah(:)))) / (Nx*Ny);
    
    
    [ke,pe,keUb]=calc_kepe(psih,Ub);
     KE(:,it)=ke;
    KEUb(:,it)=keUb;
    if Nz>1
        PE(:,it)=pe;
    end
    formstress(it)=fs;
    Ubt(it)=Ub;
    
    
    if rem(it,5000)==1
        num2str(it)
        plot_figs_ETDRK4
    end
end

% q1h = qh;q2h = qih + (-fft2(jacobmatlab))*dt;norm(q1h(:)-q2h(:))/norm(q1h(:))


plot_figs_ETDRK4
asdfas
%%
display(['did ' num2str(it) ' time steps']);

qmatlab=real(ifft2(qh));
psih = tprod(invS,[1 2 3 -1],qh,[1 2 -1]);

umatlab= real(ifft2(-1i*KY.*psih));
vmatlab= real(ifft2(+1i*KX.*psih));
load('/Users/navid/Desktop/qi.mat');
qipython=qi;


load('/Users/navid/Desktop/qpython.mat');


load('/Users/navid/Desktop/psii.mat');
psiipython=psii;
qih=fft2(qimatlab);

psiih = tprod(invS,[1 2 3 -1],qih,[1 2 -1]);
psiimatlab=real(ifft2(psiih));


% 
% for in=1:3
%     figure(334)
%     subplot(3,3,in);
%     pcolor(X,Y,squeeze(psiipython(in,:,:)));shading flat;
%     colorbar;
%     subplot(3,3,in+3);
%     pcolor(X,Y,psiimatlab(:,:,in));shading flat;
%     colorbar;
%     subplot(3,3,in+6);
%     pcolor(X,Y,abs(squeeze(psiipython(in,:,:))-psiimatlab(:,:,in)));shading flat;
%     colorbar;
% end

load('/Users/navid/Desktop/jacobpython.mat');
load('/Users/navid/Desktop/jacobmatlab.mat');

load('/Users/navid/Desktop/upython.mat');
load('/Users/navid/Desktop/vpython.mat');


load('/Users/navid/Desktop/psipython.mat');

fp = psiipython;
fp = jacobpython;
fp = qipython - jacobpython*dt;
fp = qpython;
fp = permute(fp,[2 3 1]);

fm = jacobmatlab;
fm = qmatlab;
% fm = qimatlab - jacobmatlab*dt;
% fm = qmatlab;
% fm = jacobmatlab + beta*vmatlab;

% fp = qimatlab;
maxF = max(fp(:));
% fp=fp/maxF;
% fm=fm/maxF;

for in=1:3
    figure(334)
    subplot(3,3,in);
    pcolor(X,Y,fp(:,:,in));shading flat;
    if in==1,ylabel('python');end
    colorbar;
    subplot(3,3,in+3);
    pcolor(X,Y,fm(:,:,in));shading flat;
    if in==1,ylabel('matlab');end
    colorbar;
    subplot(3,3,in+6);
    pcolor(X,Y,(fp(:,:,in)-fm(:,:,in)));shading flat;
    colorbar;
    xp=fp(:,:,in);
    xm=fm(:,:,in);
    title(num2str(norm(xp(:)-xm(:))/norm(xm(:))))
end
asdfasdf

%%
qtopomatlab = qimatlab;
qtopomatlab(:,:,Nz)=qtopomatlab(:,:,Nz) + eta;

load('/Users/navid/Desktop/qtopopython.mat');
qtopopython=permute(qtopopython,[2 3 1]);
norm(qtopomatlab(:)-qtopopython(:))/norm(qtopomatlab(:))
%%

utotmatlab = ui;
for in=1:Nz
    utotmatlab(:,:,in) = utotmatlab(:,:,in) + Uj(in);
end

load('/Users/navid/Desktop/utot.mat');
utotpython=utot;
utotpython=permute(utotpython,[2 3 1]);
norm(utotmatlab(:)-utotpython(:))/norm(utotmatlab(:))
figure(2);
for in=1:Nz,
    subplot(3,3,in);pcolor(utotmatlab(:,:,in));colorbar;
    subplot(3,3,in+3);pcolor(utotpython(:,:,in));colorbar;
    subplot(3,3,in+6);pcolor(utotpython(:,:,in)-utotmatlab(:,:,in));colorbar;
    xm=utotmatlab(:,:,in);xp=utotpython(:,:,in);
    title(num2str(norm(xm(:)-xp(:))/norm(xm(:))))
end;


utoth_m = fft2(utotmatlab);utoth_p = fft2(utotpython);

%%

vmatlab = vi;

load('/Users/navid/Desktop/vpython.mat');
vpython=permute(vpython,[2 3 1]);
norm(vmatlab(:)-vpython(:))/norm(vmatlab(:))
figure(2);
for in=1:Nz
    subplot(3,3,in);pcolor(vmatlab(:,:,in));colorbar;
    subplot(3,3,in+3);pcolor(vpython(:,:,in));colorbar;
    subplot(3,3,in+6);pcolor(vpython(:,:,in)-vmatlab(:,:,in));colorbar;
    xm=vmatlab(:,:,in);xp=vpython(:,:,in);
    title(num2str(norm(xm(:)-xp(:))/norm(xm(:))))
end;


utoth_m = fft2(utotmatlab);utoth_p = fft2(utotpython);







%%

load ~/Desktop/tests/kqsmatlab.mat;
load ~/Desktop/tests/k1qpython.mat
load ~/Desktop/tests/k2qpython.mat
load ~/Desktop/tests/k3qpython.mat
load ~/Desktop/tests/k4qpython.mat


load ~/Desktop/tests/eL2python.mat
load ~/Desktop/tests/eLpython.mat
load ~/Desktop/tests/fabpython.mat
load ~/Desktop/tests/fupython.mat
load ~/Desktop/tests/fcpython.mat
load ~/Desktop/tests/Qpython.mat
load ~/Desktop/tests/Fnbpython.mat

load ~/Desktop/tests/Fnb_t.mat;
Fnb_p = Fnb_t;

[Fnb_m,nlin_U] = Dt_qUb_ETDRK4(fft2(k2qmatlab),Ub);

k2q = fft2(k2qpython);
k2q = permute(k2q,[2 3 1]);
[Fnb_p2,nlin_U] = Dt_qUb_ETDRK4(k2q,Ub);


% fp = Fnbpython;
% 
% 
% fp = Fnb_p;
% 
% fp = permute(fp,[2 3 1]);
% fp = real(ifft2(Fnb_m));

fp = k3qpython;
fp = permute(fp,[2 3 1]);


% fm = repmat(real(ifft2(eL2)),[1 1 3]);


fm = k3qmatlab;


for in=1:3
    figure(3434)
    subplot(3,3,in);
    pcolor(X,Y,fp(:,:,in));shading flat;
    if in==1,ylabel('python');end
    colorbar;
    subplot(3,3,in+3);
    pcolor(X,Y,fm(:,:,in));shading flat;
    if in==1,ylabel('matlab');end
    colorbar;
    subplot(3,3,in+6);
    pcolor(X,Y,(fp(:,:,in)-fm(:,:,in)));shading flat;
    colorbar;
    xp=fp(:,:,in);
    xm=fm(:,:,in);
    title(num2str(norm(xp(:)-xm(:))/norm(xm(:))))
end
