clearvars -global;
clearvars;

global KX KY Uj nu4 ksq rek etah Hj H Nx Ny Nz topoflag greduced f0 taub forceqh rhoj;

figNo=100;
cfllimit=0.8;

fac=1;
nofilter=1;
Lx=2*pi;Ly=2*pi;

Lx=1e6; %m
Ly=1e6; %mah
Ld = 15000; %m


topoflag=1;


nu=0e9; h=2;

rek = 1e-7;


Nx=32*fac; Ny=Nx; Nz=3;
dx=Lx/Nx;dy=Ly/Ny;
x=0:dx:Lx-dx;x=x+dx/2;
y=0:dy:Ly-dy;y=y+dy/2;
[X,Y]=meshgrid(x,y);


Hj=[500;1750;1750]; %m
H = sum(Hj(:));
rhoj=[1025;1025.275;1025.640]; %kg*m^(-3)
% rhoj=[1025;1028;1035]; %kg*m^(-3)
rho = mean(rhoj(:));
Uj = [0.05;0.025;0]*0;

% Hj = [500:500:5000]';
% % rhoj = [1025:-.5:1020.5];
% rhoj = [1025:-.5:1020.5]';
% Uj = [.05:-.005:.005]';

beta=1.2130692965249345e-11; % m^(-1)s^(-1)
f0 = 0.0001236812857687059; % s^(-1)
g = 9.81; % m*s^(-2)

hrms=200; % m
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

betaj = [beta                       - Fp(1)*(Uj(2)-Uj(1));
         beta - Fm(1)*(Uj(1)-Uj(2)) - Fp(2)*(Uj(3)-Uj(2));
         beta - Fm(2)*(Uj(2)-Uj(3))                       ];

     




kx=2*pi/Lx*[0:Nx/2-1 -Nx/2:-1];
ky=2*pi/Ly*[0:Ny/2-1 -Ny/2:-1];

[KX,KY]=meshgrid(kx,ky);
ksq=KX.^2+KY.^2;


eta = 2*etarms*cos(5*2*pi/Lx*X).*cos(5*2*pi/Ly*Y);
eta = sqrt(2)*etarms*cos(8*2*pi/Lx*X);
etah = fft2(eta);
etah = (randn(Ny,Nx)+1i*randn(Ny,Nx)).*(ksq>(2*pi/Lx*5)^2&ksq<(2*pi/Lx*7)^2);
eta = real(ifft2(etah));
eta = 2*etarms*cos(5*2*pi/Lx*X).*cos(5*2*pi/Ly*Y);
eta = etarms*eta/sqrt(mean(eta(:).^2));
etah=fft2(eta);


figure(1);clf;
pcolor2(X,Y,eta)



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

eta=randn(Ny,Nx,Nz)+1i*randn(Ny,Nx,Nz);
eta(1,1,:)=0*eta(1,1,:);
q = real(ifft2(1e-7*eta.*FILTERsmooth));
q = q-mean(q(:));
q=0*q;
 
% for in=1:Nz
%     q(:,:,in) = 1e-7*( cos(2*pi/Lx*3*X).*cos(2*pi/Ly*2*Y) - .4*cos(2*pi/Lx*3*X).*cos(2*pi/Ly*4*Y) );
% end
qimatlab=q;


Ub=0;

qh = fft2(q);
qh(:,:,Nz)=qh(:,:,Nz)+etah;


tic
psih=PVinverse(qh,invS);

toc
tic
psih2 = tprod(invS,[1 2 3 -1],qh,[1 2 -1]);
toc


% Jpsiqh = Jac_psiq(psih,qh);
% 
% Jpsiq = real(ifft2(Jpsiqh));
% 
% figure(32)
% for in=1:Nz
%     subplot(1,Nz,in);
%     pcolor2(X,Y,Jpsiq(:,:,in));shading interp;
% end

dt=300;

 tau = .25*sin(2*pi/Ly*(Y-Ly/2)); % N*m^(-2)
taub = .05; % N*m^(-2)

forceqh = -1i*KY.*fft2(tau) / (rhoj(1)*Hj(1));


Nt=10000;

KE   = zeros(Nz,Nt);
KEUb = zeros(1,Nt);
PE   = zeros(Nz-1,Nt);

Ubt = zeros(1,Nt);
formstress = zeros(1,Nt);

qht = qh;
if topoflag==1
    qht(:,:,Nz) = qh(:,:,Nz) - etah;
end
psih = tprod(invS,[1 2 3 -1],qht,[1 2 -1]);

for it=1:Nt

    
    [k1q,k1U,psih,fs] = Dt_qUb(qh         ,Ub         ,invS,betaj);
    [k2q,k2U,psih,fs] = Dt_qUb(qh+k1q*dt/2,Ub+k1U*dt/2,invS,betaj);
    [k3q,k3U,psih,fs] = Dt_qUb(qh+k2q*dt/2,Ub+k2U*dt/2,invS,betaj);
    [k4q,k4U,psih,fs] = Dt_qUb(qh+k3q*dt  ,Ub+k3U*dt  ,invS,betaj);
    
%     
    qhnew = qh + (k1q+2*k2q+2*k3q+k4q)*dt/6;
    Ubnew = Ub + (k1U+2*k2U+2*k3U+k4U)*dt/6;
    qh = qhnew.*FILTERsmooth;
    Ub = Ubnew;
    
    qht = qh;
    if topoflag==1
        qht(:,:,Nz) = qht(:,:,Nz) - etah;
    end
    psih = tprod(invS,[1 2 3 -1],qht,[1 2 -1]);
    psih_N = psih(:,:,Nz);
    fs = real(mean(conj(psih_N(:)).*(1i*KX(:).*etah(:)))) / (Nx*Ny);

% %     qhnew = qh + k1*dt;
%     qh = qhnew.*FILTERsmooth;
%     k1 = Dt_q(qh        ,invS,ones(size(FILTERsmooth)),betaj);
%     qhnew = qh + k1*dt;

%     if it==1
%         [k1,psih] = Dt_q(qh        ,invS,ones(size(FILTERsmooth)),betaj);
%         qhnew = qh + k1*dt;
%         k1p = k1;
%     elseif it==2
%         [k1,psih]  = Dt_q(qh        ,invS,ones(size(FILTERsmooth)),betaj);
%         qhnew = qh + (3*k1 -k1p)*dt/2;
%         k1pp = k1p;
%         k1p  = k1;
%     elseif it>2
%         [k1,psih]  = Dt_q(qh        ,invS,ones(size(FILTERsmooth)),betaj);
%         qhnew = qh + (23*k1 - 16*k1p + 5*k1pp)*dt/12;
%         k1pp = k1p;
%         k1p  = k1;
%     end
%     qh = qhnew.*FILTERsmooth + Fh*dt;

    
    
%     if it==1
%         [k1q,k1U,psih,fs] = Dt_qUb(qh,Ub,invS,betaj);
%         qhnew = qh + k1q*dt;
%         Ubnew = Ub + k1U*dt;
%         k1qp = k1q;k1Up = k1U;
%     elseif it==2
%         [k1q,k1U,psih,fs] = Dt_qUb(qh,Ub,invS,betaj);
%         qhnew = qh + (3*k1q -k1qp)*dt/2;
%         Ubnew = Ub + (3*k1U -k1Up)*dt/2;
%         k1qpp = k1qp; k1Upp = k1Up;
%         k1qp  = k1q;  k1Up  = k1U;
%     elseif it>2
%         [k1q,k1U,psih,fs] = Dt_qUb(qh,Ub,invS,betaj);
%         qhnew = qh + (23*k1q - 16*k1qp + 5*k1qpp)*dt/12;
%         Ubnew = Ub + (23*k1U - 16*k1Up + 5*k1Upp)*dt/12;
%         k1qpp = k1qp; k1Upp = k1Up;
%         k1qp  = k1q;  k1Up  = k1U;
%     end
%     qh = qhnew.*FILTERsmooth;
%     Ub = Ub;% Ubnew;

    [ke,pe,keUb]=calc_kepe(psih,Ub);
     KE(:,it)=ke;
    KEUb(:,it)=keUb;
    PE(:,it)=pe;
    formstress(it)=fs;
    Ubt(it)=Ub;
    
    
    if rem(it,1000)==1
        num2str(it)
        plot_figs
    end
end

plot_figs
asdf

qmatlab=real(ifft2(qh));

%%
load('/Users/navid/Desktop/qi.mat');
qipython=qi;
%%
load('/Users/navid/Desktop/qpythonC.mat');


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

fp = qpython;
size(fp)
fp=permute(fp,[2 3 1]);
size(fp)
fm = qmatlab;
% fp = qimatlab;


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
    pcolor(X,Y,abs(fp(:,:,in)-fm(:,:,in)));shading flat;
    colorbar;
    xp=fp(:,:,in);
    xm=fm(:,:,in);
    title(num2str(norm(xp(:)-xm(:))/norm(xm(:))))
        
end
%%
asdfasdf

