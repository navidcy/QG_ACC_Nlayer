
plotU=1;
plotb=0;

np=round(Nx/24);

betaY=zeros(Ny,Nx,Nz);
   UY=zeros(Ny,Nx,Nz);
for in=1:Nz
    betaY(:,:,in) = betaj(in)*Y;
       UY(:,:,in) = -(Ub+Uj(in))*Y;
end


% psih = tprod(invS,[1 2 3 -1],qh,[1 2 -1]);

q = real(ifft2(qh)) +plotb*betaY;
psi=real(ifft2(psih)) - plotU*UY;

u = real(ifft2(-1i*KY.*psih));
v = real(ifft2(+1i*KX.*psih));



% cfl=max([max(max(abs(Ub+Uj(1)+u(:,:,1))));max(max(abs(Ub+Uj(2)+u(:,:,2))));max(max(abs(Ub+Uj(3)+u(:,:,3))))])*dt/dx;
cfl = calc_clf(Ub,Uj,u,v,dx,dt);
display(['cfl = ' num2str(cfl,'%1.1e')])
if cfl>cfllimit*10
    error('cfl>0.8')
end


% qf = [q q(:,1,:)];
% qf = [qf;qf(1,:,:)];
% 
% psif = [psi psi(:,1,:)];
% psif = [psif;psif(1,:,:)];
 
% uf = [u u(:,1,:)];
% uf = [uf;uf(1,:,:)];
% vf = [v v(:,1,:)];
% vf = [vf;vf(1,:,:)];
% 
% Xf = [X X(:,end)+dx];
% Xf = [Xf;Xf(end,:)];
% 
% Yf = [Y;Y(end,:)+dy];
% Yf = [Yf Yf(:,end)];

if Uj(1)==0
    Uscale=1;
else
    Uscale=Uj(1);
end

  
figure(figNo)
for in=1:Nz
    subplot(2,Nz,in);
    pcolor2(X/Ld,Y/Ld,q(:,:,in)/(Uscale/Ld));shading interp;
%     contourf(Xf,Yf,qf(:,:,in),21,'linestyle','none');
    colorbar;
%     colormap(gca,flipud(cbrewer('div','Spectral',64)));
    axis square;
    axis([0 Lx 0 Ly]/Ld)
    if plotb==1
        title(['$\beta_{' num2str(in) '}y + q_{' num2str(in) '}(x,y,t)$'],'interpreter','latex','fontsize',20)
    else
        title(['$ q_{' num2str(in) '}(x,y,t)$'],'interpreter','latex','fontsize',20)
    end
    set(gca,'Layer','top');box on;
    
    subplot(2,Nz,Nz+in);
    pcolor2(X/Ld,Y/Ld,psi(:,:,in));shading interp;
%     contourf(Xf,Yf,psif(:,:,in),21,'linestyle','none');    
    hold on;
    quiver(X(1:np:end,1:np:end)/Ld,Y(1:np:end,1:np:end)/Ld,plotU*(Ub+Uj(in))+u(1:np:end,1:np:end,in),v(1:np:end,1:np:end,in),1.2,'k')
    hold off;
    colorbar;
%     colormap(gca,flipud(cbrewer('div','Spectral',64)));
    axis square;
    set(gca,'Layer','top');box on;
    axis([0 Lx 0 Ly]/Ld)
    if plotU==1
        title(['$-U_{' num2str(in) '}y + \psi_{' num2str(in) '}(x,y,t)$'],'interpreter','latex','fontsize',20)
    else
        title(['$\psi_{' num2str(in) '}(x,y,t)$'],'interpreter','latex','fontsize',20)
    end
end



drawnow;

year=365*60*60*24;

PEnums = [1:Nz-1]-.5;
PEstr = {};
for in=1:Nz-1
    str=['$' num2str(2*in+1) '\big/2$'];
    PEstr = [PEstr,str];
end
KEstr = {};
for in=1:Nz
    str=['$' num2str(in) '$'];
    KEstr = [KEstr,str];
end

figure(figNo+1)
plot((1:it)*dt/year,Hj(Nz)/H * formstress(1:it)/(taub/(rhoj(1)*Hj(1))),(1:it)*dt/year,rek*Ubt(1:it)/(taub/(rhoj(1)*Hj(1))),'linewidth',2);
hold on;plot([1 it]*dt/year,[1 1],'--k',[1 it]*dt/year,[0 0],'--k','linewidth',1);hold off;
h=legend('$$\frac{\langle\psi_3\eta_x\rangle}{\langle\tau\rangle/(\varrho_1 H_1)}$$','$$\frac{\mu\,U_b}{\langle\tau\rangle/(\varrho_1 H_1)}$$');
set(h,'Orientation','horizontal','interpreter','latex','fontsize',14);
xlabel(['$t$ [year]'],'interpreter','latex','fontsize',20)



figure(figNo+2)
if Nz>1
    subplot(311)
    plot((1:it)*dt/year,KEUb(1:it),'linewidth',2);
    % hold on;plot([1 it]*dt/year,[1 1]*.5*((taub/(rhoj(1)*Hj(1)))/rek)^2,'--k','linewidth',1);hold off;
    ylabel(['$\frac1{2}U_b^2$'],'interpreter','latex','fontsize',20)
    xlabel(['$t$ [year]'],'interpreter','latex','fontsize',20)
    subplot(312)
    plot((1:it)*dt/year,KE(:,1:it),'linewidth',2);
    ylabel(['$\frac1{2}\langle|\nabla\psi_{n}|^2\rangle$'],'interpreter','latex','fontsize',20)
    xlabel(['$t$ [year]'],'interpreter','latex','fontsize',20)
    h=legend(KEstr);legtit = get(h,'title');set(legtit,'string','layer');
    set(h,'Orientation','vertical','interpreter','latex','fontsize',14,'Location','northwest')
    subplot(313)
    plot((1:it)*dt/year,PE(:,1:it),'linewidth',2);
    h=legend(PEstr);legtit = get(h,'title');set(legtit,'string','interface');
    set(h,'Orientation','vertical','interpreter','latex','fontsize',14,'Location','northwest')
    ylabel(['$\frac{f_0^2}{g''_n}\langle(\psi_{n+1}-\psi_{n})^2\rangle$'],'interpreter','latex','fontsize',20)
    xlabel(['$t$ [year]'],'interpreter','latex','fontsize',20)
else
    subplot(211)
    plot((1:it)*dt/year,KEUb(1:it),'linewidth',2);
    % hold on;plot([1 it]*dt/year,[1 1]*.5*((taub/(rhoj(1)*Hj(1)))/rek)^2,'--k','linewidth',1);hold off;
    ylabel(['$\frac1{2}U_b^2$'],'interpreter','latex','fontsize',20)
    xlabel(['$t$ [year]'],'interpreter','latex','fontsize',20)
    subplot(212)
    plot((1:it)*dt/year,KE(:,1:it),'linewidth',2);
    ylabel(['$\frac1{2}\langle|\nabla\psi_{n}|^2\rangle$'],'interpreter','latex','fontsize',20)
    xlabel(['$t$ [year]'],'interpreter','latex','fontsize',20)
end

drawnow;
