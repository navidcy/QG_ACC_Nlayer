function cfl = calc_clf(Ub,Uj,u,v,dx,dt);

global Nz

Uflow = 0*u;Vflow = 0*u;
for in=1:Nz
    Uflow(:,:,in) = Ub + Uj(in) + u(:,:,in);
    Vflow(:,:,in) =               v(:,:,in);
end


% Umag = sqrt(Uflow.^2+Vflow.^2);
Umag = abs([Uflow(:);Vflow(:)]);
cfl = max(Umag(:))*dt/dx;