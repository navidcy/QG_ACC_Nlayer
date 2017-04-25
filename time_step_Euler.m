function  [qhnew,Ubnew] = time_step_Euler(qh,Ub,dt)

global nu ksq FILTERsmooth


[nlin0_q,nlin0_U] = Dt_qUb_ETDRK4(qh,Ub);
qhnew = qh + nlin0_q*dt - nu*ksq.^2.*qh;
qhnew = qhnew.*FILTERsmooth;
% psih = tprod(invS,[1 2 3 -1],qh,[1 2 -1]);
% JAC = Jac_psiqUbUj(psih,qh,0*Ub,0*Uj);
% qhnew = qh + (-JAC-beta*(1i*KX).*psih)*dt;

Ubnew     = Ub + nlin0_U*dt;

end