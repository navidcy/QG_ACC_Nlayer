function  [qhnew,Ubnew] = time_step_ETDRK4(qh,Ub,eL,eL2,Q,fu,fab,fc,dt)

global FILTERsmooth;

[nlin0_q,nlin0_U] = Dt_qUb_ETDRK4(qh,Ub);
k1q = eL2.*qh + Q.*nlin0_q;   k1U = nlin0_U;
% k1q = k1q.*FILTERsmooth;

[nlin1_q,nlin1_U] = Dt_qUb_ETDRK4(k1q,Ub+k1U*dt/2);
k2q = eL2.*qh + Q.*nlin1_q;   k2U = nlin1_U;
% k2q = k2q.*FILTERsmooth;

[nlin2_q,nlin2_U] = Dt_qUb_ETDRK4(k2q,Ub+k2U*dt/2);
k3q = eL2.*k1q + Q.*(2*nlin2_q-nlin0_q);   k3U = nlin2_U;
% k3q = k3q.*FILTERsmooth;

[nlin3_q,nlin3_U] = Dt_qUb_ETDRK4(k3q,Ub+k3U*dt);
                                   k4U = nlin3_U;

qhnew  = eL.*qh + fu.*nlin0_q + 2*fab.*(nlin1_q+nlin2_q) +fc.*nlin3_q ;
Ubnew     = Ub + (k1U+2*k2U+2*k3U+k4U)*dt/6;
% k4qmatlab = real(ifft2(qhnew));

qhnew = qhnew.*FILTERsmooth;


% k1qmatlab = real(ifft2(k1q));
% k2qmatlab = real(ifft2(k2q));
% k3qmatlab = real(ifft2(k3q));
% save('~/Desktop/tests/kqsmatlab.mat','k1qmatlab','k2qmatlab','k3qmatlab','k4qmatlab','nlin2_q','k2q');
% display('saved');

end