function [S,Ms,rho,Sigma] = SP_Stress(rho0,n0,mu,Ndir,F,C)

[ D , W ]               = integration_points_data(Ndir) ;
inv_C                   = inv(C);

S     = 0;
Sigma = 0; 
for i = 1:size(D,1)
    %% Find the stretches in each direction of micro-sphere
    
    d = D(i,:)';
    YD_stretch                  =  sqrt( d'*(F'*F)*d );
    n                           =  n0 *exp(YD_stretch*mu) ;                               % n denotes    the  average  length  of  polymer  chains
    rho                         =  rho0*exp(-YD_stretch*mu) ;                            % rho denotes   the  average  number  of  polymer  chains  per  unit  of  volume
    Beta                        =  Inverse_langevin(YD_stretch/sqrt(n));
    Sn                          =  (1/(YD_stretch)).*rho*sqrt(n)*Beta;
    S                           =  W(i)*Sn*((d*d') - (F(3,3)*d(3))^2*inv_C) + S;
end

Sigma   = F*S*F'; 

S                  = [S(1,1);S(2,2);S(1,2)];
Ms     = [S(1) 0 S(3) 0;
    0 S(1) 0 S(3);
    S(3) 0 S(2) 0;
    0 S(3) 0 S(2)];

end

