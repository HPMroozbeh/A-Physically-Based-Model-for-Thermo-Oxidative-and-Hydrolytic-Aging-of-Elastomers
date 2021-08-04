function Cijkl  = tangantial_matrix(rho0,n0,mu,Ndir,F,C);


[ D , W ]          = integration_points_data(Ndir) ;
invC               = inv(C);
Q                  = F(3,3)^2*[-F(2,2),F(2,1),0;F(1,2),-F(1,1),0;0,0,0];
dF33_dC            = 0.5*inv(F)*Q;
DD                 = zeros(3,3,3,3) ;

for ii = 1:size(D,1)
    d = D(ii,:)';
    YD_stretch                   = sqrt( d'*(F'*F)*d );  
    n                            =  n0 *exp(YD_stretch*mu) ;                               % n denotes    the  average  length  of  polymer  chains
    rho                          =  rho0*exp(-YD_stretch*mu) ;      
    
    Beta                         =  Inverse_langevin(YD_stretch/sqrt(n));
    dBeta                        =  dInverse_langevin(YD_stretch/sqrt(n));
    A                            =  d'*C*d;
    Sn                           =  W(ii)*(1/(YD_stretch)).*rho0*sqrt(n)*Beta;
    Cn                           =  W(ii)*rho*(dBeta/(1*YD_stretch^2) - sqrt(n)*Beta/(YD_stretch^3));
    
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    DD(i,j,k,l) = ((-2/3)*Sn*(d(i)*d(j)- F(3,3)^2*d(3)^2*invC(i,j))*invC(k,l)) ...
                        + (Cn*(d(i)*d(j) -invC(i,j)*(d(3)*F(3,3))^2)*(d(k)*d(l) - (A/3)*invC(k,l))) ...
                        + Sn*(-4*F(3,3)*d(3)^2*invC(i,j)*dF33_dC(k,l)  + 2*(F(3,3)*d(3))^2*invC(i,k)*invC(l,j))...
                        + DD(i,j,k,l) ;
                end
            end
        end
    end
    
end
CC = Voigt(DD);

CCC       =  2*F(3,3)*[dF33_dC(1,1)*CC(1,3),dF33_dC(2,2)*CC(1,3),2*dF33_dC(1,2)*CC(1,3);
    dF33_dC(1,1)*CC(2,3),dF33_dC(2,2)*CC(2,3),2*dF33_dC(1,2)*CC(2,3);
    dF33_dC(1,1)*CC(6,3),dF33_dC(2,2)*CC(6,3),2*dF33_dC(1,2)*CC(6,3)];

Cijkl     = CC([1,2,6],[1,2,6]) + CCC ;


end


