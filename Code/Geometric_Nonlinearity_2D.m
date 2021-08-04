function [N, B_bar, G,j,F,C]= Geometric_Nonlinearity_2D(x,y,u_eleman,elementcourse,xx,yy)

%% Shape Functions

N1 = 1/4*(1-x)*(1-y);
N2 = 1/4*(1+x)*(1-y);
N3 = 1/4*(1+x)*(1+y);
N4 = 1/4*(1-x)*(1+y);
N         =[N1, N2, N3, N4];
dN_dkesi  = [0.25*y-0.25 0.25-0.25*y 0.25+0.25*y -0.25-0.25*y];
dN_detta  = [0.25*x-0.25 -0.25-0.25*x 0.25+0.25*x 0.25-0.25*x];
%% jacobian matrix 
j         = zeros(2,2);
F         = zeros(2,2);
for q=1:4;
    j(1,1) = xx(elementcourse(q))*dN_dkesi(q)+j(1,1);
    j(1,2) = yy(elementcourse(q))*dN_dkesi(q)+j(1,2);
    j(2,1) = xx(elementcourse(q))*dN_detta(q)+j(2,1);
    j(2,2) = yy(elementcourse(q))*dN_detta(q)+j(2,2);
end
Nx       = zeros(1,4);
Ny       = zeros(1,4);
B_Linear = zeros(3,8);
%% B_linear matrix 
for q=1:4;
    ki        = zeros(3,2);
    assemble1 = [1 2 3];
    assemble2 = [2*q-1 2*q];
    ki        = (1/det(j))*[j(2,2)*dN_dkesi(q)-j(1,2)*dN_detta(q),0;
        0,-j(2,1)*dN_dkesi(q)+j(1,1)*dN_detta(q);
        -j(2,1)*dN_dkesi(q)+j(1,1)*dN_detta(q),j(2,2)*dN_dkesi(q)-j(1,2)*dN_detta(q)];
    Nx(1,q)   = (1/det(j))*(j(2,2)*dN_dkesi(q)-j(1,2)*dN_detta(q));
    Ny(1,q)   = (1/det(j))*(-j(2,1)*dN_dkesi(q)+j(1,1)*dN_detta(q));
    B_Linear(assemble1,assemble2) = ki;
end
G = zeros(4,8);
for q=1:4;
        ki = zeros(4,2);
        assemble1 = [1 2 3 4];
        assemble2 = [2*q-1 2*q];
        ki        = (1/det(j))*[j(2,2)*dN_dkesi(q)-j(1,2)*dN_detta(q),0;
            0,j(2,2)*dN_dkesi(q)-j(1,2)*dN_detta(q);
            -j(2,1)*dN_dkesi(q)+j(1,1)*dN_detta(q),0;
            0,-j(2,1)*dN_dkesi(q)+j(1,1)*dN_detta(q)];
        G(assemble1,assemble2) = ki;
end

ux=0 ; uy=0 ;vx=0 ; vy=0 ;
for q=1:4
    ux = ux +Nx(q)*u_eleman(2*q-1) ;
    uy = uy +Ny(q)*u_eleman(2*q-1) ;
    vx = vx +Nx(q)*u_eleman(2*q) ;
    vy = vy +Ny(q)*u_eleman(2*q) ;
end
EL    = [ux;vy;uy+vx] ;
%% Deformation Gradiant
F      = [ux,uy;vx,vy] + eye(2,2);
F(3,3) = 1/det(F);

A_teta = [ux vx 0 0;
    0 0 uy vy;
    uy vy ux vx];
teta   = [ux;vx;uy;vy];
ENL    = 0.5*A_teta*teta;
BNL    = A_teta*G;

E      = EL + ENL ;
B_bar  = B_Linear+BNL ;
%% Left cauchy tensor
C      = F'*F;

end

