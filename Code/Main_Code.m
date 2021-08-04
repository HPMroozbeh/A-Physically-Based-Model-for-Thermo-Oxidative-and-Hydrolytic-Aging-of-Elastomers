% A Physically-Based Model for Thermo-Oxidative and Hydrolytic Aging of Elastomers

% Nonlinear Finite element code is written by : AMIR BAHROLOLUMIE %
% bahrolol@msu.edu

%% Read Data
clc
clear all
close all
fid = fopen('coarse.msh','r');
fgetl(fid);
fgetl(fid);
coord = fscanf(fid,'%e',[4,inf])';
fgetl(fid);
fgetl(fid);
fgetl(fid);
connect = fscanf(fid,'%e',[5,inf])';
xx      = coord(:,2);
yy      = coord(:,3);
nn      = max(coord(:,1));
ne      = max(connect(:,1));

%% Constitutive Model Material Parameters
n0       = 2.497;
rho_0    =  1.4416;       % rho_0 denotes the initial average  number  of  polymer  chains  per  unit  of  volume
rho_inf  =  2.2072;
mu       = 0.52 ;
pho      =  4.098748165807843e+18;
E        =  1.614741099118090e+04;
%% define Aging parameter

t        = 30;    % the unit is day
T        = 353;   % the unit is K

% Please not that here here rho0 is for thermo-oxidative aging, if you are wiling to do it for hydrolytic aging you need to change it to rho0     = rho0 * exp(-pho.*t.*exp(-E/T));
% rho0     = rho_inf + rho_0*exp(-pho.*t.*exp(-E/T));
rho0     = rho_inf - rho_0*exp(-pho.*t.*exp(-E/T));

%% Number of Gauss points in each direction

NGP            = 2 ;
Gauss_number   = NGP*NGP;
[r, w]         = gaussian_points(NGP) ;
Ndir           = 8;                                     % Determine the number of gauss point( i.e. 8,42,45, and 400)
PP_Quantity    = zeros(Gauss_number*ne,1);



%% Implement Loading and B.C
ID      = fopen('tecplot_data.dat','wt');
step    = 50;
tlr     = 1e-7;
UX      = max(xx);
UY      = 0;
%%% Ddof & Adof
% Loading Boundary Conditions degree of freedom
LDOF       = find(coord(:,2)== max(xx))';
D_LDOF_H   = 2.*LDOF - 1;
D_LDOF_V   = 2.*LDOF ;
% Geometric Boundary Conditions degree of freedom
GDOF       = find(coord(:,2)== min(xx))';
D_GDOF_H   = 2.*GDOF - 1;
D_GDOF_V   = 2.*GDOF ;
Ddof       = sort([D_LDOF_H,D_LDOF_V,D_GDOF_H,D_GDOF_V]);
Adof       = setdiff((1:2*nn)',Ddof);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U  = zeros(2*nn,1);
dU = zeros(2*nn,1);
for i=1:step + 1
    counter = 0;
    converge  = 1 ;
    sai_gameghabl = zeros(2*nn,1);
    while converge >tlr
        counter  = counter+1;
        Kmat       = zeros(2*nn,2*nn);
        Kgeo       = zeros(2*nn,2*nn);
        KT         = zeros(2*nn,2*nn);
        dU         = zeros(2*nn,1);
        dF         = zeros(2*nn,1);
        F_internal = zeros(2*nn,1);
        
        for t=1:ne
            Gaus_counter     = 0;
            elementcourse    = connect(t,:);
            elementcourse(1) = [];
            ElementDof       = [[2*connect(t,2)-1] [2*connect(t,2)] [2*connect(t,3)-1] [2*connect(t,3)] [2*connect(t,4)-1] [2*connect(t,4)] [2*connect(t,5)-1]  [2*connect(t,5)]];
            Kmat_el          = zeros(8,8);
            Kgeo_el         = zeros(8,8);
            KT_el           = zeros(8,8);
            F_internal_el   = zeros(8,1);
            u_eleman        = U(ElementDof,1);
            for m=1:NGP
                x = r(m);
                for k=1:NGP
                    y = r(k);
                    Gaus_counter = Gaus_counter +1;
                    
                    
                    [N, B_bar, G,j,F,C]                                  = Geometric_Nonlinearity_2D(x,y,u_eleman,elementcourse,xx,yy) ;
                    [S,Ms,rho,Sigma]                                     = SP_Stress(rho0,n0,mu,Ndir,F,C);
                    PP_Quantity((Gauss_number*(t-1) + Gaus_counter),:)   = S(1);
                    
                    Cijkl                                                = tangantial_matrix(rho0,n0,mu,Ndir,F,C);
                    Kmat_el                                              = Kmat_el + B_bar'*Cijkl*B_bar*w(m)*w(k)*det(j) ;
                    Kgeo_el                                              = Kgeo_el + G'*Ms*G*w(m)*w(k)*det(j) ;
                    F_internal_el                                        = F_internal_el + B_bar'*S*w(m)*w(k)*det(j) ;
                end
            end
            
            KT_el                     = Kmat_el + Kgeo_el ;
            KT(ElementDof,ElementDof) = KT_el + KT(ElementDof,ElementDof);
            F_internal(ElementDof,1)  = F_internal_el + F_internal(ElementDof,1);
        end
        converge      = abs((norm(F_internal)-norm(sai_gameghabl))/(norm(sai_gameghabl)));
        sai_gameghabl = F_internal;
        
        if counter==1
            dU(Ddof)=0;
            du_Horizental = UX/step;
            du_Vertical   = UY/step;
            
            dU(D_LDOF_H,1)= du_Horizental;
            dU(D_LDOF_V,1)= du_Vertical;
            
            dU(Adof)=KT(Adof,Adof)\(-F_internal(Adof)-KT(Adof,Ddof)*dU(Ddof)) ;
            U=U+dU;
            
        else
            dU(Ddof)=0;
            dU(Adof)=KT(Adof,Adof)\(-F_internal(Adof)-KT(Adof,Ddof)*dU(Ddof)) ;
            
            U=U+dU;
        end
    end
    
    iterations(i) = counter;
    F_internal    = zeros(2*nn,1);
    for t=1:ne
        F_internal_el=zeros(8,1);
        elementcourse    = connect(t,:);
        elementcourse(1) = [];
        ElementDof       = [[2*connect(t,2)-1] [2*connect(t,2)] [2*connect(t,3)-1] [2*connect(t,3)] [2*connect(t,4)-1] [2*connect(t,4)] [2*connect(t,5)-1]  [2*connect(t,5)]];
        for m=1:NGP
            x=r(m);
            for k=1:NGP
                y=r(k);
                [N, B_bar, G,j,F,C]  = Geometric_Nonlinearity_2D(x,y,u_eleman,elementcourse,xx,yy) ;
                [S,Ms,rho,Sigma]     = SP_Stress(rho0,n0,mu,Ndir,F,C);
                F_internal_el        = F_internal_el + B_bar'*S*w(m)*w(k)*det(j) ;
            end
        end
        F_internal(ElementDof,1) = F_internal_el + F_internal(ElementDof,1);
    end
    
    
    
    for t=1:nn;
        Pho  = 0 ;
        m    = 0;
        frgp = find(connect(:,2)==t)';
        if frgp > 0;
            fr = size(frgp,2);
            m  = fr + m;
            for z = 1 : fr
                Pho = PP_Quantity(Gauss_number*(frgp(1,z)-1)+1,:)'+ Pho;
            end
        end
        
        scgp = find(connect(:,3)==t)';
        if scgp>0;
            sc = size(scgp,2);
            m      = sc + m;
            for z=1:sc
                Pho = PP_Quantity(Gauss_number*(scgp(1,z)-1)+2,:)'+ Pho;
            end
        end
        
        thrgp       = find(connect(:,4)==t)';
        if thrgp>0;
            thr = size(thrgp,2);
            m       = thr + m;
            for z=1:thr
                Pho = PP_Quantity(Gauss_number*(thrgp(1,z)-1)+3,:)'+ Pho;
            end
        end
        
        forgp = find(connect(:,5)==t)';
        if forgp>0;
            four = size(forgp,2);
            m         =   four + m;
            for z=1:four
                Pho = PP_Quantity(Gauss_number*(forgp(1,z)-1)+4,:)'+ Pho;
            end
        end
        
        if m==0;
            Pho     = 0;
        else
            
            Pho = (1/m)*Pho;
        end
        
        Density(t) = Pho;
    end
    
    
    
    
    
    
    
    
    
    Force(i) = 4*sum(F_internal(D_LDOF_H)) ;
    ndofn = 2*nn ;
    X     = coord(:,2)+U(1:2:ndofn) ;
    Y     = coord(:,3)+U(2:2:ndofn) ;
    elemf = connect ;
    fprintf (ID,'VARIABLES = "X","Y","S_{XX}","U(1:2:ndofn)","U(2:2:ndofn)" \n');
    fprintf(ID,['ZONE ,N =  ',num2str(size(X,1)),',E =  ',num2str(size(elemf,1)),',F = FEPOINT ,ET = QUADRILATERAL \n']);
    fprintf(ID,'%6.4f % 6.4f % 10.8f % 10.8f % 10.8f  \n',[X';Y';Density;U(1:2:ndofn)';U(2:2:ndofn)']);
    fprintf(ID,'%6.0f % 6.0f % 6.0f % 6.0f \n',[elemf(:,2)';elemf(:,3)';elemf(:,4)';elemf(:,5)']);
end



figure(2)
xx = coord(:,2) ;
yy =  coord(:,3) ;
ux =  U([1:2:size(U)],1) ;
[xi, yi] = meshgrid(...
    linspace(min(xx),max(xx),200),...
    linspace(min(yy),max(yy),200));
zi = griddata(xx,yy,ux, xi,yi);
contourf(xi,yi,zi,200, 'LineStyle','none')
daspect([1 1 1])
col = colorbar('vert');
title(col,'u_{x}(mm)','FontSize',10)
colormap(jet(50))

