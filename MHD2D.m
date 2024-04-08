%MHD2D.m we use LDF basis
%2dproblem,we can set uz and Bz to zero,so there are 6 variables.

global gamma point weight phi phix phiy hh m M phiL phiR phiU phiD psi psix psiy psiL psiR psiU psiD
gamma=5/3;
[point,weight]=fourpoint_Gauss();

%Legendre basis
%on [-1,1]*[-1,1]
%%P2
phi1=@(x,y) 1;
phi2=@(x,y) x;
phi3=@(x,y) y;

phi4=@(x,y) x.^2-1/3;
phi5=@(x,y) y.^2-1/3;
phi6=@(x,y) x.*y;
m=[4,4/3,4/3,16/45,16/45,4/9];

phix1=@(x,y) 0.*x;
phix2=@(x,y) 1;
phix3=@(x,y) 0.*x;

phix4=@(x,y) 2*x;
phix5=@(x,y) 0.*x;
phix6=@(x,y) y;

phiy1=@(x,y) 0.*x;
phiy2=@(x,y) 0.*x;
phiy3=@(x,y) 1;

phiy4=@(x,y) 0.*y;
phiy5=@(x,y) 2*y;
phiy6=@(x,y) x;

phi=zeros(4,4,6);
phix=zeros(4,4,6);
phiy=zeros(4,4,6);
for i=1:4
    for j=1:4
        phi(i,j,1)=phi1(point(i),point(j));
        phi(i,j,2)=phi2(point(i),point(j));
        phi(i,j,3)=phi3(point(i),point(j));
        phi(i,j,4)=phi4(point(i),point(j));
        phi(i,j,5)=phi5(point(i),point(j));
        phi(i,j,6)=phi6(point(i),point(j));

        phix(i,j,1)=phix1(point(i),point(j));
        phix(i,j,2)=phix2(point(i),point(j));
        phix(i,j,3)=phix3(point(i),point(j));
        phix(i,j,4)=phix4(point(i),point(j));
        phix(i,j,5)=phix5(point(i),point(j));
        phix(i,j,6)=phix6(point(i),point(j));

        phiy(i,j,1)=phiy1(point(i),point(j));
        phiy(i,j,2)=phiy2(point(i),point(j));
        phiy(i,j,3)=phiy3(point(i),point(j));
        phiy(i,j,4)=phiy4(point(i),point(j));
        phiy(i,j,5)=phiy5(point(i),point(j));
        phiy(i,j,6)=phiy6(point(i),point(j));
    end
end

%%value on edge
%left (-1,point(q))
phiL=zeros(4,6);
for q=1:4
    phiL(q,1)=phi1(-1,point(q));
    phiL(q,2)=phi2(-1,point(q));
    phiL(q,3)=phi3(-1,point(q));
    phiL(q,4)=phi4(-1,point(q));
    phiL(q,5)=phi5(-1,point(q));
    phiL(q,6)=phi6(-1,point(q));
end

%right (1,point(p))
phiR=zeros(4,6);
for q=1:4
    phiR(q,1)=phi1(1,point(q));
    phiR(q,2)=phi2(1,point(q));
    phiR(q,3)=phi3(1,point(q));
    phiR(q,4)=phi4(1,point(q));
    phiR(q,5)=phi5(1,point(q));
    phiR(q,6)=phi6(1,point(q));
end

phiU=zeros(4,6);
for q=1:4
    phiU(q,1)=phi1(point(q),1);
    phiU(q,2)=phi2(point(q),1);
    phiU(q,3)=phi3(point(q),1);
    phiU(q,4)=phi4(point(q),1);
    phiU(q,5)=phi5(point(q),1);
    phiU(q,6)=phi6(point(q),1);
end

phiD=zeros(4,6);
for q=1:4
    phiD(q,1)=phi1(point(q),-1);
    phiD(q,2)=phi2(point(q),-1);
    phiD(q,3)=phi3(point(q),-1);
    phiD(q,4)=phi4(point(q),-1);
    phiD(q,5)=phi5(point(q),-1);
    phiD(q,6)=phi6(point(q),-1);
end


%%LDF basis
%LDF basis on [-1,1]*[-1,1]
psi1=@(x,y) [1;0];
psi2=@(x,y) [0;1];

psi3=@(x,y) [y;0];
psi4=@(x,y) [0;x];
psi5=@(x,y) [x;-y];

psi6=@(x,y) [y.^2-1/3;0];
psi7=@(x,y) [0;x.^2-1/3];

psi8=@(x,y) [x.^2-1/3;-2*x.*y];
psi9=@(x,y) [-2*x.*y;y.^2-1/3];

M=[4,4,4/3,4/3,8/3,16/45,16/45,32/15,32/15];

psi=zeros(4,4,9,2);
for p=1:4
    for q=1:4
        temp=psi1(point(p),point(q));
        psi(p,q,1,1)=temp(1);
        psi(p,q,1,2)=temp(2);

        temp=psi2(point(p),point(q));
        psi(p,q,2,1)=temp(1);
        psi(p,q,2,2)=temp(2);

        temp=psi3(point(p),point(q));
        psi(p,q,3,1)=temp(1);
        psi(p,q,3,2)=temp(2);

        temp=psi4(point(p),point(q));
        psi(p,q,4,1)=temp(1);
        psi(p,q,4,2)=temp(2);

        temp=psi5(point(p),point(q));
        psi(p,q,5,1)=temp(1);
        psi(p,q,5,2)=temp(2);

        temp=psi6(point(p),point(q));
        psi(p,q,6,1)=temp(1);
        psi(p,q,6,2)=temp(2);

        temp=psi7(point(p),point(q));
        psi(p,q,7,1)=temp(1);
        psi(p,q,7,2)=temp(2);

        temp=psi8(point(p),point(q));
        psi(p,q,8,1)=temp(1);
        psi(p,q,8,2)=temp(2);

        temp=psi9(point(p),point(q));
        psi(p,q,9,1)=temp(1);
        psi(p,q,9,2)=temp(2);

    end
end

%psix
psix1=@(x,y) [0.*x;0.*x];
psix2=@(x,y) [0.*x;0.*x];

psix3=@(x,y) [0.*x;0];
psix4=@(x,y) [0;1];
psix5=@(x,y) [1;0];

psix6=@(x,y) [0.*x;0];
psix7=@(x,y) [0;2*x];

psix8=@(x,y) [2*x;-2*y];
psix9=@(x,y) [-2*y;0];
%psiy
psiy1=@(x,y) [0.*y;0];
psiy2=@(x,y) [0;0.*y];

psiy3=@(x,y) [1;0];
psiy4=@(x,y) [0;0.*x];
psiy5=@(x,y) [0;-1];

psiy6=@(x,y) [2*y;0];
psiy7=@(x,y) [0;0.*x];

psiy8=@(x,y) [0;-2*x];
psiy9=@(x,y) [-2*x;2*y];

psix=zeros(4,4,9,2);
psiy=zeros(4,4,9,2);
for p=1:4
    for q=1:4
        temp=psix1(point(p),point(q));
        psix(p,q,1,1)=temp(1);
        psix(p,q,1,2)=temp(2);

        temp=psix2(point(p),point(q));
        psix(p,q,2,1)=temp(1);
        psix(p,q,2,2)=temp(2);

        temp=psix3(point(p),point(q));
        psix(p,q,3,1)=temp(1);
        psix(p,q,3,2)=temp(2);

        temp=psix4(point(p),point(q));
        psix(p,q,4,1)=temp(1);
        psix(p,q,4,2)=temp(2);

        temp=psix5(point(p),point(q));
        psix(p,q,5,1)=temp(1);
        psix(p,q,5,2)=temp(2);

        temp=psix6(point(p),point(q));
        psix(p,q,6,1)=temp(1);
        psix(p,q,6,2)=temp(2);

        temp=psix7(point(p),point(q));
        psix(p,q,7,1)=temp(1);
        psix(p,q,7,2)=temp(2);

        temp=psix8(point(p),point(q));
        psix(p,q,8,1)=temp(1);
        psix(p,q,8,2)=temp(2);

        temp=psix9(point(p),point(q));
        psix(p,q,9,1)=temp(1);
        psix(p,q,9,2)=temp(2);


        temp=psiy1(point(p),point(q));
        psiy(p,q,1,1)=temp(1);
        psiy(p,q,1,2)=temp(2);

        temp=psiy2(point(p),point(q));
        psiy(p,q,2,1)=temp(1);
        psiy(p,q,2,2)=temp(2);

        temp=psiy3(point(p),point(q));
        psiy(p,q,3,1)=temp(1);
        psiy(p,q,3,2)=temp(2);

        temp=psiy4(point(p),point(q));
        psiy(p,q,4,1)=temp(1);
        psiy(p,q,4,2)=temp(2);

        temp=psiy5(point(p),point(q));
        psiy(p,q,5,1)=temp(1);
        psiy(p,q,5,2)=temp(2);

        temp=psiy6(point(p),point(q));
        psiy(p,q,6,1)=temp(1);
        psiy(p,q,6,2)=temp(2);

        temp=psiy7(point(p),point(q));
        psiy(p,q,7,1)=temp(1);
        psiy(p,q,7,2)=temp(2);

        temp=psiy8(point(p),point(q));
        psiy(p,q,8,1)=temp(1);
        psiy(p,q,8,2)=temp(2);

        temp=psiy9(point(p),point(q));
        psiy(p,q,9,1)=temp(1);
        psiy(p,q,9,2)=temp(2);
    end
end
% hh=1;


%value on edge of LDF basis
psiL=zeros(4,9,2);
psiR=zeros(4,9,2);
psiU=zeros(4,9,2);
psiD=zeros(4,9,2);
for p=1:4
    %psiL
    temp=psi1(-1,point(p));
    psiL(p,1,1)=temp(1);
    psiL(p,1,2)=temp(2);

    temp=psi2(-1,point(p));
    psiL(p,2,1)=temp(1);
    psiL(p,2,2)=temp(2);

    temp=psi3(-1,point(p));
    psiL(p,3,1)=temp(1);
    psiL(p,3,2)=temp(2);

    temp=psi4(-1,point(p));
    psiL(p,4,1)=temp(1);
    psiL(p,4,2)=temp(2);

    temp=psi5(-1,point(p));
    psiL(p,5,1)=temp(1);
    psiL(p,5,2)=temp(2);

    temp=psi6(-1,point(p));
    psiL(p,6,1)=temp(1);
    psiL(p,6,2)=temp(2);

    temp=psi7(-1,point(p));
    psiL(p,7,1)=temp(1);
    psiL(p,7,2)=temp(2);

    temp=psi8(-1,point(p));
    psiL(p,8,1)=temp(1);
    psiL(p,8,2)=temp(2);

    temp=psi9(-1,point(p));
    psiL(p,9,1)=temp(1);
    psiL(p,9,2)=temp(2);
    %psiR
    temp=psi1(1,point(p));
    psiR(p,1,1)=temp(1);
    psiR(p,1,2)=temp(2);

    temp=psi2(1,point(p));
    psiR(p,2,1)=temp(1);
    psiR(p,2,2)=temp(2);

    temp=psi3(1,point(p));
    psiR(p,3,1)=temp(1);
    psiR(p,3,2)=temp(2);

    temp=psi4(1,point(p));
    psiR(p,4,1)=temp(1);
    psiR(p,4,2)=temp(2);

    temp=psi5(1,point(p));
    psiR(p,5,1)=temp(1);
    psiR(p,5,2)=temp(2);

    temp=psi6(1,point(p));
    psiR(p,6,1)=temp(1);
    psiR(p,6,2)=temp(2);

    temp=psi7(1,point(p));
    psiR(p,7,1)=temp(1);
    psiR(p,7,2)=temp(2);

    temp=psi8(1,point(p));
    psiR(p,8,1)=temp(1);
    psiR(p,8,2)=temp(2);

    temp=psi9(1,point(p));
    psiR(p,9,1)=temp(1);
    psiR(p,9,2)=temp(2);
    %psiU
    temp=psi1(point(p),1);
    psiU(p,1,1)=temp(1);
    psiU(p,1,2)=temp(2);

    temp=psi2(point(p),1);
    psiU(p,2,1)=temp(1);
    psiU(p,2,2)=temp(2);

    temp=psi3(point(p),1);
    psiU(p,3,1)=temp(1);
    psiU(p,3,2)=temp(2);

    temp=psi4(point(p),1);
    psiU(p,4,1)=temp(1);
    psiU(p,4,2)=temp(2);

    temp=psi5(point(p),1);
    psiU(p,5,1)=temp(1);
    psiU(p,5,2)=temp(2);

    temp=psi6(point(p),1);
    psiU(p,6,1)=temp(1);
    psiU(p,6,2)=temp(2);

    temp=psi7(point(p),1);
    psiU(p,7,1)=temp(1);
    psiU(p,7,2)=temp(2);

    temp=psi8(point(p),1);
    psiU(p,8,1)=temp(1);
    psiU(p,8,2)=temp(2);

    temp=psi9(point(p),1);
    psiU(p,9,1)=temp(1);
    psiU(p,9,2)=temp(2);

    %psiD
    temp=psi1(point(p),-1);
    psiD(p,1,1)=temp(1);
    psiD(p,1,2)=temp(2);

    temp=psi2(point(p),-1);
    psiD(p,2,1)=temp(1);
    psiD(p,2,2)=temp(2);

    temp=psi3(point(p),-1);
    psiD(p,3,1)=temp(1);
    psiD(p,3,2)=temp(2);

    temp=psi4(point(p),-1);
    psiD(p,4,1)=temp(1);
    psiD(p,4,2)=temp(2);

    temp=psi5(point(p),-1);
    psiD(p,5,1)=temp(1);
    psiD(p,5,2)=temp(2);

    temp=psi6(point(p),-1);
    psiD(p,6,1)=temp(1);
    psiD(p,6,2)=temp(2);

    temp=psi7(point(p),-1);
    psiD(p,7,1)=temp(1);
    psiD(p,7,2)=temp(2);

    temp=psi8(point(p),-1);
    psiD(p,8,1)=temp(1);
    psiD(p,8,2)=temp(2);

    temp=psi9(point(p),-1);
    psiD(p,9,1)=temp(1);
    psiD(p,9,2)=temp(2);
end


%%
% initial condition
p = @(x,y)  1 - (x.^2 + y.^2)./(8*pi.^2)*exp(1 - x.^2 - y.^2);
v1 = @(x,y) 1 - 1./(2.*pi).*exp(0.5.*(1 - (x.^2 + y.^2))).*y;
v2 = @(x,y) 1 + 1./(2.*pi).*exp(0.5.*(1 - (x.^2 + y.^2))).*x;
%rho
f1 = @(x,y) 1 + 0.*x;
%rhoux
f2 = @(x,y) f1(x,y).*v1(x,y);
%rhouy
f3 = @(x,y) f1(x,y).*v2(x,y);

%magnetic field
f5 = @(x,y) - 1./(2.*pi).*exp(0.5.*(1 - (x.^2 + y.^2))).*y;
f6 = @(x,y) 1./(2.*pi).*exp(0.5.*(1 - (x.^2 + y.^2))).*x;

%energy
f4 = @(x,y) p(x,y)/(gamma-1) + 0.5*(f2(x,y).^2 + f3(x,y).^2)./f1(x,y) + 0.5*(f5(x,y).^2 + f6(x,y).^2);

xa = -10; xb = 10;
ya = -10; yb = 10;

%%
% % scalar smooth example
%  p = @(x,y) 5;
%  v1 = @(x,y) 1;
%  v2 = @(x,y) 1;
% 
% %rho
%  f1 = @(x,y) 2 + sin(x+y);
% %rhoux
%  f2 = @(x,y) f1(x,y).*v1(x,y);
% %rhouy
%  f3 = @(x,y) f1(x,y).*v2(x,y);
% 
% %magnetic field
%  f5 = @(x,y) 0;
%  f6 = @(x,y) 0;
% 
% %energy
%  f4 = @(x,y) p(x,y)/(gamma-1) + 0.5*(f2(x,y).^2 + f3(x,y).^2)./f1(x,y) + 0.5*(f5(x,y).^2 + f6(x,y).^2);
%  xa = 0; xb = 2*pi;
%  ya = 0; yb = 2*pi; 

%%
% Orszag-Tang Vortex
gamma=5/3;
p = @(x,y)  gamma;
v1 = @(x,y) -sin(y);
v2 = @(x,y) sin(x);
%rho
f1 = @(x,y) gamma^2;
%rhoux
f2 = @(x,y) f1(x,y).*v1(x,y);
%rhouy
f3 = @(x,y) f1(x,y).*v2(x,y);

%magnetic field
f5 = @(x,y) -sin(y);
f6 = @(x,y) sin(2*x);

%energy
f4 = @(x,y) p(x,y)/(gamma-1) + 0.5*(f2(x,y).^2 + f3(x,y).^2)./f1(x,y) + 0.5*(f5(x,y).^2 + f6(x,y).^2);

xa = 0; xb = 2*pi;
ya = 0; yb = 2*pi;

%%
tFinal = 10;
N = 64; 
CFL=0.18;
tic

%mesh
h = (xb - xa)/N;
hh = h/2;

phix=phix/hh;
phiy=phiy/hh;
psix=psix/hh;
psiy=psiy/hh;

xc=linspace(xa,xb,N+1);
yc=linspace(ya,yb,N+1);


xmid = (xc(1:end-1) + xc(2:end))/2;
ymid = (yc(1:end-1) + yc(2:end))/2;
[Xmid,Ymid] = meshgrid(xmid,ymid);


uinit=zeros(N,N,4,4,4);
uinitB=zeros(N,N,4,4,2);

%initial value on gauss point
for i=1:N
    for j=1:N
        pointx=hh*point+Xmid(i,j);
        pointy=hh*point+Ymid(i,j);
        for p=1:4
            for q=1:4
                uinit(i,j,p,q,1)=f1(pointx(p),pointy(q));
                uinit(i,j,p,q,2)=f2(pointx(p),pointy(q));
                uinit(i,j,p,q,3)=f3(pointx(p),pointy(q));
                uinit(i,j,p,q,4)=f4(pointx(p),pointy(q));

                uinitB(i,j,p,q,1)=f5(pointx(p),pointy(q));
                uinitB(i,j,p,q,2)=f6(pointx(p),pointy(q));
            end
        end
    end
end

%initial coefficient
uh=zeros(N,N,6,4);
uhB=zeros(N,N,9,1);
for i=1:N
    for j=1:N
        for k=1:6
            for p=1:4
                for q=1:4
                    uh(i,j,k,:)=uh(i,j,k,:)+weight(p)*weight(q)*reshape(uinit(i,j,p,q,:),[1,1,1,4]).*phi(p,q,k);
                end
            end
        end
    end
end
for k=1:6
    uh(:,:,k,:)=uh(:,:,k,:)/m(k);
end

for i=1:N
    for j=1:N
        for k=1:9
            for p=1:4
                for q=1:4
                    uhB(i,j,k,1)=uhB(i,j,k,1)+weight(p)*weight(q)*(uinitB(i,j,p,q,1)*psi(p,q,k,1)+uinitB(i,j,p,q,2)*psi(p,q,k,2));
                end
            end
        end
    end
end

for k=1:9
    uhB(:,:,k,:)=uhB(:,:,k,:)/M(k);
end


t=0;
%RK3
while t<tFinal
    uhG=ValueGausspoint(uh);
    uhGB=ValueGausspointB(uhB);

    [alphax,alphay]=Alpha(uhG,uhGB,N,4);

    dt=CFL/(alphax/h+alphay/h);

    if t+dt>=tFinal
        dt=tFinal-t;
        t=tFinal;
    else
        t=t+dt;
    end

    [du,duB]=Lh(uh,uhB);
    uh1=uh+dt*du;
    uhB1=uhB+dt*duB;

    [du,duB]=Lh(uh1,uhB1);
    uh2=3/4*uh+1/4*uh1+1/4*dt*du;
    uhB2=3/4*uhB+1/4*uhB1+1/4*dt*duB;

    [du,duB]=Lh(uh2,uhB2);
    uh=1/3*uh+2/3*uh2+2/3*dt*du;
    uhB=1/3*uhB+2/3*uhB2+2/3*dt*duB;
    hold on
    uref=zeros(N,N);

    for i=1:N
    for j=1:N
        uref(i,j)=uh(i,j,1,1)-1/3*uh(i,j,4,1)-1/3*uh(i,j,5,1);%-f4(Xmid(i,j)-tFinal,Ymid(i,j)-tFinal);
    end
    end
    %
    contour(Xmid,Ymid,uref);
    xlabel('x');
    pause(0.0001)
    clf
end
%refrence solution



%Bxerr
uerr=zeros(N,N,4,4);
error=0;
uhG=ValueGausspoint(uh);
uhGB=ValueGausspointB(uhB);
for i=1:N
    for j=1:N
        for p=1:4
            for q=1:4
                x0=hh*point(p)+Xmid(i,j);
                y0=hh*point(q)+Ymid(i,j);
                x1=x0-tFinal;
                y1=y0-tFinal;

                %uerr(i,j,p,q)=abs(f5(x1,y1)-uhGB(i,j,p,q,1));
                %uerr(i,j,p,q)=abs(f4(x1,y1)-uhG(i,j,p,q,1));
                %error=error+weight(p)*weight(q)*uerr(i,j,p,q)^2*hh^2;
            end
        end
    end
end
toc
disp(sqrt(error));
axis on

