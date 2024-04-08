%mhd2dnormaldg.m
%dont use LDF basis 

%2dproblem,we can set uz and Bz to zero,so there are 6 variables.
%clear,clc;
%%
global gamma point weight phi phix phiy hh m phiL phiR phiU phiD
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

%%
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
% f5=@(x,y) 0.*x;
% f6=@(x,y) 0.*x;
%energy
f4 = @(x,y) p(x,y)/(gamma-1) + 0.5*(f2(x,y).^2 + f3(x,y).^2)./f1(x,y) + 0.5*(f5(x,y).^2 + f6(x,y).^2);

xa = -10; xb = 10;
ya = -10; yb = 10; 

%%scalar smooth example
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



tFinal = 0.1;
N = 64;
CFL = 0.18;

%mesh

h = (xb - xa)/N;
hh=h/2;
xc=linspace(xa,xb,N+1);
yc=linspace(ya,yb,N+1);

phix=phix/hh;
phiy=phiy/hh;

xmid = (xc(1:end - 1) + xc(2:end))/2;
ymid = (yc(1:end - 1) + yc(2:end))/2;
[Xmid,Ymid] = meshgrid(xmid,ymid);

%%
uinit=zeros(N,N,4,4,6);

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

                uinit(i,j,p,q,5)=f5(pointx(p),pointy(q));
                uinit(i,j,p,q,6)=f6(pointx(p),pointy(q));
            end
        end
    end
end

%initial coefficient
uh=zeros(N,N,6,6);

for i=1:N
    for j=1:N
        for k=1:6
            for p=1:4
                for q=1:4
                    for dim=1:6
                    uh(i,j,k,dim)=uh(i,j,k,dim)+weight(p)*weight(q)*uinit(i,j,p,q,dim)*phi(p,q,k);
                    end
                end
            end
        end
    end
end
for k=1:6
    uh(:,:,k,:)=uh(:,:,k,:)/m(k);
end

dx=h;
t=0;
%RK3
tic
while t<tFinal
    uhG=ValueGausspoint(uh);
    [alphax,alphay]=Alphaa(uhG,N,4);
    %alphax=3;alphay=3;
    dt=CFL/(alphax/h+alphay/h);

    if t+dt>=tFinal
        dt=tFinal-t;
        t=tFinal;
    else
        t=t+dt;
    end

    du=Lh2(uh);
    uh1=uh+dt*du;

    du=Lh2(uh1);
    uh2=3/4*uh+1/4*uh1+1/4*dt*du;

    du=Lh2(uh2);
    uh=1/3*uh+2/3*uh2+2/3*dt*du;


end

%reference solution  energy
uref=zeros(N,N);
for i=1:N
    for j=1:N
        uref(i,j)=uh(i,j,1,4)-1/3*uh(i,j,4,4)-1/3*uh(i,j,5,4)-f4(Xmid(i,j)-tFinal,Ymid(i,j)-tFinal);
    end
end
uerr=zeros(N,N,4,4);
error=0;

uhG=ValueGausspoint(uh);
for i=1:N
    for j=1:N
        for p=1:4
            for q=1:4
                x0=hh*point(p)+Xmid(i,j);
                y0=hh*point(q)+Ymid(i,j);
                x1=x0-tFinal;
                y1=y0-tFinal;
         
                uerr(i,j,p,q)=abs(f5(x1,y1)-uhG(i,j,p,q,5));
                error=error+weight(p)*weight(q)*uerr(i,j,p,q)^2*hh^2;
            end
        end
    end
end
toc
disp(sqrt(error));
axis on

mesh(Xmid,Ymid,uref);
xlabel('x');
%err=abs(uref);
%disp(max(max(err)))