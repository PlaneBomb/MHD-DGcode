% FortranMHD2D.m

%设定初值
clear,clc

xa = -10; xb = 10; ya = -10; yb = 10; %求解域

tb = 0.00; % final time

N = 32; % 网格数量
k = 1; % 多项式最高次数
dim1 = (k + 1)*(k + 2)/2; % 基的数目
dim0 = 9;
dim2 = 4; % 分量的数目
kt = 3; % RK的阶数

CFL = 0.08;

quad = 4; % 积分点的数量

plot = 1; % 画图

QF = zeros(N,N,1,dim2);
QF0 = zeros(N,N,1,2);
T = 0;

% 边界条件
p = @(x,y)  1 - (x.^2 + y.^2)./(8*pi.^2)*exp(1 - x.^2 - y.^2);
v1 = @(x,y) 1 - 1./(2.*pi).*exp(0.5.*(1 - (x.^2 + y.^2))).*y;
v2 = @(x,y) 1 + 1./(2.*pi).*exp(0.5.*(1 - (x.^2 + y.^2))).*x;

f1 = @(x,y) 1 + 0.*x;
f2 = @(x,y) f1(x,y).*v1(x,y);
f3 = @(x,y) f1(x,y).*v2(x,y);
f5 = @(x,y) - 1./(2.*pi).*exp(0.5.*(1 - (x.^2 + y.^2))).*y;
f6 = @(x,y) 1./(2.*pi).*exp(0.5.*(1 - (x.^2 + y.^2))).*x;
f4 = @(x,y) 1.5*p(x,y) + 0.5*f1(x,y).*(v1(x,y).^2 + v2(x,y).^2) + 0.5*(f5(x,y).^2 + f6(x,y).^2);

% 网格
h = (xb - xa)/N;
h1 = h/2;

for i = 1:N + 1
    X(i) = xa + (i - 1)*h;
    Y(i) = ya + (i - 1)*h;
end

% 整格点

Xc = (X(1:end - 1) + X(2:end))/2;
Yc = (Y(1:end - 1) + Y(2:end))/2;
[xc,yc] = meshgrid(Xc,Yc);

% Legendre多项式
P0 = @(x) sqrt(1/h) + 0.*x;
DP0 = @(x) 0.*x;
P1 = @(x) sqrt(3/h)*(x/h1);
DP1 = @(x) sqrt(3/h)*(1/h1) + 0.*x;
P2 = @(x) sqrt(5/h)*(3*(x/h1)^2 - 1)/2;
DP2 = @(x) sqrt(5/h)*3*x/(h1^2);

% LDF多项式
PB01 = @(x,y) [1/h;0];
PB02 = @(x,y) [0;1/h];

PB11 = @(x,y) [sqrt(1/h).*P1(y);0];
PB12 = @(x,y) [0;sqrt(1/h).*P1(x)];
PB13 = @(x,y) [sqrt(1/(2*h)).*P1(x);-sqrt(1/(2*h)).*P1(y)];

PB21 = @(x,y) [sqrt(1/h).*P2(y);0];
PB22 = @(x,y) [0;sqrt(1/h).*P2(x)];
PB23 = @(x,y) sqrt(5/6)*[P1(x).*P1(y);-1/(sqrt(5*h))*P2(y)];
PB24 = @(x,y) sqrt(5/6)*[1/sqrt(5*h)*P2(x);-P1(x).*P1(y)];


% 积分点和权重
switch quad
    case 1
        A = [0      2.0000000000000000000000000];
        
    case 2
        A = [0.5773502691896257645091488     1.0000000000000000000000000
            -0.5773502691896257645091488     1.0000000000000000000000000];
        
    case 3
        A = [0     0.8888888888888888888888889
            0.7745966692414833770358531     0.5555555555555555555555556
            -0.7745966692414833770358531     0.5555555555555555555555556];
        
    case 4
        A = [0.3399810435848562648026658     0.6521451548625461426269361
            0.8611363115940525752239465     0.3478548451374538573730639
            -0.3399810435848562648026658     0.6521451548625461426269361
            -0.8611363115940525752239465     0.3478548451374538573730639];
        
    case 5
        A = [0                                 0.5688888888888888888888889
            0.5384693101056830910363144     0.4786286704993664680412915
            0.9061798459386639927976269     0.2369268850561890875142640
            -0.5384693101056830910363144     0.4786286704993664680412915
            -0.9061798459386639927976269     0.2369268850561890875142640];
        
    case 6
        A = [0.2386191860831969086305017     0.4679139345726910473898703
            0.6612093864662645136613996     0.3607615730481386075698335
            0.9324695142031520278123016     0.1713244923791703450402961
            -0.2386191860831969086305017     0.4679139345726910473898703
            -0.6612093864662645136613996     0.3607615730481386075698335
            -0.9324695142031520278123016     0.1713244923791703450402961];
        
    case 7
        A = [0                                 0.4179591836734693877551020
            0.4058451513773971669066064     0.3818300505051189449503698
            0.7415311855993944398638648     0.2797053914892766679014678
            0.9491079123427585245261897     0.1294849661688696932706114
            -0.4058451513773971669066064     0.3818300505051189449503698
            -0.7415311855993944398638648     0.2797053914892766679014678
            -0.9491079123427585245261897     0.1294849661688696932706114];
        
    case 8
        A = [0.1834346424956498049394761     0.3626837833783619829651504
            0.5255324099163289858177390     0.3137066458778872873379622
            0.7966664774136267395915539     0.2223810344533744705443560
            0.9602898564975362316835609     0.1012285362903762591525314
            -0.1834346424956498049394761     0.3626837833783619829651504
            -0.5255324099163289858177390     0.3137066458778872873379622
            -0.7966664774136267395915539     0.2223810344533744705443560
            -0.9602898564975362316835609     0.1012285362903762591525314];
        
    case 9
        A = [0                                 0.3302393550012597631645251
            0.3242534234038089290385380     0.3123470770400028400686304
            0.6133714327005903973087020     0.2606106964029354623187429
            0.8360311073266357942994298     0.1806481606948574040584720
            0.9681602395076260898355762     0.0812743883615744119718922
            -0.3242534234038089290385380     0.3123470770400028400686304
            -0.6133714327005903973087020     0.2606106964029354623187429
            -0.8360311073266357942994298     0.1806481606948574040584720
            -0.9681602395076260898355762     0.0812743883615744119718922];
        
    case 10
        A = [0.1488743389816312108848260     0.2955242247147528701738930
            0.4333953941292471907992659     0.2692667193099963550912269
            0.6794095682990244062343274     0.2190863625159820439955349
            0.8650633666889845107320967     0.1494513491505805931457763
            0.9739065285171717200779640     0.0666713443086881375935688
            -0.1488743389816312108848260     0.2955242247147528701738930
            -0.4333953941292471907992659     0.2692667193099963550912269
            -0.6794095682990244062343274     0.2190863625159820439955349
            -0.8650633666889845107320967     0.1494513491505805931457763
            -0.9739065285171717200779640     0.0666713443086881375935688];
end

lambda = A(:,1);
weight = A(:,2);

lambdai = h1*lambda;

% 将多项式的值储存起来

VP = zeros(quad,quad,dim1,dim2); % 多项式在积分点的值
VPB = zeros(quad,quad,dim0,2); % LDF多项式在积分点的值
DxVPB = zeros(quad,quad,dim0,2); % LDF多项式的导数
DyVPB = zeros(quad,quad,dim0,2); % LDF多项式的导数
VP00 = zeros(dim1,1,1,dim2); % 多项式在整格点的值
VPB0 = zeros(dim0,1,1,2); % LDF多项式在整格点的值
VPdx = zeros(quad,quad,dim1,dim2); % 对x的导数
VPdy = zeros(quad,quad,dim1,dim2); % 对y的导数
VPpd = zeros(4,quad,dim1,dim2); % 多项式在边缘上的值，第一个索引的顺序是右、左、上、下
VPBpd = zeros(4,quad,dim0,2); % LDF多项式在边缘上的值，第一个索引的顺序是右、左、上、下


for p = 1:dim2
    VP00(1,1,1,p) = P0(0)*P0(0);
    VP00(2,1,1,p) = P1(0)*P0(0);
    VP00(3,1,1,p) = P0(0)*P1(0);
    VP00(4,1,1,p) = P2(0)*P0(0);
    VP00(5,1,1,p) = P1(0)*P1(0);
    VP00(6,1,1,p) = P0(0)*P2(0);
end

for i = 1:quad
    for j = 1:quad
        for p = 1:dim2
            VP(i,j,1,p) = P0(lambdai(i))*P0(lambdai(j));
            VP(i,j,2,p) = P1(lambdai(i))*P0(lambdai(j));
            VP(i,j,3,p) = P0(lambdai(i))*P1(lambdai(j));
            VP(i,j,4,p) = P2(lambdai(i))*P0(lambdai(j));
            VP(i,j,5,p) = P1(lambdai(i))*P1(lambdai(j));
            VP(i,j,6,p) = P0(lambdai(i))*P2(lambdai(j));
            
            VPdx(i,j,1,p) = DP0(lambdai(i))*P0(lambdai(j));
            VPdx(i,j,2,p) = DP1(lambdai(i))*P0(lambdai(j));
            VPdx(i,j,3,p) = DP0(lambdai(i))*P1(lambdai(j));
            VPdx(i,j,4,p) = DP2(lambdai(i))*P0(lambdai(j));
            VPdx(i,j,5,p) = DP1(lambdai(i))*P1(lambdai(j));
            VPdx(i,j,6,p) = DP0(lambdai(i))*P2(lambdai(j));
            
            VPdy(i,j,1,p) = P0(lambdai(i))*DP0(lambdai(j));
            VPdy(i,j,2,p) = P1(lambdai(i))*DP0(lambdai(j));
            VPdy(i,j,3,p) = P0(lambdai(i))*DP1(lambdai(j));
            VPdy(i,j,4,p) = P2(lambdai(i))*DP0(lambdai(j));
            VPdy(i,j,5,p) = P1(lambdai(i))*DP1(lambdai(j));
            VPdy(i,j,6,p) = P0(lambdai(i))*DP2(lambdai(j));
        end
    end
end

for i = 1:quad
    for j = 1:quad
        VPB(i,j,1,1) = 1/h;
        VPB(i,j,1,2) = 0;
        VPB(i,j,2,1) = 0;
        VPB(i,j,2,2) = 1/h;
        VPB(i,j,3,1) = (1/h)^0.5*P1(lambdai(j));
        VPB(i,j,3,2) = 0;
        VPB(i,j,4,1) = 0;
        VPB(i,j,4,2) = (1/h)^0.5*P1(lambdai(i));
        VPB(i,j,5,1) = (1/(2*h))^0.5*P1(lambdai(i));
        VPB(i,j,5,2) = -(1/(2*h))^0.5*P1(lambdai(j));
        VPB(i,j,6,1) = (1/h)^0.5*P2(lambdai(j));
        VPB(i,j,6,2) = 0;
        VPB(i,j,7,1) = 0;
        VPB(i,j,7,2) = (1/h)^0.5*P2(lambdai(i));
        VPB(i,j,8,1) = (5.0/6.0)^0.5*P1(lambdai(i))*P1(lambdai(j));
        VPB(i,j,8,2) = -(1.0/(6.0*h))^0.5*P2(lambdai(j));
        VPB(i,j,9,1) = (1.0/(6.0*h))^0.5*P2(lambdai(i));
        VPB(i,j,9,2) = -(5.0/6.0)^0.5*P1(lambdai(i))*P1(lambdai(j));
        
        DxVPB(i,j,1,1) = 0;
        DxVPB(i,j,1,2) = 0;
        DxVPB(i,j,2,1) = 0;
        DxVPB(i,j,2,2) = 0;
        DxVPB(i,j,3,1) = 0;
        DxVPB(i,j,3,2) = 0;
        DxVPB(i,j,4,1) = 0;
        DxVPB(i,j,4,2) = (1/h)^0.5*DP1(lambdai(i));
        DxVPB(i,j,5,1) = (1/(2*h))^0.5*DP1(lambdai(i));
        DxVPB(i,j,5,2) = 0;
        DxVPB(i,j,6,1) = 0;
        DxVPB(i,j,6,2) = 0;
        DxVPB(i,j,7,1) = 0;
        DxVPB(i,j,7,2) = (1/h)^0.5*DP2(lambdai(i));
        DxVPB(i,j,8,1) = (5.0/6.0)^0.5*DP1(lambdai(i))*P1(lambdai(j));
        DxVPB(i,j,8,2) = 0;
        DxVPB(i,j,9,1) = (1/(6*h))^0.5*DP2(lambdai(i));
        DxVPB(i,j,9,2) = -(5.0/6.0)^0.5*DP1(lambdai(i))*P1(lambdai(j));
        
        DyVPB(i,j,1,1) = 0;
        DyVPB(i,j,1,2) = 0;
        DyVPB(i,j,2,1) = 0;
        DyVPB(i,j,2,2) = 0;
        DyVPB(i,j,3,1) = (1/h)^0.5*DP1(lambdai(j));
        DyVPB(i,j,3,2) = 0;
        DyVPB(i,j,4,1) = 0;
        DyVPB(i,j,4,2) = 0;
        DyVPB(i,j,5,1) = 0;
        DyVPB(i,j,5,2) = -(1/(2*h))^0.5*DP1(lambdai(j));
        DyVPB(i,j,6,1) = (1/h)^0.5*DP2(lambdai(j));
        DyVPB(i,j,6,2) = 0;
        DyVPB(i,j,7,1) = 0;
        DyVPB(i,j,7,2) = 0;
        DyVPB(i,j,8,1) = (5.0/6.0)^0.5*P1(lambdai(i))*DP1(lambdai(j));
        DyVPB(i,j,8,2) = -(1/(6*h))^0.5*DP2(lambdai(j));
        DyVPB(i,j,9,1) = 0;
        DyVPB(i,j,9,2) = -(5.0/6.0)^0.5*P1(lambdai(i))*DP1(lambdai(j));
    end
end

VPB0(1,1,1,1) = 1/h;
VPB0(1,1,1,2) = 0;
VPB0(2,1,1,1) = 0;
VPB0(2,1,1,2) = 1/h;
VPB0(3,1,1,1) = (1/h)^0.5*P1(0);
VPB0(3,1,1,2) = 0;
VPB0(4,1,1,1) = 0;
VPB0(4,1,1,2) = (1/h)^0.5*P1(0);
VPB0(5,1,1,1) = (1/(2*h))^0.5*P1(0);
VPB0(5,1,1,2) = -(1/(2*h))^0.5*P1(0);
VPB0(6,1,1,1) = (1/h)^0.5*P2(0);
VPB0(6,1,1,2) = 0;
VPB0(7,1,1,1) = 0;
VPB0(7,1,1,2) = (1/h)^0.5*P2(0);
VPB0(8,1,1,1) = (5.0/6.0)^0.5*P1(0)*P1(0);
VPB0(8,1,1,2) = -(1/(6*h))^0.5*P2(0);
VPB0(9,1,1,1) = (1/(6*h))^0.5*P2(0);
VPB0(9,1,1,2) = -(5.0/6.0)^0.5*P1(0)*P1(0);

for i = 1:quad
    for p = 1:dim2
        VPpd(1,i,1,p) = P0(h1)*P0(lambdai(i));
        VPpd(1,i,2,p) = P1(h1)*P0(lambdai(i));
        VPpd(1,i,3,p) = P0(h1)*P1(lambdai(i));
        VPpd(1,i,4,p) = P2(h1)*P0(lambdai(i));
        VPpd(1,i,5,p) = P1(h1)*P1(lambdai(i));
        VPpd(1,i,6,p) = P0(h1)*P2(lambdai(i));
        
        VPpd(2,i,1,p) = P0(-h1)*P0(lambdai(i));
        VPpd(2,i,2,p) = P1(-h1)*P0(lambdai(i));
        VPpd(2,i,3,p) = P0(-h1)*P1(lambdai(i));
        VPpd(2,i,4,p) = P2(-h1)*P0(lambdai(i));
        VPpd(2,i,5,p) = P1(-h1)*P1(lambdai(i));
        VPpd(2,i,6,p) = P0(-h1)*P2(lambdai(i));
        
        VPpd(3,i,1,p) = P0(lambdai(i))*P0(h1);
        VPpd(3,i,2,p) = P1(lambdai(i))*P0(h1);
        VPpd(3,i,3,p) = P0(lambdai(i))*P1(h1);
        VPpd(3,i,4,p) = P2(lambdai(i))*P0(h1);
        VPpd(3,i,5,p) = P1(lambdai(i))*P1(h1);
        VPpd(3,i,6,p) = P0(lambdai(i))*P2(h1);
        
        VPpd(4,i,1,p) = P0(lambdai(i))*P0(-h1);
        VPpd(4,i,2,p) = P1(lambdai(i))*P0(-h1);
        VPpd(4,i,3,p) = P0(lambdai(i))*P1(-h1);
        VPpd(4,i,4,p) = P2(lambdai(i))*P0(-h1);
        VPpd(4,i,5,p) = P1(lambdai(i))*P1(-h1);
        VPpd(4,i,6,p) = P0(lambdai(i))*P2(-h1);
    end
end

for j = 1:quad
    VPBpd(1,j,1,1) = 1/h;
    VPBpd(1,j,1,2) = 0;
    VPBpd(1,j,2,1) = 0;
    VPBpd(1,j,2,2) = 1/h;
    VPBpd(1,j,3,1) = (1/h)^0.5*P1(lambdai(j));
    VPBpd(1,j,3,2) = 0;
    VPBpd(1,j,4,1) = 0;
    VPBpd(1,j,4,2) = (1/h)^0.5*P1(h1);
    VPBpd(1,j,5,1) = (1/(2*h))^0.5*P1(h1);
    VPBpd(1,j,5,2) = -(1/(2*h))^0.5*P1(lambdai(j));
    VPBpd(1,j,6,1) = (1/h)^0.5*P2(lambdai(j));
    VPBpd(1,j,6,2) = 0;
    VPBpd(1,j,7,1) = 0;
    VPBpd(1,j,7,2) = (1/h)^0.5*P2(h1);
    VPBpd(1,j,8,1) = (5.0/6.0)^0.5*P1(h1)*P1(lambdai(j));
    VPBpd(1,j,8,2) = -(1/(6*h))^0.5*P2(lambdai(j));
    VPBpd(1,j,9,1) = (1/(6*h))^0.5*P2(h1);
    VPBpd(1,j,9,2) = -(5.0/6.0)^0.5*P1(h1)*P1(lambdai(j));
    
    VPBpd(2,j,1,1) = 1/h;
    VPBpd(2,j,1,2) = 0;
    VPBpd(2,j,2,1) = 0;
    VPBpd(2,j,2,2) = 1/h;
    VPBpd(2,j,3,1) = (1/h)^0.5*P1(lambdai(j));
    VPBpd(2,j,3,2) = 0;
    VPBpd(2,j,4,1) = 0;
    VPBpd(2,j,4,2) = (1/h)^0.5*P1(-h1);
    VPBpd(2,j,5,1) = (1/(2*h))^0.5*P1(-h1);
    VPBpd(2,j,5,2) = -(1/(2*h))^0.5*P1(lambdai(j));
    VPBpd(2,j,6,1) = (1/h)^0.5*P2(lambdai(j));
    VPBpd(2,j,6,2) = 0;
    VPBpd(2,j,7,1) = 0;
    VPBpd(2,j,7,2) = (1/h)^0.5*P2(-h1);
    VPBpd(2,j,8,1) = (5.0/6.0)^0.5*P1(-h1)*P1(lambdai(j));
    VPBpd(2,j,8,2) = -(1/(6*h))^0.5*P2(lambdai(j));
    VPBpd(2,j,9,1) = (1/(6*h))^0.5*P2(-h1);
    VPBpd(2,j,9,2) = -(5.0/6.0)^0.5*P1(-h1)*P1(lambdai(j));
    
    VPBpd(3,j,1,1) = 1/h;
    VPBpd(3,j,1,2) = 0;
    VPBpd(3,j,2,1) = 0;
    VPBpd(3,j,2,2) = 1/h;
    VPBpd(3,j,3,1) = (1/h)^0.5*P1(h1);
    VPBpd(3,j,3,2) = 0;
    VPBpd(3,j,4,1) = 0;
    VPBpd(3,j,4,2) = (1/h)^0.5*P1(lambdai(j));
    VPBpd(3,j,5,1) = (1/(2*h))^0.5*P1(lambdai(j));
    VPBpd(3,j,5,2) = -(1/(2*h))^0.5*P1(h1);
    VPBpd(3,j,6,1) = (1/h)^0.5*P2(h1);
    VPBpd(3,j,6,2) = 0;
    VPBpd(3,j,7,1) = 0;
    VPBpd(3,j,7,2) = (1/h)^0.5*P2(lambdai(j));
    VPBpd(3,j,8,1) = (5.0/6.0)^0.5*P1(lambdai(j))*P1(h1);
    VPBpd(3,j,8,2) = -(1/(6*h))^0.5*P2(h1);
    VPBpd(3,j,9,1) = (1/(6*h))^0.5*P2(lambdai(j));
    VPBpd(3,j,9,2) = -(5.0/6.0)^0.5*P1(lambdai(j))*P1(h1);
    
    VPBpd(4,j,1,1) = 1/h;
    VPBpd(4,j,1,2) = 0;
    VPBpd(4,j,2,1) = 0;
    VPBpd(4,j,2,2) = 1/h;
    VPBpd(4,j,3,1) = (1/h)^0.5*P1(-h1);
    VPBpd(4,j,3,2) = 0;
    VPBpd(4,j,4,1) = 0;
    VPBpd(4,j,4,2) = (1/h)^0.5*P1(lambdai(j));
    VPBpd(4,j,5,1) = (1/(2*h))^0.5*P1(lambdai(j));
    VPBpd(4,j,5,2) = -(1/(2*h))^0.5*P1(-h1);
    VPBpd(4,j,6,1) = (1/h)^0.5*P2(-h1);
    VPBpd(4,j,6,2) = 0;
    VPBpd(4,j,7,1) = 0;
    VPBpd(4,j,7,2) = (1/h)^0.5*P2(lambdai(j));
    VPBpd(4,j,8,1) = (5.0/6.0)^0.5*P1(lambdai(j))*P1(-h1);
    VPBpd(4,j,8,2) = -(1/(6*h))^0.5*P2(-h1);
    VPBpd(4,j,9,1) = (1/(6*h))^0.5*P2(lambdai(j));
    VPBpd(4,j,9,2) = -(5.0/6.0)^0.5*P1(lambdai(j))*P1(-h1);
end

C = zeros(N + 2,N + 2,dim1,dim2);
C0 = zeros(N + 2,N + 2,dim0);
ff = zeros(quad,quad,1,dim2);
ffB = zeros(quad,quad,1,2);

% 求初始时间层的值
for i = 1:N
    for j = 1:N
        lambdai = (h*lambda + X(i) + X(i + 1))/2;
        lambdaj = (h*lambda + Y(j) + Y(j + 1))/2;
        
        % 积分点处的函数值
        for i1 = 1:quad
            for j1 = 1:quad
                ff(i1,j1,1,1) = f1(lambdai(i1),lambdaj(j1));
                ff(i1,j1,1,2) = f2(lambdai(i1),lambdaj(j1));
                ff(i1,j1,1,3) = f3(lambdai(i1),lambdaj(j1));
                ff(i1,j1,1,4) = f4(lambdai(i1),lambdaj(j1));
                ffB(i1,j1,1,1) = f5(lambdai(i1),lambdaj(j1));
                ffB(i1,j1,1,2) = f6(lambdai(i1),lambdaj(j1));
            end
        end
        
        % Gauss积分
        for deg = 1:dim1
            for i1 = 1:quad
                for j1 = 1:quad
                    C(i,j,deg,:) = C(i,j,deg,:) + weight(i1).*weight(j1).*ff(i1,j1,1,:).*VP(i1,j1,deg,:);
                end
            end
            C(i,j,deg,:) = C(i,j,deg,:)*h1^2;
        end
        
        for deg = 1:dim0
            for i1 = 1:quad
                for j1 = 1:quad
                    C0(i,j,deg) = C0(i,j,deg) + weight(i1).*weight(j1).*(ffB(i1,j1,1,1).*VPB(i1,j1,deg,1) + ffB(i1,j1,1,2).*VPB(i1,j1,deg,2));
                end
            end
            C0(i,j,deg) = C0(i,j,deg)*h1^2;
        end
        
    end
end

Q = zeros(N,N,1,dim2);
Q0 = zeros(N,N,1,2);

for i = 1:N
    for j = 1:N
        for deg = 1:dim1
            Q(i,j,1,:) = Q(i,j,1,:) + C(i,j,deg,:).*VP00(deg,1,1,:);
        end
        for deg = 1:dim0
            Q0(i,j,1,:) = Q0(i,j,1,:) + C0(i,j,deg).*VPB0(deg,1,1,:);
        end
    end
end

Qr = zeros(N,N,1,dim2);
Qr0 = zeros(N,N,1,2);
Qr = zeros(N,N,1,dim2);
for i=1:N
    for j=1:N
        Qr(i,j,1,2)=f2(Xc(i),Yc(j));
    end
end
for i = 1:N
    for j = 1:N
        Qr(i,j,1,1) = f1(Xc(i),Yc(j));
        Qr(i,j,1,2) = f2(Xc(i),Yc(j));
        Qr(i,j,1,3) = f3(Xc(i),Yc(j));
        Qr(i,j,1,4) = f4(Xc(i),Yc(j));
        Qr0(i,j,1,1) = f5(Xc(i),Yc(j));
        Qr0(i,j,1,2) = f6(Xc(i),Yc(j));
    end
end

E = abs(Q - Qr);
E0 = abs(Q0 - Qr0);

Linfty = max(max(max(max(E))));
Linfty0 = max(max(max(max(E0))));
testquadgk = Linfty;
testquadgk0 = Linfty0;

QF(:,:,1,:) = Q;
QF0(:,:,1,:) = Q0;

C(N + 1,:,:,:) = C(1,:,:,:);
C(:,N + 1,:,:) = C(:,1,:,:);
C(N + 2,:,:,:) = C(2,:,:,:);
C(:,N + 2,:,:) = C(:,2,:,:);

C0(N + 1,:,:) = C0(1,:,:);
C0(:,N + 1,:) = C0(:,1,:);
C0(N + 2,:,:) = C0(2,:,:);
C0(:,N + 2,:) = C0(:,2,:);

fprintf('初始时间层的误差为%d\n',Linfty)
    %s = [Xc(1),Xc(end),Yc(1),Yc(end),-0.5,0.5];
    %axis(s);
    axis on
    
    mesh(xc,yc,Qr(:,:,1,2));
    xlabel('x')
    %pause(0.1)

tend = 0;
lambdai = h1*lambda;



% 用RK3推进时间层
while tend < tb
    
    alphax = 0;
    alphay = 0;
    for i = 1:N
        for j = 1:N
            if cxB(Q(i,j,1,:),Q0(i,j,1,:)) > alphax
                alphax = cxB(Q(i,j,1,:),Q0(i,j,1,:));
            end
            if cyB(Q(i,j,1,:),Q0(i,j,1,:)) > alphay
                alphay = cyB(Q(i,j,1,:),Q0(i,j,1,:));
            end
        end
    end
    
    t = CFL*(h/alphax + h/alphay)^((k + 1)/kt);
    
    tend = tend + t;
    T = [T;tend];
    
    if tend > tb
        t = tb + t - tend;
        tend = tb;
    end
    
    
    ts = tic;
    if kt == 1
        [dC,dC0] = EulerB(C,C0,h,VP,VPdx,VPdy,VPpd,VPB,DxVPB,DyVPB,VPBpd,lambdai,weight);
        C = C + t*dC;
        C0 = C0 + t*dC0;
    elseif kt == 3
        [dC,dC0] = EulerB(C,C0,h,VP,VPdx,VPdy,VPpd,VPB,DxVPB,DyVPB,VPBpd,lambdai,weight);
        K1 = C + t*dC;
        K10 = C0 + t*dC0;
        testK1 = sum(sum(sum(sum(abs(K1)))));
        [LK1,LK10] = EulerB(K1,K10,h,VP,VPdx,VPdy,VPpd,VPB,DxVPB,DyVPB,VPBpd,lambdai,weight);
        testLK1 = sum(sum(sum(sum(abs(LK1)))));
        K2 = (3/4)*C + (1/4)*K1 + (1/4)*t*LK1;
        K20 = (3/4)*C0 + (1/4)*K10 + (1/4)*t*LK10;
        testK2 = sum(sum(sum(sum(abs(K2)))));
        [LK2,LK20] = EulerB(K2,K20,h,VP,VPdx,VPdy,VPpd,VPB,DxVPB,DyVPB,VPBpd,lambdai,weight);
        testLK2 = sum(sum(sum(sum(abs(LK2)))));
        C = (1/3)*C + (2/3)*K2 + (2/3)*t*LK2;
        C0 = (1/3)*C0 + (2/3)*K20 + (2/3)*t*LK20;
        testC = sum(sum(sum(sum(abs(C)))));
    end
    t0 = toc(ts);
    
    Q = zeros(N,N,1,dim2);
    Q0 = zeros(N,N,1,2);
    for i = 1:N
        for j = 1:N
            for deg = 1:dim1
                Q(i,j,1,:) = Q(i,j,1,:) + C(i,j,deg,:).*VP00(deg,1,1,:);
            end
            for deg = 1:dim0
                Q0(i,j,1,:) = Q0(i,j,1,:) + C0(i,j,deg).*VPB0(deg,1,1,:);
            end
        end
    end
    
    for i = 1:N
        for j = 1:N
            Qr(i,j,1,1) = f1(Xc(i) - tend,Yc(j) - tend);
            Qr(i,j,1,2) = f2(Xc(i) - tend,Yc(j) - tend);
            Qr(i,j,1,3) = f3(Xc(i) - tend,Yc(j) - tend);
            Qr(i,j,1,4) = f4(Xc(i) - tend,Yc(j) - tend);
            Qr0(i,j,1,1) = f5(Xc(i) - tend,Yc(j) - tend);
            Qr0(i,j,1,2) = f6(Xc(i) - tend,Yc(j) - tend);
        end
    end
    
    E = abs(Q - Qr);
    E0 = abs(Q0 - Qr0);
    QF(:,:,end + 1,:) = Q;
    QF0(:,:,end + 1,:) = Q0;
    
    for i = 1:N
        for j = 1:N
            if X(i) < -5 + tend || X(i) > 5 + tend || Y(j) < -5 + tend || Y(j) > 5 + tend
                E(i,j,1,:) = 0;
                E0(i,j,1,:) = 0;
            end
        end
    end
    
    %Linfty = max(max(max(max(E))));
    L2 = sqrt(sum(sum(sum(sum(E0(:,:,1,1).^2)))))/N;
    
    fprintf('%d  %d  %d  %d %d\n',tend,tb,L2,t0,(tb - tend)/t*t0)
    
end


Q = zeros(N,N,1,dim2);
Q0 = zeros(N,N,1,2);

for i = 1:N
    for j = 1:N
        for deg = 1:dim1
            Q(i,j,1,:) = Q(i,j,1,:) + C(i,j,deg,:).*VP00(deg,1,1,:);
        end
        for deg = 1:dim0
            Q0(i,j,1,:) = Q0(i,j,1,:) + C0(i,j,deg).*VPB0(deg,1,1,:);
        end
        Qr(i,j,1,1) = f1(Xc(i) - tb,Yc(j) - tb);
        Qr(i,j,1,2) = f2(Xc(i) - tb,Yc(j) - tb);
        Qr(i,j,1,3) = f3(Xc(i) - tb,Yc(j) - tb);
        Qr(i,j,1,4) = f4(Xc(i) - tb,Yc(j) - tb);
        Qr0(i,j,1,1) = f5(Xc(i) - tb,Yc(j) - tb);
        Qr0(i,j,1,2) = f6(Xc(i) - tb,Yc(j) - tb);
    end
end

E = abs(Q - Qr);
E0 = abs(Q0 - Qr0);

for i = 1:N
    for j = 1:N
        if X(i) < -5 + tb || X(i) > 5 + tb || Y(j) < -5 + tb || Y(j) > 5 + tb
            E(i,j,1,:) = 0;
            E0(i,j,1,:) = 0;
        end
    end
end

flash2DB
erB
