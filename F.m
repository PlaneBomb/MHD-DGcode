function Fu=F(uhG,uhGB)
global gamma
    F1=zeros(1,6);
    F2=zeros(1,6);
    u=zeros(1,4);
    uB=zeros(1,2);
    for i=1:4
        u(i)=uhG(1,1,1,1,i);
    end
    
    for i=1:2
        uB(i)=uhGB(1,1,1,1,i);
    end
    
    pre=(gamma-1)*(u(4)-0.5*(u(2)^2+u(3)^2)/u(1)-0.5*(uB(1)^2+uB(2)^2));
    rho=u(1);
    Bnorm=sum(uB.^2);
    Bx=uB(1);
    By=uB(2);
    
    F1(1)=u(2);
    F1(2)=pre+0.5*Bnorm-Bx^2+u(2)^2/u(1);
    F1(3)=u(2)*u(3)/u(1)-Bx*By;
    F1(4)=u(2)/u(1)*(u(4)+0.5*Bnorm+pre)-Bx*(u(2)*Bx+u(3)*By)/rho;
    
    F1(5)=0;
    F1(6)=u(2)*By/rho-u(3)*Bx/rho;

    F2(1)=u(3);
    F2(2)=u(2)*u(3)/u(1)-Bx*By;
    F2(3)=pre+0.5*Bnorm-By^2+u(3)^2/u(1);
    F2(4)=u(3)/u(1)*(u(4)+0.5*Bnorm+pre)-By*(u(2)*Bx+u(3)*By)/rho;
    
    F2(5)=u(3)*Bx/rho-u(2)*By/rho;
    F2(6)=0;

    Fu=[F1;F2]';
end