function [alphax,alphay]=Alphaa(uhG,N,quad)
global gamma
alphax=0;
alphay=0;
for i=1:N
    for j=1:N
        for p=1:quad
            for q=1:quad
                u =zeros(1,6);
                for k = 1:6
                    u(k) = uhG(i,j,p,q,k);
                end
        Bnorm=sum(u(5:6).^2);
        pre=(gamma-1)*(u(4)-0.5*(u(2)^2+u(3)^2)/u(1)-0.5*Bnorm);
        if pre<0
            exit
        end
        a=sqrt(gamma*pre/u(1));
        %disp(u(1))
        rho=u(1);
        
        Bx=u(5);
        By=u(6);
        
        cfx=sqrt(0.5*(a^2+Bnorm/rho+sqrt((a^2+Bnorm/rho)^2-4*a^2*Bx^2/rho)));
        cfy=sqrt(0.5*(a^2+Bnorm/rho+sqrt((a^2+Bnorm/rho)^2-4*a^2*By^2/rho)));

        alphax=max(alphax,abs(u(2)/u(1))+cfx);
        alphay=max(alphay,abs(u(3)/u(1))+cfy);     
            end
        end
    end
end
end