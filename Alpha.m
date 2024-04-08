function [alphax,alphay]=Alpha(uhG,uhGB,N,quad)
global gamma
alphax=0;
alphay=0;
for i=1:N
    for j=1:N
        for p=1:quad
            for q=1:quad
                u =zeros(1,4);
                uB =zeros(1,2);
                for k = 1:4
                    u(k) = uhG(i,j,p,q,k);
                end
                for k = 1:2
                    uB(k) = uhGB(i,j,p,q,k);
                end

                pre=(gamma-1)*(u(4)-0.5*(u(2)^2+u(3)^2)/u(1)-0.5*(uB(1)^2+uB(2)^2));
                a=sqrt(gamma*pre/u(1));
                rho=u(1);
                Bnorm=sum(uB.^2);
                Bx=uB(1);
                By=uB(2);

                cfx=sqrt(0.5*(a^2+Bnorm/rho+sqrt((a^2+Bnorm/rho)^2-4*a^2*Bx^2/rho)));
                cfy=sqrt(0.5*(a^2+Bnorm/rho+sqrt((a^2+Bnorm/rho)^2-4*a^2*By^2/rho)));

                alphax=max(alphax,abs(u(2)/u(1))+cfx);
                alphay=max(alphay,abs(u(3)/u(1))+cfy);
            end
        end
    end
end
end