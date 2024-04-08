function [du,duB]=Lh(uh,uhB)
global  weight phiD phiU phiL phiR phix phiy hh psix psiy psiD psiU psiL psiR m M
%%
% period BC
[N,~,~,~]=size(uh);
uhh=zeros(N+2,N+2,6,4);
uhhB=zeros(N+2,N+2,9,1);

uhh(2:end-1,2:end-1,:,:)=uh;
uhhB(2:end-1,2:end-1,:,:)=uhB;

uhh(1,2:end-1,:,:)=uh(end,:,:,:);
uhh(end,2:end-1,:,:)=uh(1,:,:,:);
uhh(2:end-1,1,:,:)=uh(:,end,:,:);
uhh(2:end-1,end,:,:)=uh(:,1,:,:);

uhhB(1,2:end-1,:,:)=uhB(end,:,:,:);
uhhB(end,2:end-1,:,:)=uhB(1,:,:,:);
uhhB(2:end-1,1,:,:)=uhB(:,end,:,:);
uhhB(2:end-1,end,:,:)=uhB(:,1,:,:);

%value on gauss point
uhG=ValueGausspoint(uh);
uhGB=ValueGausspointB(uhB);

du=zeros(size(uh));
duB=zeros(size(uhB));


%integral on faces
for i=1:N
    for j=1:N
        for p=1:4
            for q=1:4
                FF=F(uhG(i,j,p,q,:),uhGB(i,j,p,q,:));
                %disp(uhG(i,j,p,q,:))
                F1=FF(1:4,1);
                F2=FF(1:4,2);
                for k=1:6
                    for dim=1:4
                        du(i,j,k,dim)=du(i,j,k,dim)+weight(p)*weight(q)*(F1(dim)*phix(p,q,k)+F2(dim)*phiy(p,q,k))*hh^2;
                    end
                end
                FB1=FF(5:6,1);
                FB2=FF(5:6,2);
                for k=1:9
                    for dim=1:2
                        duB(i,j,k,1)=duB(i,j,k,1)+weight(p)*weight(q)*(FB1(dim)*psix(p,q,k,dim)+FB2(dim)*psiy(p,q,k,dim))*hh^2;
                    end
                end
            end
        end
    end
end


%integral on edges  with LF-flux

%above n=[0;1];
%uhh(i+2,j+1)
% _____
%|     |
%|  ^  |
%|__|__|
%|  |  |
%|     |
%|_____|
%
%uhh(i+1,j+1)=uh(i,j)
%

for i=1:N
    for j=1:N
        uint=zeros(1,1,1,4,4);
        uext=zeros(1,1,1,4,4);
        
        uintB=zeros(1,1,1,4,2);
        uextB=zeros(1,1,1,4,2);
        for p=1:4
            for dim=1:4
                for k=1:6
                    uint(1,1,1,p,dim)=uint(1,1,1,p,dim)+uh(i,j,k,dim)*phiU(p,k);
                    uext(1,1,1,p,dim)=uext(1,1,1,p,dim)+uhh(i+2,j+1,k,dim)*phiD(p,k);
                end
            end
        end
        for p=1:4
            for dim=1:2
                for k=1:9
                    uintB(1,1,1,p,dim)=uintB(1,1,1,p,dim)+uhB(i,j,k,1)*psiU(p,k,dim);
                    uextB(1,1,1,p,dim)=uextB(1,1,1,p,dim)+uhhB(i+2,j+1,k,1)*psiD(p,k,dim);
                end
            end
        end
        for p=1:4
            %F2
            FF=F(uint(1,1,1,p,:),uintB(1,1,1,p,:));
            F2=FF(1:4,2);
            FB2=FF(5:6,2);

            FF=F(uext(1,1,1,p,:),uextB(1,1,1,p,:));
            F22=FF(1:4,2);
            FB22=FF(5:6,2);
            [~,alphay]=Alpha(uint(1,1,1,p,:),uintB(1,1,1,p,:),1,1);
            [~,alphayy]=Alpha(uext(1,1,1,p,:),uextB(1,1,1,p,:),1,1);
            alphay=max(alphay,alphayy);
            for k=1:6
                for dim=1:4
                    du(i,j,k,dim)=du(i,j,k,dim)-weight(p)*(F2(dim)+F22(dim)-alphay*(uext(1,1,1,p,dim)-uint(1,1,1,p,dim)))*phiU(p,k)/2*hh;
                end
            end
            for k=1:9
                for dim=1:2
                    duB(i,j,k,1)=duB(i,j,k,1)-weight(p)*(FB2(dim)+FB22(dim)-alphay*(uextB(1,1,1,p,dim)-uintB(1,1,1,p,dim)))*psiU(p,k,dim)/2*hh;
                end
            end
        end
    end
end




%below  n=[0;-1];
%uhh(i+1,j+1)=uh(i,j)
% _____
%|     |
%|     |
%|__|__|
%|  |  |
%|  v  |
%|_____|
%
%uhh(i,j+1)
%

for i=1:N
    for j=1:N
        %values of gauss point on edge
        %interior
        uint=zeros(1,1,1,4,4);
        uext=zeros(1,1,1,4,4);
        uintB=zeros(1,1,1,4,2);
        uextB=zeros(1,1,1,4,2);
        for p=1:4
            for dim=1:4
                for k=1:6
                    uint(1,1,1,p,dim)=uint(1,1,1,p,dim)+uh(i,j,k,dim)*phiD(p,k);
                    uext(1,1,1,p,dim)=uext(1,1,1,p,dim)+uhh(i,j+1,k,dim)*phiU(p,k);
                end
            end
        end
        for p=1:4
            for dim=1:2
                for k=1:9
                    uintB(1,1,1,p,dim)=uintB(1,1,1,p,dim)+uhB(i,j,k,1)*psiD(p,k,dim);
                    uextB(1,1,1,p,dim)=uextB(1,1,1,p,dim)+uhhB(i,j+1,k,1)*psiU(p,k,dim);
                end
            end
        end
        for p=1:4
            %F2
            FF=F(uint(1,1,1,p,:),uintB(1,1,1,p,:));
            F2=FF(1:4,2);
            FB2=FF(5:6,2);

            FF=F(uext(1,1,1,p,:),uextB(1,1,1,p,:));
            F22=FF(1:4,2);
            FB22=FF(5:6,2);
            [~,alphay]=Alpha(uint(1,1,1,p,:),uintB(1,1,1,p,:),1,1);
            [~,alphayy]=Alpha(uext(1,1,1,p,:),uextB(1,1,1,p,:),1,1);
            alphay=max(alphay,alphayy);
            for k=1:6
                for dim=1:4
                    du(i,j,k,dim)=du(i,j,k,dim)-weight(p)*(-F2(dim)-F22(dim)-alphay*(uext(1,1,1,p,dim)-uint(1,1,1,p,dim)))*phiD(p,k)/2*hh;
                end
            end
            for k=1:9
                for dim=1:2
                    duB(i,j,k,1)=duB(i,j,k,1)-weight(p)*(-FB2(dim)-FB22(dim)-alphay*(uextB(1,1,1,p,dim)-uintB(1,1,1,p,dim)))*psiD(p,k,dim)/2*hh;
                end
            end
        end
    end
end

%left n=[-1;0];
%uhh(i+1,j)
% ____ ____
%|    |    |
%|  <-|--  |
%|____|____|
%
%    uhh(i+1,j+1)=uh(i,j)
%

for i=1:N
    for j=1:N
        %values of gauss point on edge
        %interior
        uint=zeros(1,1,1,4,4);
        uext=zeros(1,1,1,4,4);
        uintB=zeros(1,1,1,4,2);
        uextB=zeros(1,1,1,4,2);
        for p=1:4
            for dim=1:4
                 for k=1:6
                    uint(1,1,1,p,dim)=uint(1,1,1,p,dim)+uh(i,j,k,dim)*phiL(p,k);
                    uext(1,1,1,p,dim)=uext(1,1,1,p,dim)+uhh(i+1,j,k,dim)*phiR(p,k);
                 end
            end
        end
        for p=1:4
            for dim=1:2
                for k=1:9
                    uintB(1,1,1,p,dim)=uintB(1,1,1,p,dim)+uhB(i,j,k,1)*psiL(p,k,dim);
                    uextB(1,1,1,p,dim)=uextB(1,1,1,p,dim)+uhhB(i+1,j,k,1)*psiR(p,k,dim);
                end
            end
        end
        for p=1:4
            %F1
            FF=F(uint(1,1,1,p,:),uintB(1,1,1,p,:));
            F1=FF(1:4,1);
            FB1=FF(5:6,1);
            FF=F(uext(1,1,1,p,:),uextB(1,1,1,p,:));

            F11=FF(1:4,1);
            FB11=FF(5:6,1);
            [alphax,~]=Alpha(uint(1,1,1,p,:),uintB(1,1,1,p,:),1,1);
            [alphaxx,~]=Alpha(uext(1,1,1,p,:),uextB(1,1,1,p,:),1,1);
            alphax=max(alphax,alphaxx);
            for k=1:6
                for dim=1:4
                    du(i,j,k,dim)=du(i,j,k,dim)-weight(p)*(-F1(dim)-F11(dim)-alphax*(uext(1,1,1,p,dim)-uint(1,1,1,p,dim)))*phiL(p,k)/2*hh;
                end
            end
            for k=1:9
                for dim=1:2
                    duB(i,j,k,1)=duB(i,j,k,1)-weight(p)*(-FB1(dim)-FB11(dim)-alphax*(uextB(1,1,1,p,dim)-uintB(1,1,1,p,dim)))*psiL(p,k,dim)/2*hh;
                end
            end
        end
    end
end


%right n=[1;0];
%uhh(i+1,j+1)
% ____ ____
%|    |    |
%|  --|->  |
%|____|____|
%
%    uhh(i+1,j+2)
%


for i=1:N
    for j=1:N
        %values of gauss point on edge
        %interior
        uint=zeros(1,1,1,4,4);
        uext=zeros(1,1,1,4,4);
        uintB=zeros(1,1,1,4,2);
        uextB=zeros(1,1,1,4,2);
        for p=1:4
            for dim=1:4
                 for k=1:6
                    uint(1,1,1,p,dim)=uint(1,1,1,p,dim)+uh(i,j,k,dim)*phiR(p,k);
                    uext(1,1,1,p,dim)=uext(1,1,1,p,dim)+uhh(i+1,j+2,k,dim)*phiL(p,k);
                 end
            end
        end
        for p=1:4
            for dim=1:2
                for k=1:9
                    uintB(1,1,1,p,dim)=uintB(1,1,1,p,dim)+uhB(i,j,k,1)*psiR(p,k,dim);
                    uextB(1,1,1,p,dim)=uextB(1,1,1,p,dim)+uhhB(i+1,j+2,k,1)*psiL(p,k,dim);
                end
            end
        end
        for p=1:4
            %F1
            FF=F(uint(1,1,1,p,:),uintB(1,1,1,p,:));
            F1=FF(1:4,1);
            FB1=FF(5:6,1);

            FF=F(uext(1,1,1,p,:),uextB(1,1,1,p,:));
            F11=FF(1:4,1);
            FB11=FF(5:6,1);
            [alphax,~]=Alpha(uint(1,1,1,p,:),uintB(1,1,1,p,:),1,1);
            [alphaxx,~]=Alpha(uext(1,1,1,p,:),uextB(1,1,1,p,:),1,1);
            alphax=max(alphax,alphaxx);
            for k=1:6
                for dim=1:4
                    du(i,j,k,dim)=du(i,j,k,dim)-weight(p)*(F1(dim)+F11(dim)-alphax*(uext(1,1,1,p,dim)-uint(1,1,1,p,dim)))*phiR(p,k)/2*hh;
                end
            end
            for k=1:9
                for dim=1:2
                    duB(i,j,k,1)=duB(i,j,k,1)-weight(p)*(FB1(dim)+FB11(dim)-alphax*(uextB(1,1,1,p,dim)-uintB(1,1,1,p,dim)))*psiR(p,k,dim)/2*hh;
                end
            end
        end
    end
end

for i=1:6
    du(:,:,i,:)=du(:,:,i,:)/m(i);
end
for i=1:9
    duB(:,:,i,:)=duB(:,:,i,:)/M(i);
end
du=du/hh^2;
duB=duB/hh^2;
