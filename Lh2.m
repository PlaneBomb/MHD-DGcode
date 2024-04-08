function du=Lh2(uh)
global point weight phi phix phiy hh m phiL phiR phiU phiD

%%
%period BC
[N,~,~,~,]=size(uh);
uhh=zeros(N+2,N+2,6,6);

uhh(2:end-1,2:end-1,:,:)=uh;
uhh(1,2:end-1,:,:)=uh(end,:,:,:);
uhh(end,2:end-1,:,:)=uh(1,:,:,:);
uhh(2:end-1,1,:,:)=uh(:,end,:,:);
uhh(2:end-1,end,:,:)=uh(:,1,:,:);



%value on gauss point
uhG=ValueGausspoint(uh);

du=zeros(size(uh));
%integral on faces
for i=1:N
    for j=1:N
        for p=1:4
            for q=1:4
                FF=Ff(uhG(i,j,p,q,:));
                F1=FF(:,1);
                F2=FF(:,2);
                for k=1:6
                    for dim=1:6
                        du(i,j,k,dim)=du(i,j,k,dim)+weight(p)*weight(q)*(F1(dim)*phix(p,q,k)+F2(dim)*phiy(p,q,k))*hh^2;
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
        uint=zeros(1,1,1,4,6);
        uext=zeros(1,1,1,4,6);

        for p=1:4
            for dim=1:6
            uint(1,1,1,p,dim)=uh(i,j,1,dim)*phiU(p,1)+uh(i,j,2,dim)*phiU(p,2)+uh(i,j,3,dim)*phiU(p,3)+uh(i,j,4,dim)*phiU(p,4)+uh(i,j,5,dim)*phiU(p,5)+uh(i,j,6,dim)*phiU(p,6);
            uext(1,1,1,p,dim)=uhh(i+2,j+1,1,dim)*phiD(p,1)+uhh(i+2,j+1,2,dim)*phiD(p,2)+uhh(i+2,j+1,3,dim)*phiD(p,3)+uhh(i+2,j+1,4,dim)*phiD(p,4)+uhh(i+2,j+1,5,dim)*phiD(p,5)+uhh(i+2,j+1,6,dim)*phiD(p,6);
            end
        end
        for p=1:4
            %F2
            FF=Ff(uint(1,1,1,p,:));
            F2=FF(:,2);

            FF=Ff(uext(1,1,1,p,:));
            F22=FF(:,2);
            [~,alphay]=Alphaa(uint(1,1,1,p,:),1,1);
            [~,alphayy]=Alphaa(uext(1,1,1,p,:),1,1);
            alphay=max(alphay,alphayy);
            for k=1:6
                for dim=1:6
                    du(i,j,k,dim)=du(i,j,k,dim)-weight(p)*(F2(dim)+F22(dim)-alphay*(uext(1,1,1,p,dim)-uint(1,1,1,p,dim)))*phiU(p,k)/2*hh;
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
        uint=zeros(1,1,1,4,6);
        uext=zeros(1,1,1,4,6);

        for p=1:4
            for dim=1:6
            uint(1,1,1,p,dim)=uh(i,j,1,dim)*phiD(p,1)+uh(i,j,2,dim)*phiD(p,2)+uh(i,j,3,dim)*phiD(p,3)+uh(i,j,4,dim)*phiD(p,4)+uh(i,j,5,dim)*phiD(p,5)+uh(i,j,6,dim)*phiD(p,6);
            uext(1,1,1,p,dim)=uhh(i,j+1,1,dim)*phiU(p,1)+uhh(i,j+1,2,dim)*phiU(p,2)+uhh(i,j+1,3,dim)*phiU(p,3)+uhh(i,j+1,4,dim)*phiU(p,4)+uhh(i,j+1,5,dim)*phiU(p,5)+uhh(i,j+1,6,dim)*phiU(p,6);
            end
        end
        for p=1:4
            %F2
            FF=Ff(uint(1,1,1,p,:));
            F2=FF(:,2);

            FF=Ff(uext(1,1,1,p,:));
            F22=FF(:,2);
            [~,alphay]=Alphaa(uint(1,1,1,p,:),1,1);
            [~,alphayy]=Alphaa(uext(1,1,1,p,:),1,1);
            alphay=max(alphay,alphayy);
            for k=1:6
                for dim=1:6
                du(i,j,k,dim)=du(i,j,k,dim)-weight(p)*(-F2(dim)-F22(dim)-alphay*(uext(1,1,1,p,dim)-uint(1,1,1,p,dim)))*phiD(p,k)/2*hh;
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
        uint=zeros(1,1,1,4,6);
        uext=zeros(1,1,1,4,6);

        for p=1:4
            for dim=1:6
            uint(1,1,1,p,dim)=uh(i,j,1,dim)*phiL(p,1)+uh(i,j,2,dim)*phiL(p,2)+uh(i,j,3,dim)*phiL(p,3)+uh(i,j,4,dim)*phiL(p,4)+uh(i,j,5,dim)*phiL(p,5)+uh(i,j,6,dim)*phiL(p,6);
            uext(1,1,1,p,dim)=uhh(i+1,j,1,dim)*phiR(p,1)+uhh(i+1,j,2,dim)*phiR(p,2)+uhh(i+1,j,3,dim)*phiR(p,3)+uhh(i+1,j,4,dim)*phiR(p,4)+uhh(i+1,j,5,dim)*phiR(p,5)+uhh(i+1,j,6,dim)*phiR(p,6);
            end
        end
        for p=1:4
            %F1
            FF=Ff(uint(1,1,1,p,:));
            F1=FF(:,1);

            FF=Ff(uext(1,1,1,p,:));
            F11=FF(:,1);
            [alphax,~]=Alphaa(uint(1,1,1,p,:),1,1);
            [alphaxx,~]=Alphaa(uext(1,1,1,p,:),1,1);
            alphax=max(alphax,alphaxx);

            for k=1:6
                for dim=1:6
                du(i,j,k,dim)=du(i,j,k,dim)-weight(p)*(-F1(dim)-F11(dim)-alphax*(uext(1,1,1,p,dim)-uint(1,1,1,p,dim)))*phiL(p,k)/2*hh;
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
        uint=zeros(1,1,1,4,6);
        uext=zeros(1,1,1,4,6);

        for p=1:4
            for dim=1:6
            uint(1,1,1,p,dim)=uh(i,j,1,dim)*phiR(p,1)+uh(i,j,2,dim)*phiR(p,2)+uh(i,j,3,dim)*phiR(p,3)+uh(i,j,4,dim)*phiR(p,4)+uh(i,j,5,dim)*phiR(p,5)+uh(i,j,6,dim)*phiR(p,6);
            uext(1,1,1,p,dim)=uhh(i+1,j+2,1,dim)*phiL(p,1)+uhh(i+1,j+2,2,dim)*phiL(p,2)+uhh(i+1,j+2,3,dim)*phiL(p,3)+uhh(i+1,j+2,4,dim)*phiL(p,4)+uhh(i+1,j+2,5,dim)*phiL(p,5)+uhh(i+1,j+2,6,dim)*phiL(p,6);
            end
        end
        for p=1:4
            %F1
            FF=Ff(uint(1,1,1,p,:));
            F1=FF(:,1);

            FF=Ff(uext(1,1,1,p,:));
            F11=FF(:,1);
            [alphax,~]=Alphaa(uint(1,1,1,p,:),1,1);
            [alphaxx,~]=Alphaa(uext(1,1,1,p,:),1,1);
            alphax=max(alphax,alphaxx);
            for k=1:6
                for dim=1:6
                du(i,j,k,dim)=du(i,j,k,dim)-weight(p)*(F1(dim)+F11(dim)-alphax*(uext(1,1,1,p,dim)-uint(1,1,1,p,dim)))*phiR(p,k)/2*hh;
                end
            end
        end
    end
end

for i=1:6
    du(:,:,i,:)=du(:,:,i,:)/m(i);
end
du=du/hh^2;