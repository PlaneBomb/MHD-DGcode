function uhG = ValueGausspoint(uh)
global phi
[N,~,dim1,dim2]=size(uh);
uhG=zeros(N,N,4,4,dim2);
for i=1:N
    for j=1:N
        for p=1:4
            for q=1:4
                for k=1:dim2
                    for dim=1:dim1
                        uhG(i,j,p,q,k)=uhG(i,j,p,q,k)+uh(i,j,dim,k)*phi(p,q,dim);  
                    end
                end
            end
        end
    end
end    
end

