function uhGB = ValueGausspointB(uhB)
global psi
[N,~,dim1,~]=size(uhB);
uhGB=zeros(N,N,4,4,2);
for i=1:N
    for j=1:N
        for p=1:4
            for q=1:4
                for dim=1:2
                    for k=1:dim1
                        uhGB(i,j,p,q,dim)=uhGB(i,j,p,q,dim)+uhB(i,j,k,1)*psi(p,q,k,dim);
                    end
                end
            end
        end
    end
end    
end