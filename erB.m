%误差分析
Q = QF(:,:,end,:);
Q0 = QF0(:,:,end,:);
Qr = zeros(N,N,1,dim2);
Qr0 = zeros(N,N,1,2);

if tb == 20
    
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
    
else
    
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
    
end

E = abs(Q - Qr);
E0 = abs(Q0 - Qr0);

if tb == 20
    
    for i = 1:N
        for j = 1:N
            if X(i) < -5 || X(i) > 5 || Y(j) < -5 || Y(j) > 5
                E(i,j,1,:) = 0;
                E0(i,j,1,:) = 0;
            end
        end
    end
    
else
    
    for i = 1:N
        for j = 1:N
            if X(i) < -5 + tb || X(i) > 5 + tb || Y(j) < -5 + tb || Y(j) > 5 + tb
                E(i,j,1,:) = 0;
                E0(i,j,1,:) = 0;
            end
        end
    end
    
end

Linfty1 = max(max(max(max(E(:,:,1,1)))));
Linfty2 = max(max(max(max(E(:,:,1,2)))));
Linfty3 = max(max(max(max(E(:,:,1,3)))));
Linfty4 = max(max(max(max(E(:,:,1,4)))));
Linfty5 = max(max(max(max(E0(:,:,1,1)))));
Linfty6 = max(max(max(max(E0(:,:,1,2)))));
L21 = sqrt(sum(sum(sum(sum(E(:,:,1,1).^2)))))/N;
L22 = sqrt(sum(sum(sum(sum(E(:,:,1,2).^2)))))/N;
L23 = sqrt(sum(sum(sum(sum(E(:,:,1,3).^2)))))/N;
L24 = sqrt(sum(sum(sum(sum(E(:,:,1,4).^2)))))/N;
L25 = sqrt(sum(sum(sum(sum(E0(:,:,1,1).^2)))))/N;
L26 = sqrt(sum(sum(sum(sum(E0(:,:,1,2).^2)))))/N;