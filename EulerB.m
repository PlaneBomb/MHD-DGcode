function [DC,DC0] = EulerB(C,C0,h,VP,VPdx,VPdy,VPpd,VPB,DxVPB,DyVPB,VPBpd,lambdai,weight)

[N,~,dim1,dim2] = size(C);
[~,~,dim0] = size(C0);
N = N - 2;
quad = length(lambdai);
DC = zeros(N + 2,N + 2,dim1,dim2);
DC0 = zeros(N + 2,N + 2,dim0);
h1 = h/2;

for i = 2:N + 1
    for j = 2:N + 1
        for deg = 1:dim1
            
            % 体积分
            II = zeros(1,1,1,dim2);
            II0 = 0;
            for i1 = 1:quad
                for j1 = 1:quad
                    uint = zeros(1,1,1,dim2);
                    uBint = zeros(1,1,1,2);
                    for degg = 1:dim1
                        uint = uint + C(i,j,degg,:).*VP(i1,j1,degg,:);
                    end
                    for degg = 1:dim0
                        uBint = uBint + C0(i,j,degg).*VPB(i1,j1,degg,:);
                    end
                    [y1,y2,yB1,yB2] = f(uint,uBint);
                    II = II + weight(i1)*weight(j1)*(y1.*VPdx(i1,j1,deg,:) + y2.*VPdy(i1,j1,deg,:));
                    II0 = II0 + weight(i1)*weight(j1)*(yB1(1,1,1,1)*DxVPB(i1,j1,deg,1) + yB1(1,1,1,2)*DxVPB(i1,j1,deg,2) + yB2(1,1,1,1)*DyVPB(i1,j1,deg,1) + yB2(1,1,1,2)*DyVPB(i1,j1,deg,2));
                    
                end
            end
            II = II*h1^2;
            II0 = II0*h1^2;
            % 面积分
            % 右
            IR = zeros(1,1,1,dim2);
            IR0 = 0;
            for j1 = 1:quad
                uint = zeros(1,1,1,dim2);
                uext = zeros(1,1,1,dim2);
                uBint = zeros(1,1,1,2);
                uBext = zeros(1,1,1,2);
                for degg = 1:dim1
                    uint = uint + C(i,j,degg,:).*VPpd(1,j1,degg,:);
                    uext = uext + C(i + 1,j,degg,:).*VPpd(2,j1,degg,:);
                end
                for degg = 1:dim0
                    uBint = uBint + C0(i,j,degg).*VPBpd(1,j1,degg,:);
                    uBext = uBext + C0(i + 1,j,degg).*VPBpd(2,j1,degg,:);
                end
                alphax = cxB(uext,uBext);
                [y11,~,yB11,~] = f(uint,uBint);
                [y21,~,yB21,~] = f(uext,uBext);
                IR = IR + weight(j1)*(y11 + y21 - alphax*(uext - uint)).*VPpd(1,j1,deg,:)/2;
                IR0 = IR0 + weight(j1)*(yB11 + yB21 - alphax*(uBext - uBint)).*VPBpd(1,j1,deg,:)/2;
            end
            IR = IR*h1;
            IR0 = IR0(1,1,1,1) + IR0(1,1,1,2);
            IR0 = IR0*h1;
            
            % 左
            IL = zeros(1,1,1,dim2);
            IL0 = 0;
            for j1 = 1:quad
                uint = zeros(1,1,1,dim2);
                uext = zeros(1,1,1,dim2);
                uBint = zeros(1,1,1,2);
                uBext = zeros(1,1,1,2);
                for degg = 1:dim1
                    uint = uint + C(i,j,degg,:).*VPpd(2,j1,degg,:);
                    uext = uext + C(i - 1,j,degg,:).*VPpd(1,j1,degg,:);
                end
                for degg = 1:dim0
                    uBint = uBint + C0(i,j,degg).*VPBpd(2,j1,degg,:);
                    uBext = uBext + C0(i - 1,j,degg).*VPBpd(1,j1,degg,:);
                end
                alphax = cxB(uext,uBext);
                [y11,~,yB11,~] = f(uint,uBint);
                [y21,~,yB21,~] = f(uext,uBext);
                IL = IL + weight(j1)*(-y11 - y21 - alphax*(uext - uint)).*VPpd(2,j1,deg,:)/2;
                IL0 = IL0 + weight(j1)*(-yB11 - yB21 - alphax*(uBext - uBint)).*VPBpd(2,j1,deg,:)/2;
            end
            IL = IL*h1;
            IL0 = IL0(1,1,1,1) + IL0(1,1,1,2);
            IL0 = IL0*h1;
            
            % 上
            IU = zeros(1,1,1,dim2);
            IU0 = 0;
            for i1 = 1:quad
                uint = zeros(1,1,1,dim2);
                uext = zeros(1,1,1,dim2);
                uBint = zeros(1,1,1,2);
                uBext = zeros(1,1,1,2);
                for degg = 1:dim1
                    uint = uint + C(i,j,degg,:).*VPpd(3,i1,degg,:);
                    uext = uext + C(i,j + 1,degg,:).*VPpd(4,i1,degg,:);
                end
                for degg = 1:dim0
                    uBint = uBint + C0(i,j,degg).*VPBpd(3,i1,degg,:);
                    uBext = uBext + C0(i,j + 1,degg).*VPBpd(4,i1,degg,:);
                end
                alphay = cyB(uext,uBext);
                [~,y12,~,yB12] = f(uint,uBint);
                [~,y22,~,yB22] = f(uext,uBext);
                IU = IU + weight(i1)*(y12 + y22 - alphay*(uext - uint)).*VPpd(3,i1,deg,:)/2;
                IU0 = IU0 + weight(i1)*(yB12 + yB22 - alphay*(uBext - uBint)).*VPBpd(3,i1,deg,:)/2;
            end
            IU = IU*h1;
            IU0 = IU0(1,1,1,1) + IU0(1,1,1,2);
            IU0 = IU0*h1;
            
            % 下
            ID = zeros(1,1,1,dim2);
            ID0 = 0;
            for i1 = 1:quad
                uint = zeros(1,1,1,dim2);
                uext = zeros(1,1,1,dim2);
                uBint = zeros(1,1,1,2);
                uBext = zeros(1,1,1,2);
                for degg = 1:dim1
                    uint = uint + C(i,j,degg,:).*VPpd(4,i1,degg,:);
                    uext = uext + C(i,j - 1,degg,:).*VPpd(3,i1,degg,:);
                end
                for degg = 1:dim0
                    uBint = uBint + C0(i,j,degg).*VPBpd(4,i1,degg,:);
                    uBext = uBext + C0(i,j - 1,degg).*VPBpd(3,i1,degg,:);
                end
                alphay = cyB(uext,uBext);
                [~,y12,~,yB12] = f(uint,uBint);
                [~,y22,~,yB22] = f(uext,uBext);
                ID = ID + weight(i1)*(-y12 - y22 - alphay*(uext - uint)).*VPpd(4,i1,deg,:)/2;
                ID0 = ID0 + weight(i1)*(-yB12 - yB22 - alphay*(uBext - uBint)).*VPBpd(4,i1,deg,:)/2;
            end
            ID = ID*h1;
            ID0 = ID0(1,1,1,1) + ID0(1,1,1,2);
            ID0 = ID0*h1;
            
            DC(i,j,deg,:) = II - IR - IL - IU - ID;
            DC0(i,j,deg) = II0 - IR0 - IL0 - IU0 - ID0;
        end
        
        
        
        for deg = dim1 + 1:dim0
            
            II0 = 0;
            % 体积分
            for i1 = 1:quad
                for j1 = 1:quad
                    uint = zeros(1,1,1,dim2);
                    uBint = zeros(1,1,1,2);
                    for degg = 1:dim1
                        uint = uint + C(i,j,degg,:).*VP(i1,j1,degg,:);
                    end
                    for degg = 1:dim0
                        uBint = uBint + C0(i,j,degg).*VPB(i1,j1,degg,:);
                    end
                    [~,~,yB1,yB2] = f(uint,uBint);
                    II0 = II0 + weight(i1)*weight(j1)*(yB1(1,1,1,1)*DxVPB(i1,j1,deg,1) + yB1(1,1,1,2)*DxVPB(i1,j1,deg,2) + yB2(1,1,1,1)*DyVPB(i1,j1,deg,1) + yB2(1,1,1,2)*DyVPB(i1,j1,deg,2));
                end
            end
            II0 = II0*h1^2;
            % 面积分
            % 右
            IR0 = 0;
            for j1 = 1:quad
                uint = zeros(1,1,1,dim2);
                uext = zeros(1,1,1,dim2);
                uBint = zeros(1,1,1,2);
                uBext = zeros(1,1,1,2);
                for degg = 1:dim1
                    uint = uint + C(i,j,degg,:).*VPpd(1,j1,degg,:);
                    uext = uext + C(i + 1,j,degg,:).*VPpd(2,j1,degg,:);
                end
                for degg = 1:dim0
                    uBint = uBint + C0(i,j,degg).*VPBpd(1,j1,degg,:);
                    uBext = uBext + C0(i + 1,j,degg).*VPBpd(2,j1,degg,:);
                end
                alphax = cxB(uext,uBext);
                [~,~,yB11,~] = f(uint,uBint);
                [~,~,yB21,~] = f(uext,uBext);
                IR0 = IR0 + weight(j1)*(yB11 + yB21 - alphax*(uBext - uBint)).*VPBpd(1,j1,deg,:)/2;
                
            end
            IR0 = IR0(1,1,1,1) + IR0(1,1,1,2);
            IR0 = IR0*h1;
            
            
            % 左
            IL0 = 0;
            for j1 = 1:quad
                uint = zeros(1,1,1,dim2);
                uext = zeros(1,1,1,dim2);
                uBint = zeros(1,1,1,2);
                uBext = zeros(1,1,1,2);
                for degg = 1:dim1
                    uint = uint + C(i,j,degg,:).*VPpd(2,j1,degg,:);
                    uext = uext + C(i - 1,j,degg,:).*VPpd(1,j1,degg,:);
                end
                for degg = 1:dim0
                    uBint = uBint + C0(i,j,degg).*VPBpd(2,j1,degg,:);
                    uBext = uBext + C0(i - 1,j,degg).*VPBpd(1,j1,degg,:);
                end
                alphax = cxB(uext,uBext);
                [~,~,yB11,~] = f(uint,uBint);
                [~,~,yB21,~] = f(uext,uBext);
                IL0 = IL0 + weight(j1)*(-yB11 - yB21 - alphax*(uBext - uBint)).*VPBpd(2,j1,deg,:)/2;
                
            end
            IL0 = IL0(1,1,1,1) + IL0(1,1,1,2);
            IL0 = IL0*h1;
            
            % 上
            IU0 = 0;
            for i1 = 1:quad
                uint = zeros(1,1,1,dim2);
                uext = zeros(1,1,1,dim2);
                uBint = zeros(1,1,1,2);
                uBext = zeros(1,1,1,2);
                for degg = 1:dim1
                    uint = uint + C(i,j,degg,:).*VPpd(3,i1,degg,:);
                    uext = uext + C(i,j + 1,degg,:).*VPpd(4,i1,degg,:);
                end
                for degg = 1:dim0
                    uBint = uBint + C0(i,j,degg).*VPBpd(3,i1,degg,:);
                    uBext = uBext + C0(i,j + 1,degg).*VPBpd(4,i1,degg,:);
                end
                alphay = cyB(uext,uBext);
                [~,~,~,yB12] = f(uint,uBint);
                [~,~,~,yB22] = f(uext,uBext);
                IU0 = IU0 + weight(i1)*(yB12 + yB22 - alphay*(uBext - uBint)).*VPBpd(3,i1,deg,:)/2;
                
            end
            IU0 = IU0(1,1,1,1) + IU0(1,1,1,2);
            IU0 = IU0*h1;
            
            % 下
            ID0 = 0;
            for i1 = 1:quad
                uint = zeros(1,1,1,dim2);
                uext = zeros(1,1,1,dim2);
                uBint = zeros(1,1,1,2);
                uBext = zeros(1,1,1,2);
                for degg = 1:dim1
                    uint = uint + C(i,j,degg,:).*VPpd(4,i1,degg,:);
                    uext = uext + C(i,j - 1,degg,:).*VPpd(3,i1,degg,:);
                end
                for degg = 1:dim0
                    uBint = uBint + C0(i,j,degg).*VPBpd(4,i1,degg,:);
                    uBext = uBext + C0(i,j - 1,degg).*VPBpd(3,i1,degg,:);
                end
                alphay = cyB(uext,uBext);
                [~,~,~,yB12] = f(uint,uBint);
                [~,~,~,yB22] = f(uext,uBext);
                ID0 = ID0 + weight(i1)*(-yB12 - yB22 - alphay*(uBext - uBint)).*VPBpd(4,i1,deg,:)/2;
                
            end
            ID0 = ID0(1,1,1,1) + ID0(1,1,1,2);
            ID0 = ID0*h1;
            
            DC0(i,j,deg) = II0 - IR0 - IL0 - IU0 - ID0;
            
        end
    end
end

DC(1,:,:,:) = DC(N + 1,:,:,:);
DC(:,1,:,:) = DC(:,N + 1,:,:);
DC(N + 2,:,:,:) = DC(2,:,:,:);
DC(:,N + 2,:,:) = DC(:,2,:,:);

DC0(1,:,:) = DC0(N + 1,:,:);
DC0(:,1,:) = DC0(:,N + 1,:);
DC0(N + 2,:,:) = DC0(2,:,:);
DC0(:,N + 2,:) = DC0(:,2,:);

end



    