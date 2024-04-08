%显示解随时间变化规律,flash2DB.m

h = figure();				% 创建图形窗口
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');	% 关闭相关的警告提示（因为调用了非公开接口）
jFrame = get(h,'JavaFrame');	% 获取底层 Java 结构相关句柄
pause(0.1);					% 在 Win 10，Matlab 2017b 环境下不加停顿会报 Java 底层错误。各人根据需要可以进行实验验证
set(jFrame,'Maximized',1);	%设置其最大化为真（0 为假）
pause(0.1);					% 个人实践中发现如果不停顿，窗口可能来不及变化，所获取的窗口大小还是原来的尺寸。各人根据需要可以进行实验验证
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');		% 打开相关警告设置


[xc,yc] = meshgrid(Xc,Yc);
M = size(QF0,3);
xlabel('X');
ylabel('Y');
zlabel('Q(X,Y)');
grid on;
%s = [X(1),X(end),Y(1),Y(end),min(min(min(Q(:,:,:,4)))) - 0.1,max(max(max(Q(:,:,:,4)))) + 0.1];
s = [Xc(1),Xc(end),Yc(1),Yc(end),2,3];
axis(s);
TT = 150;
t0 = T(end)/TT;
for i = 1:TT + 1
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
    mesh(xc,yc,QF(:,:,j,4));
    if max(QF0(:,:,j,1)) > 10
        break
    end
    axis(s);
    pause(0.0001);
end