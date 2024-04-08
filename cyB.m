
function alpha = cyB(U,UB)
x = [0;0;0;0];
xB = [0;0];
for i = 1:4
    x(i) = U(1,1,1,i);
end
for i = 1:2
    xB(i) = UB(1,1,1,i);
end
p = (x(4) - 0.5*(x(2).^2 + x(3).^2)./x(1) - 0.5*(xB(1).^2 + xB(2).^2))*(2/3);
a = sqrt((5/3)*p/x(1));
cf = (1/2)*sqrt(a.^2 + (xB(1).^2 + xB(2).^2)/x(1) + sqrt((a.^2 + (xB(1).^2 + xB(2).^2)/x(1)).^2 - 4*a.^2*(xB(2)).^2/x(1)));
alpha = abs(x(3)/x(1)) + abs(cf);
end