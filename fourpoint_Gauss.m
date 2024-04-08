function [point,weight]=fourpoint_Gauss(~)
    point=[-sqrt(525+70*sqrt(30))/35,-sqrt(525-70*sqrt(30))/35,sqrt(525-70*sqrt(30))/35,sqrt(525+70*sqrt(30))/35];
    weight=[(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];
end