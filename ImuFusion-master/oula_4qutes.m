% clear all;
close all;
% clc;

%欧拉角
x = 30/180*pi;
y = 0;
z = 0;
Ang1 = [x y z];

%欧拉角转旋转矩阵
Rx = [1      0      0;
    0 cos(x) -sin(x);
    0 sin(x) cos(x)];
Ry = [cos(y)  0 sin(y);
    0       1      0;
    -sin(y) 0 cos(y)];
Rz = [cos(z) -sin(z) 0;
    sin(z) cos(z)  0;
    0      0       1];
R = Rz*Ry*Rx;
R1 = R;

%旋转矩阵转欧拉角
x = atan2(R(3,2),R(3,3));
y = atan2(-R(3,1),sqrt(R(3,2)^2+R(3,3)^2));
z = atan2(R(2,1),R(1,1));
Ang2 = [x y z];

%旋转矩阵转四元数
t=sqrt(1+R(1,1)+R(2,2)+R(3,3))/2;
q=[t (R(3,2)-R(2,3))/(4*t) (R(1,3)-R(3,1))/(4*t) (R(2,1)-R(1,2))/(4*t)];
Q1 = q;


%四元数转旋转矩阵
R=[ 2*q(1).^2-1+2*q(2)^2    2*(q(2)*q(3)-q(1)*q(4)) 2*(q(2)*q(4)+q(1)*q(3));
    2*(q(2)*q(3)+q(1)*q(4)) 2*q(1)^2-1+2*q(3)^2     2*(q(3)*q(4)-q(1)*q(2));
    2*(q(2)*q(4)-q(1)*q(3)) 2*(q(3)*q(4)+q(1)*q(2)) 2*q(1)^2-1+2*q(4)^2];
R2 = R;

%欧拉角转四元数
q = [cos(x/2)*cos(y/2)*cos(z/2) + sin(x/2)*sin(y/2)*sin(z/2) ...
    sin(x/2)*cos(y/2)*cos(z/2) - cos(x/2)*sin(y/2)*sin(z/2) ...
    cos(x/2)*sin(y/2)*cos(z/2) + sin(x/2)*cos(y/2)*sin(z/2) ...
    cos(x/2)*cos(y/2)*sin(z/2) - sin(x/2)*sin(y/2)*cos(z/2)];
Q2 = q;


%四元数转欧拉角
x = atan2(2*(q(1)*q(2)+q(3)*q(4)),1 - 2*(q(2)^2+q(3)^2));
y = asin(2*(q(1)*q(3) - q(2)*q(4)));
z = atan2(2*(q(1)*q(4)+q(2)*q(3)),1 - 2*(q(3)^2+q(4)^2));
Ang3 = [x y z];


