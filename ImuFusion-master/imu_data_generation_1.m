% =========================================================================
%
%                  IMU_kalman_filter
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 4月14日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.生成一组数据 [m/s^2] 与 [rad/s]  ，并且用欧拉积分与中值积分
%--------------------------------------------------------------------------


clear all;
close all;
dt = 1/256;
M_PI = 3.1415926;
gn = [0, 0, -9.81];
k = 0;
time = 0:dt:20;
for t = 0:dt:20
    k = k+1;
    ellipse_x = 15;
    ellipse_y = 20;
    z = 1; 
    K1 = 10; 
    K = M_PI / 10;
    position = [ellipse_x * cos(K * t) + 5, ellipse_y * sin(K * t) + 5, z * sin(K1 * K * t) + 5];
    dp = [-K * ellipse_x * sin(K * t), K * ellipse_y * cos(K * t), z * K1 * K * cos(K1 * K * t)];
    K2 = K * K;
    ddp = [-K2 * ellipse_x * cos(K * t), -K2 * ellipse_y * sin(K * t), -z * K1 * K1 * K2 * sin(K1 * K * t)];
    %%Rotation
    k_roll = 0.1;
    k_pitch = 0.2;
    eulerAngles = [k_roll * cos(t), k_pitch * sin(t), K * t];
    eulerAnglesRates = [-k_roll * sin(t), k_pitch * cos(t), K];
    
    Rwb = euler2Rotation(eulerAngles);
    
    imu_gyro = eulerRates2bodyRates(eulerAngles) * eulerAnglesRates';
    imu_acc = Rwb'* (ddp - gn)';
    
    
    %%
    vel(k,1:6) = [dp,ddp];
    pos(k,1:6) = [position , eulerAngles];
    imu_data(1:6,k) = [imu_acc;imu_gyro];
end

figure;plot(pos(:,1:3));title('位置')
figure;plot(pos(:,4:6));title('角度')


%% way1 
    gyroX = imu_data(4,:);
    gyroY = imu_data(5,:);
    gyroZ = imu_data(6,:);
    gyro =  [gyroX ; gyroY ; gyroZ];
    gyro = gyro';
    angle = pos(1,4:6);
    for t = 2:length(gyroX)
        angle(t,:) = angle(t-1,:)+ dt*gyro(t,:);
        for ind = 1:3
            if angle(t,ind)>pi
                angle(t,ind) = angle(t,ind)-2*pi;
            elseif angle(t,ind)<-pi
                angle(t,ind) = angle(t,ind)+2*pi;
            end
        end
    end
    
figure;plot(angle-pos(:,4:6));
legend 1 2

    

%% 中值积分
gw = [0, 0, -9.81]';
    % init_a = mean(imu_data(1:3,1),2);
    % init_a = init_a / norm(init_a);
    % init_psi =  0;
    % init_theta = -asin(init_a(1));
    % init_phi = atan2(init_a(2), init_a(3));
    % init_quat = angle2quat(init_psi, init_theta, init_phi); %%初始角度z,y,x
init_ang = pos(1,4:6);
init_quat = angle2quat(init_ang(3),init_ang(2),init_ang(1));
 
init_quat_tp = init_quat;     %%初始角度赋值
quat = init_quat_tp';
N = length(imu_data);
Vw = vel(1,1:3)';
Pwb = pos(1,1:3)';
for k = 2:N
    w = (imu_data(4:6, k-1)    +    imu_data(4:6, k))/2;
    a = imu_data(1:3, k-1);
    a1 = imu_data(1:3, k);
    Ow = [0     -w(1)   -w(2)    -w(3);...
          w(1)   0       w(3)    -w(2);...
          w(2)  -w(3)    0        w(1);...
          w(3)   w(2)   -w(1)     0  ];
    Fc = 0.5*Ow;
    F = eye(4) + Fc*dt;%%或者在这里改变积分的方法(对于矩阵的积分方法)
    quat1 = F* quat;
    
    
    %%
    R_S_n = quat2dcm(quat');
    R_S_n_1 = quat2dcm(quat1');
    acc_w = (R_S_n' * a + gw   +   R_S_n_1' * a1 + gw)/2;

    quat = F* quat;
    Pwb = Pwb + Vw * dt + 0.5 * dt * dt * acc_w;
    Vw = Vw + acc_w * dt;
%%
    quatAll(k,:) = quat';
    oula_new(k,:) = qtpoula(quat)';
    
    PwbSav(k,:) = Pwb';
end

%% 欧拉积分
% for k = 1:N
%     w =   imu_data(4:6, k);
%     a =   imu_data(1:3, k);
%     Ow = [0     -w(1)   -w(2)    -w(3);...
%           w(1)   0       w(3)    -w(2);...
%           w(2)  -w(3)    0        w(1);...
%           w(3)   w(2)   -w(1)     0  ];
%     Fc = 0.5*Ow;
%     F = eye(4) + Fc*dt;%%或者在这里改变积分的方法(对于矩阵的积分方法)
% 
%     %%
%     R_S_n = quat2dcm(quat');%%Rbw
%     
%     acc_w = R_S_n' * a + gw;
%     quat = F* quat;
%     quat = quat/norm(quat);
%     Vw = Vw + acc_w * dt;
%     Pwb = Pwb + Vw * dt + 0.5 * dt * dt * acc_w;
%     
% %%
%     quatAll(k,:) = quat';
%     oula_new(k,:) = qtpoula(quat)';
%     PwbSav(k,:) = Pwb';
% end


figure;plot(time , PwbSav-pos(:,1:3));title('中值积分的误差对比')
figure;plot(time, real(oula_new));hold on;
plot(time , pos(:,4:6))





%% fuction else
function Ang3 = qtpoula(q)
%四元数转欧拉角
x = atan2(2*(q(1)*q(2)+q(3)*q(4)),1 - 2*(q(2)^2+q(3)^2));
y = asin(2*(q(1)*q(3) - q(2)*q(4)));
z = atan2(2*(q(1)*q(4)+q(2)*q(3)),1 - 2*(q(3)^2+q(4)^2));
Ang3 = [x y z]';
end



% function getdata(t)

% end


function out = euler2Rotation(in) %%Rwb
x = in(1);
y = in(2);
z = in(3);
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
out = R;

end

function R  =  eulerRates2bodyRates( eulerAngles)

 roll = eulerAngles(1);
 pitch = eulerAngles(2);

 cr = cos(roll);
 sr = sin(roll);
 cp = cos(pitch);
 sp = sin(pitch);

R = [1, 0, -sp;
    0, cr, sr * cp;
    0, -sr, cr * cp];


end



