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
%  功能： 1.生成一组数据 [m/s^2] 与 [rad/s]
%        2.加噪声之后
% in_accel:
%    noise噪声: [ 0.016925432397973516, 0.016735310195561025, 0.017452487504590969 ]
%    bias 偏差: [ 0.00019031356589714596, 0.00016996777864587261, 0.00054490537096493644 ]
% in_gyro:
%    noise: [ 0.0010848026158819934, 0.0012466367883501759, 0.0011003229919806443 ]
%    bias: [ 0.000023404834136742844, 0.000023596771567764949, 0.000014970418056326829 ]
%       3. EKF
%       4. 增加位置的观测量
%       5. 减少数据倍数，改变噪声的刻度
%       5.2 中值积分+EKF+LK
%--------------------------------------------------------------------------


clear all;
close all;
end_t = 80;
dt = 1/256;
M_PI = 3.1415926;
gn = [0, 0, -9.81e3];
k = 0;
time = 0:dt:end_t;
for t = 0:dt:end_t
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
PLOT = 1;
if PLOT
    figure;plot(pos(:,1:3));title('位置')
    figure;plot(pos(:,4:6)/pi*180);title('角度')
end

%% 加噪声
add_nosie = 1;
if add_nosie == 1
     sigma_acc  = 1;           % accelerometer noise (mm/s^2)
    sigma_acc_randwalk = 0.1;
    sigma_gyro = 0.01*pi/180;    % gyroscope noise (rad/s)
    sigma_gyro_randwalk = 0.0002;
    sigma_vel = 0.01;
    sigma_pos = 0.3;
else
    sigma_acc  = 0.0;           % accelerometer noise (m/s^2)
    sigma_acc_randwalk = 0.000;
    sigma_gyro = 0.0*pi/180;    % gyroscope noise (rad/s)
    sigma_gyro_randwalk = 0.0000;
    sigma_vel = 0.0;
    sigma_pos = 0.00;
end
%%complex add noise :randn(2,4).*[100;100000]
imu_data(1:3,:) = imu_data(1:3,:) + randn(size(imu_data(1:3,:)))*sigma_acc;
imu_data(4:6,:) = imu_data(4:6,:) + randn(size(imu_data(4:6,:)))*sigma_gyro;

walkbias_acc = 0;
walkbias_gyro = 0;
for i = 1:length(imu_data)
    walkbias_acc = walkbias_acc + sigma_acc_randwalk * randn(3,1);
    imu_data(1:3,i) = walkbias_acc + imu_data(1:3,i);
    walkbias_gyro = walkbias_gyro + sigma_gyro_randwalk * randn(3,1);
    imu_data(4:6,i) = walkbias_gyro + imu_data(4:6,i);
end

%% 中值积分
gw        = [0, 0, -9.81e3]';
init_ang  = pos(1,4:6);  %%初始欧拉角
init_quat = angle2quat(init_ang(3),init_ang(2),init_ang(1));%%初始四元数
quat      = init_quat';     %%初始角度赋值
N         = length(imu_data);
Vw        = vel(1,1:3)';
Pwb       = pos(1,1:3)';
for k = 2:N
    w = (imu_data(4:6, k-1)    +    imu_data(4:6, k))/2;
    ap =  imu_data(1:3, k-1);
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
    acc_w = (R_S_n' * ap + gw   +   R_S_n_1' * a1 + gw)/2;

    quat = F* quat;
    Pwb = Pwb + Vw * dt + 0.5 * dt * dt * acc_w;
    Vw = Vw + acc_w * dt;
%% data save
    quatAll(k,:) = quat';
    oula_new(k,:) = qtpoula(quat)';
    PwbSav(k,:) = Pwb';
end

%% EKF
 
% length of data
N = length(imu_data);
g = 9.81;

% set the initial state vector
x = zeros(10,1);
x(1:6) = [ pos(1,1:3)' ; vel(1,1:3)' ];
x(7:10,1) = init_quat';

% set the initial covariance
P = diag([1e-10*ones(1,6), 1e-6*ones(1,4)]);
%
x_r = zeros(10,N);
x_r(:,1) = x;

% measurement matrix
vel_measure = vel(:,1:3) + sigma_vel*randn(size(vel(:,1:3)));
pos_measure = pos(:,1:3) + sigma_pos*randn(size(pos(:,1:3)));
H = [eye(3),zeros(3,7);
    zeros(3) , eye(3), zeros(3,4)];
R =  [sigma_pos^2 * eye(3)  , zeros(3,3);
    zeros(3,3) , sigma_vel^2 * eye(3)];


%%
xx = zeros(6,1);
xxsav = [];
%%
%=========================================================================%
%==                             Main  Loop                               =%
%=========================================================================%
Ksave = [];
for k = 2:N
    %% compute state transition matrix F and covariance Q
    wp = imu_data(4:6, k-1); % - bias_w;
    w1 = imu_data(4:6, k);
    w = (wp + w1)/2;
    quat = x(7:10,1);
    ap = imu_data(1:3, k-1); % - bias_a;
    a1 = imu_data(1:3, k);
    % continuous state transition matrix
    %%连续状态转移矩阵
    Ow = [0     -w(1)   -w(2)    -w(3);...
          w(1)   0       w(3)    -w(2);...
          w(2)  -w(3)    0        w(1);...
          w(3)   w(2)   -w(1)     0  ];   %%q的变化速度
    Vq = compVq(quat, ap);%%3,4 加速度变形之后，对于位置的速度
    
    Fc = zeros(10);
    Fc(1:3, 4:6) = eye(3);
    Fc(4:10,7:10)= [Vq; 0.5*Ow];%%10,10 连续Fc
    
    % continuous process covariance
    Gq = 0.5* [-quat(2)  -quat(3)   -quat(4); ...
                quat(1)  -quat(4)    quat(3); ...
                quat(4)   quat(1)   -quat(2); ...
               -quat(3)   quat(2)    quat(1)];       %%利用状态值更新Gq
    Qc = zeros(10);                             %%不断更新Q
    Qc(4:6, 4:6)  =  sigma_acc^2*eye(3);        %%与加速度相关的值(这里始终是相对于过程噪声来讲的)
    Qc(7:10,7:10) =  sigma_gyro^2*(Gq*Gq');     %%获得Q值的方法
    %%Gq*Gq' 4，4
    % discretilization
    F = eye(10) + Fc* dt;       %%不用管？
    Q = Qc* dt;                 %%设定Q值
    
    %% 下一个q
%     quat1 = (eye(4) + 0.5*Ow* dt)*quat;%%这里的处理与加速度的不一样
    quat1 = attitude_update_RK4(quat,dt,wp,w1);
    quat1 = quat1/norm(quat1);
    R_S_n1 = quat2dcm(quat1');
    
    %% state propagation
    R_S_n = quat2dcm(quat');        %%guan to ben
    acc = (R_S_n' * ap + gw    + R_S_n1' * a1 + gw)/2;  %%guan
%     Vq*quat/2 - acc + gw
    x(1:3) = x(1:3) + x(4:6)* dt + 0.5*acc* dt^2;
    x(4:6) = x(4:6) + acc* dt;
    
%     quat = (eye(4) + 0.5*Ow* dt)*quat;%%这里的处理与加速度的不一样
    quat = attitude_update_RK4(quat,dt,wp,w1);
    quat = quat/norm(quat);
    x(7:10) = quat;
    
    %% covariance propagation
    P = F*P*F' + Q;
    
    %% zero-velocity update
    K = (P*H')/(H*P*H'+R);
    mea = [pos_measure(k,:) , vel_measure(k,:)];
    y = mea'-x(1:6);        %%这里没看懂，观测速度为0

    x = x + K*y;
    x(7:10) = x(7:10)/norm(x(7:10));

    P = (eye(10)-K*H)*P;
    P = (P+P')/2;
    x_r(:,k) = x;
    x_oula(:,k) = qtpoula(quat);
end
%%position
figure;plot3(x_r(1,:),x_r(2,:),x_r(3,:));
hold on;plot3(pos(:,1),pos(:,2),pos(:,3));
%%err
figure;plot(time,pos(:,1:3)- x_r(1:3,:)');
% hold on;plot(PwbSav-pos(:,1:3));



%% plot
if PLOT
%     figure;plot(time , PwbSav-pos(:,1:3));title('中值积分的误差对比')
%     figure;plot(time, real(oula_new));hold on;
%     plot(time , pos(:,4:6))
%     figure;plot(pos(:,4:6)-real(oula_new))
end
%% fuction else
function Ang3 = qtpoula(q)
%四元数转欧拉角
x = atan2(2*(q(1)*q(2)+q(3)*q(4)),1 - 2*(q(2)^2+q(3)^2));
y = asin(2*(q(1)*q(3) - q(2)*q(4)));
z = atan2(2*(q(1)*q(4)+q(2)*q(3)),1 - 2*(q(3)^2+q(4)^2));
Ang3 = [x y z]';
end

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

function [Qk_plus1] = attitude_update_RK4(Qk,dt,gyro0,gyro1)
    % RK4
    % conference: A Robust and Easy to implement method for imu
    % calibration without External Equipments

    q_1=Qk;
    k1=(1/2)*omegaMatrix(gyro0)*q_1;
    q_2=Qk+dt*(1/2)*k1;
    k2=(1/2)*omegaMatrix((1/2)*(gyro0+gyro1))*q_2;
    q_3=Qk+dt*(1/2)*k2;
    k3=(1/2)*omegaMatrix((1/2)*(gyro0+gyro1))*q_3;
    q_4=Qk+dt*k3;
    k4=(1/2)*omegaMatrix(gyro1)*q_4;
    Qk_plus1=Qk+dt*(k1/6+k2/3+k3/3+k4/6);
    Qk_plus1=Qk_plus1/norm(Qk_plus1);
end
function [omega]=omegaMatrix(data)

    % wx=data(1)*pi/180;
    % wy=data(2)*pi/180;
    % wz=data(3)*pi/180;
    wx=data(1);
    wy=data(2);
    wz=data(3);

    omega=[0  , -wx , -wy , -wz ;...
        wx ,  0  ,  wz , -wy ;...
        wy , -wz ,  0  ,  wx ;...
        wz ,  wy , -wx ,  0   ];

end


