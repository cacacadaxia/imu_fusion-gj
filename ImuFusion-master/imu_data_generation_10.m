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
%       6. 改变运动的方式，变成限界的旋转
%       7.
%       8. 改变观测量，只有位移
%       9. 直接对四元数进行观测
%       10.使用真实的角速度行不行？不太好
%--------------------------------------------------------------------------



clear all;
close all;
dt = 1/256;
t_end = 20;
gw        = [0, 0, -9.81e3]';
rollpoint = [0,100,-1000]';
r = norm(rollpoint);
fai0 = -atan(rollpoint(2)/rollpoint(3));
time  = 0:dt:t_end;
K = pi/5;       %10s周期
A = 0.017;
fai = fai0+ A*sin(K*time);
dfai = A*K*cos(K*time);
ddfai = -A*K^2*sin(K*time);
k = 0;
for i = time
   k = k+1;
    angle = fai(k);
    oula = [angle;0;0];
    zeta = [0;0;r];%%本体系
    zeta0 = rollpoint;
    Rwc = euler2Rotation(oula);%%bw
    position = zeta0 + Rwc*zeta;
    
    
    %%
    v_c = [0; - dfai(k)*r;0];
    a_c = [0; - ddfai(k)*r ;0];
    v_w = Rwc*v_c;
    a_w = Rwc*a_c;
    imu_acc = Rwc'*(a_w - gw);
%     imu_acc = Rwc'*(a_w);
    imu_gyro = [dfai(k);0;0];
    imu_data(1:6,k) = [imu_acc;imu_gyro];
    %%
    eulerAngles = oula;
    pos(k,1:6) = [position ; eulerAngles];
    vel(k,1:6) = [v_w ; a_w];
    quatsaveTrue(k,:) = oula2quat(oula);
end

%% 加噪声
add_nosie = 1;
if add_nosie == 1
     sigma_acc  = 0.01;           % accelerometer noise (mm/s^2)
    sigma_acc_randwalk = 0.01;
    sigma_gyro = 0.01;    % gyroscope noise (rad/s)
    sigma_gyro_randwalk = 0.0002;
    sigma_vel = 0.01;
    sigma_pos = 1;
    sigma_angle_x = 2.8284e-04;
    sigma_quat1 = 2.8284e-06;
    sigma_quat2 = 2.8284e-04;
    sigma_quat3 = 1e-8;
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
    a =  imu_data(1:3, k-1);
    a1 = imu_data(1:3, k);
    Ow = [0     -w(1)   -w(2)    -w(3);...
          w(1)   0       w(3)    -w(2);...
          w(2)  -w(3)    0        w(1);...
          w(3)   w(2)   -w(1)     0  ];
    Fc = 0.5*Ow;
    F = eye(4) + Fc*dt;%%或者在这里改变积分的方法(对于矩阵的积分方法)
    quat1 = F* quat;
  
    %%
    Rbw = quat2dcm(quat');
    Rbw1 = quat2dcm(quat1');
    acc_w = (Rbw' * a + gw   +   Rbw1' * a1 + gw)/2;

    quat = F* quat;
    Pwb = Pwb + Vw * dt + 0.5 * dt * dt * acc_w;
    Vw = Vw + acc_w * dt;
%% data save
    quatAll(k,:) = quat';
    oula_new(k,:) = qtpoula(quat)';
    PwbSav(k,:) = Pwb';
end

figure;plot(time , PwbSav-pos(:,1:3));title('中值积分的误差对比')
figure;plot(pos(:,1:3));
%% EKF
    %  
    % % length of data
    % N = length(imu_data);
    % 
    % % set the initial state vector
    % x = zeros(5,1);
    % x(1:4) = [ pos(1,2:3)' ; vel(1,2:3)' ];
    % x(5,1) = init_ang(1);
    % 
    % % set the initial covariance
    % P = diag([1e-6*ones(1,5)]);
    % Q = diag([sigma_acc^2 , sigma_acc^2 , sigma_gyro^2]);
    % %
    % x_r = zeros(5,N);
    % x_r(:,1) = x;
    % 
    % % measurement matrix
    % pos_measure = pos(:,2:3) + sigma_pos*randn(size(pos(:,2:3)));
    % angle_measure = pos(:,4) + sigma_angle_x*randn(size(pos(:,3)));
    % H = [1,0,0,0,0;
    %     0,1,0,0,0;
    %     0,0,0,0,1];
    % R =  [sigma_pos^2.*eye(2),zeros(2,1);
    %     0,0,sigma_angle_x^2];
    % 
    % 
    % 
    % 
    % %%
    % xx = zeros(6,1);
    % xxsav = [];
    % %%
    % %=========================================================================%
    % %==                             Main  Loop    （简化之后的卡尔曼滤波）                           =%
    % %=========================================================================%
    % 
    % Ksave = [];
    % for k = 2:N
    % 
    %     %% compute state transition matrix F and covariance Q
    %     w = imu_data(4, k-1); % - bias_w;
    %     a = imu_data(1:2, k-1); % - bias_a;
    %     
    %     fai = x(5);
    %     Rtp = [cos(fai) , -sin(fai) , 0;
    %         sin(fai) , cos(fai) , 0;
    %         0,0,1];
    %     
    %     F = eye(5);
    %     F(1,3) = dt;
    %     F(2,4) = dt;
    %     
    %     Bu = [zeros(2,3) ; dt*eye(3)];
    %     
    %     %%
    %     Rwb = [cos(fai) , -sin(fai);
    %         sin(fai) , cos(fai)];
    %     acc = Rwb*a + [0;0];
    %     x(1:2) = x(1:2) + x(3:4)* dt + 0.5*acc* dt^2;
    %     x(3:4) = x(3:4) + acc* dt;
    %     x(5) = x(5) + w;
    % 
    %     P = F*P*F' + Bu*Rtp*Q*Rtp'*Bu';
    %     %%
    %     K = (P*H')/(H*P*H' + R);
    % 
    %     mea = [pos_measure(k,:)' ; angle_measure(k) ];
    %     y = mea - H*x; 
    %     x = x + K*y;
    %     
    %     P = (eye(5)-K*H)*P;
    %     P = (P + P')/2;
    %     x_r(:,k) = x;
    % end
    % figure;plot(x_r(1,:));hold on;plot(pos(:,2))
    % figure;plot(x_r(2,:));hold on;plot(pos(:,3))
    % figure;plot(x_r(5,:));hold on;plot(pos(:,4))
%% EKF
 
% length of data
N = length(imu_data);

% set the initial state vector
x = zeros(6,1);
x(1:6) = [ pos(1,1:3)' ; vel(1,1:3)' ];
% x(7:10,1) = init_quat';

% set the initial covariance
P = diag([1e-10*ones(1,6)]);
%
x_r = zeros(6,N);
x_r(:,1) = x;

% measurement matrix
vel_measure = vel(:,1:3) + sigma_vel*randn(size(vel(:,1:3)));
pos_measure = pos(:,1:3) + sigma_pos*randn(size(pos(:,1:3)));

quat_measure = quatsaveTrue + randn(size(quatsaveTrue)).*[sigma_quat1,sigma_quat2,0,0];
H = [eye(3)    , zeros(3,3)];
R =  [ sigma_pos^2 * eye(3)];

%%
%=========================================================================%
%==                             Main  Loop                               =%
%=========================================================================%
Ksave = [];
for k = 2:N
    %% compute state transition matrix F and covariance Q
%     w = imu_data(4:6, k-1); % - bias_w;
%     quat = x(7:10,1);
    quat = quat_measure(k,:);
    a = imu_data(1:3, k-1); % - bias_a;
    
    % continuous state transition matrix
    %%连续状态转移矩阵
    Ow = [0     -w(1)   -w(2)    -w(3);...
          w(1)   0       w(3)    -w(2);...
          w(2)  -w(3)    0        w(1);...
          w(3)   w(2)   -w(1)     0  ];   %%q的变化速度
    Vq = compVq(quat, a);%%3,4 加速度变形之后，对于位置的速度
    
    Fc = zeros(10);
    Fc(1:3, 4:6) = eye(3);
    Fc(4:6,7:10)= [Vq; 0.5*Ow];%%10,10 连续Fc
    
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
    Q = Qc * dt;                 %%设定Q值
    
    %% state propagation
    R_S_n = quat2dcm(quat');        %%guan to ben
    RbwTrue = quat2dcm(quatsaveTrue(k,:));%%没用
    acc = R_S_n' * a +gw;  %%guan
%     Vq*quat/2 - acc + gw
    x(1:3) = x(1:3) + x(4:6)* dt + 0.5*acc* dt^2;
    x(4:6) = x(4:6) + acc* dt;
    
    quat = (eye(4) + 0.5*Ow* dt)*quat;%%这里的处理与加速度的不一样
    quat = quat/norm(quat);
    x(7:10) = quat;
    
    %% covariance propagation
    P = F*P*F' + Q;

    %% zero-velocity update
    K = (P*H')/(H*P*H'+R);
    mea = [pos_measure(k,:) , quat_measure(k,:)];
    y = mea'- H*x; 

    x = x + K*y;
    x(7:10) = x(7:10)/norm(x(7:10));

    P = (eye(10)-K*H)*P;
    P = (P+P')/2;
    x_r(:,k) = x;
    x_oula(:,k) = qtpoula(quat);
    quatsave(:,k) = quat;
end
%%position
figure;plot3(x_r(1,:),x_r(2,:),x_r(3,:));
hold on;plot3(pos(:,1),pos(:,2),pos(:,3));
%%err
figure;plot(time,pos(:,1:3)- x_r(1:3,:)');

%%angle
figure;plot(x_oula'- pos(:,4:6));


std(pos(:,1:3)- x_r(1:3,:)')
std(x_oula'- pos(:,4:6))/pi*180

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
    sin(z) cos(z)    0;
    0      0         1];
R = Rz*Ry*Rx;
out = R;
end

function Ang3 = qtpoula(q)
%四元数转欧拉角
x = atan2(2*(q(1)*q(2)+q(3)*q(4)),1 - 2*(q(2)^2+q(3)^2));
y = asin(2*(q(1)*q(3) - q(2)*q(4)));
z = atan2(2*(q(1)*q(4)+q(2)*q(3)),1 - 2*(q(3)^2+q(4)^2));
Ang3 = [x y z]';
end

function Q2 = oula2quat(in)
x = in(1);
y = in(2);
z = in(3);
%欧拉角转四元数
q = [cos(x/2)*cos(y/2)*cos(z/2) + sin(x/2)*sin(y/2)*sin(z/2) ...
    sin(x/2)*cos(y/2)*cos(z/2) - cos(x/2)*sin(y/2)*sin(z/2) ...
    cos(x/2)*sin(y/2)*cos(z/2) + sin(x/2)*cos(y/2)*sin(z/2) ...
    cos(x/2)*cos(y/2)*sin(z/2) - sin(x/2)*sin(y/2)*cos(z/2)];
Q2 = q;

end
