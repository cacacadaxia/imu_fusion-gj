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
% add at 0915
%       5. 讨论位置的积分问题，对比不同的积分
%       6. 
%--------------------------------------------------------------------------


clear all;
close all;
end_t = 80;
[imu_data ,  pos , vel , time , dt] = data_gener(end_t);
PLOT = 1;
if PLOT
%     figure;plot(pos(:,1:3));title('位置')
%     figure;plot(pos(:,4:6)/pi*180);title('角度')
end

%% 加噪声
add_nosie = 1;
if add_nosie == 1
     sigma_acc  = 0.0016;           % accelerometer noise (m/s^2)
    sigma_acc_randwalk = 0.0001;
    sigma_gyro = 0.01*pi/180;    % gyroscope noise (rad/s)
    sigma_gyro_randwalk = 0.00002;
    sigma_vel = 0.01;
    sigma_pos = 0.0001;
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

%% 中值积分 (第一种积分方法)（常规）
gw        = [0, 0, 9.81]';
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
    R_S_n = quat2dcm(quat');
    R_S_n_1 = quat2dcm(quat1');
    acc_w = (R_S_n' * a + gw   +   R_S_n_1' * a1 + gw)/2;

    quat = F* quat;
    Pwb = Pwb + Vw * dt + 0.5 * dt * dt * acc_w;
    Vw = Vw + acc_w * dt;
%% data save
    quatAll(k,:) = quat';
    oula_new(k,:) = qtpoula(quat)';
    PwbSav(k,:) = Pwb';
end
figure;plot(PwbSav - pos(:,1:3))

PwbSav1 = PwbSav;


%% 中值积分 (第二种积分方法)
gw        = [0, 0, 9.81]';
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
%     R_b_w = quat2dcm(quat');
%     R_b_w_1 = quat2dcm(quat1');
%     acc_w = (R_b_w' * a + gw   +   R_b_w_1' * a1 + gw)/2;
% 
%     quat = F* quat;
%     Pwb = Pwb + Vw * dt + 0.5 * dt * dt * acc_w;
%     Vw = Vw + acc_w * dt;
    
    %%
    quat = F* quat;
    R_b_w = quat2dcm(quat');
    B_a = a + R_b_w*gw;     %%去除g之后的在b系下的加速度
    W_a = R_b_w'*B_a;
    tmp = W_a*dt*dt;    %%二阶差分
    
    
    % 下面是积分环节
    % 这个应该就是二阶差分的意思吧，实际上和之前也没什么区别
    Vw = Vw + W_a*dt;
    Pwb = Pwb + Vw*dt;
    
%% data save 
    quatAll(k,:) = quat';
    oula_new(k,:) = qtpoula(quat)';
    PwbSav(k,:) = Pwb';
end
figure;plot(PwbSav - pos(:,1:3))
hold on;
plot(PwbSav1- pos(:,1:3))


%% 中值积分 (第三种积分方法)
gw        = [0, 0, 9.81]';
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
%     R_b_w = quat2dcm(quat');
%     R_b_w_1 = quat2dcm(quat1');
%     acc_w = (R_b_w' * a + gw   +   R_b_w_1' * a1 + gw)/2;
% 
%     quat = F* quat;
%     Pwb = Pwb + Vw * dt + 0.5 * dt * dt * acc_w;
%     Vw = Vw + acc_w * dt;
    
    %%
    quat = F* quat;
    R_b_w = quat2dcm(quat');
    B_a = a + R_b_w*gw;     %%去除g之后的在b系下的加速度
    W_a = R_b_w'*B_a;
    tmp = W_a*dt*dt;    %%二阶差分
    %% 这个积分方法始终不对
    
    % 下面是积分环节
    % 这个应该就是二阶差分的意思吧，实际上和之前也没什么区别
    Vw = Vw + tmp;
    Pwb = Pwb + Vw;
    
    
%% data save 
    quatAll(k,:) = quat';
    oula_new(k,:) = qtpoula(quat)';
    PwbSav(k,:) = Pwb';
end
figure;plot(PwbSav - pos(:,1:3))



%%

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
    w = imu_data(4:6, k-1); % - bias_w;
    quat = x(7:10,1);
    a = imu_data(1:3, k-1); % - bias_a;
    
    % continuous state transition matrix
    %%连续状态转移矩阵
    Ow = [0     -w(1)   -w(2)    -w(3);...
          w(1)   0       w(3)    -w(2);...
          w(2)  -w(3)    0        w(1);...
          w(3)   w(2)   -w(1)     0  ];   %%q的变化速度
    Vq = compVq(quat, a)/3;%%3,4 加速度变形之后，对于位置的速度
    
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
    
    %% state propagation
    R_S_n = quat2dcm(quat');        %%guan to ben
    acc = R_S_n' * a + gw;  %%guan
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
% figure;plot3(x_r(1,:),x_r(2,:),x_r(3,:));
% hold on;plot3(pos(:,1),pos(:,2),pos(:,3));
%%err
figure;plot(pos(:,1:3) -  x_r(1:3,:)');
figure;plot(pos(:,4:6) -  x_oula');title('EKF角度误差对比')
% hold on;plot(PwbSav-pos(:,1:3));

PLOT = 0;
if PLOT
    figure;plot(time , PwbSav - pos(:,1:3));title('中值积分的位置误差对比');
%     figure;plot(time, real(oula_new));hold on;
    figure;plot(pos(:,4:6) - real(oula_new));title('中值积分的角度误差对比');
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



