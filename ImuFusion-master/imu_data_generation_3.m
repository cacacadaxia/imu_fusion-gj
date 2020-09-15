% =========================================================================
%
%                  IMU_kalman_filter
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 4��14��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.����һ������ [m/s^2] �� [rad/s]
%        2.������֮��
% in_accel:
%    noise����: [ 0.016925432397973516, 0.016735310195561025, 0.017452487504590969 ]
%    bias ƫ��: [ 0.00019031356589714596, 0.00016996777864587261, 0.00054490537096493644 ]
% in_gyro:
%    noise: [ 0.0010848026158819934, 0.0012466367883501759, 0.0011003229919806443 ]
%    bias: [ 0.000023404834136742844, 0.000023596771567764949, 0.000014970418056326829 ]
% 3. EKF
%--------------------------------------------------------------------------


clear all;
close all;
end_t = 80;
[imu_data ,  pos , vel , time , dt]= data_gener(end_t);
PLOT = 0;
if PLOT
    figure;plot(pos(:,1:3));title('λ��')
    figure;plot(pos(:,4:6));title('�Ƕ�')
end

%% ������

sigma_acc  = 0.001;           % accelerometer noise (m/s^2)
sigma_gyro = 0.01*pi/180;    % gyroscope noise (rad/s)
sigma_vel = 0.01;

%%complex add noise :randn(2,4).*[100;100000]
imu_data(1:3,:) = imu_data(1:3,:) + randn(size(imu_data(1:3,:)))*sigma_acc;
imu_data(4:6,:) = imu_data(4:6,:) + randn(size(imu_data(4:6,:)))*sigma_gyro;
%% ��ֵ����
gw        = [0, 0, -9.81]';
init_ang  = pos(1,4:6);  %%��ʼŷ����
init_quat = angle2quat(init_ang(3),init_ang(2),init_ang(1));%%��ʼ��Ԫ��
quat      = init_quat';     %%��ʼ�Ƕȸ�ֵ
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
    F = eye(4) + Fc*dt;%%����������ı���ֵķ���(���ھ���Ļ��ַ���)
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

%%
% �۲�
vel_measure = vel(:,1:3) + sigma_vel*randn(size(vel(:,1:3)));
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
%%����۲��ǲ��ǲ�̫�׵�
H = [zeros(3), eye(3), zeros(3,4)];
R =  sigma_vel^2 * eye(3);


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
    %%����״̬ת�ƾ���
    Ow = [0     -w(1)   -w(2)    -w(3);...
          w(1)   0       w(3)    -w(2);...
          w(2)  -w(3)    0        w(1);...
          w(3)   w(2)   -w(1)     0  ];   %%q�ı仯�ٶ�
    Vq = compVq(quat, a);%%3,4 ���ٶȱ���֮�󣬶���λ�õ��ٶ�
    
    Fc = zeros(10);
    Fc(1:3, 4:6) = eye(3);
    Fc(4:10,7:10)= [Vq; 0.5*Ow];%%10,10 ����Fc
    
    % continuous process covariance
    Gq = 0.5* [-quat(2)  -quat(3)   -quat(4); ...
                quat(1)  -quat(4)    quat(3); ...
                quat(4)   quat(1)   -quat(2); ...
               -quat(3)   quat(2)    quat(1)];       %%����״ֵ̬����Gq
    Qc = zeros(10);                             %%���ϸ���Q
    Qc(4:6, 4:6)  =  sigma_acc^2*eye(3);        %%����ٶ���ص�ֵ(����ʼ��������ڹ�������������)
    Qc(7:10,7:10) =  sigma_gyro^2*(Gq*Gq');     %%���Qֵ�ķ���
    %%Gq*Gq' 4��4
    % discretilization
    F = eye(10) + Fc* dt;       %%���ùܣ�
    Q = Qc* dt;                 %%�趨Qֵ
    
    %% state propagation
    R_S_n = quat2dcm(quat');        %%guan to ben
    acc = R_S_n' * a + gw;  %%guan
    x(1:3) = x(1:3) + x(4:6)* dt + 0.5*acc* dt^2;
    x(4:6) = x(4:6) + acc* dt;
    
    quat = (eye(4) + 0.5*Ow* dt)*quat;%%����Ĵ�������ٶȵĲ�һ��
    quat = quat/norm(quat);
    x(7:10) = quat;
    
    %% covariance propagation
    P = F*P*F' + Q;
    
    %% zero-velocity update
    K = (P*H')/(H*P*H'+R);
    mea = vel_measure(k,:);
    y = mea'-x(4:6);        %%����û�������۲��ٶ�Ϊ0

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
figure;plot(pos(:,1:3)- x_r(1:3,:)');
hold on;plot(PwbSav-pos(:,1:3));


if PLOT
    figure;plot(time , PwbSav-pos(:,1:3));title('��ֵ���ֵ����Ա�')
%     figure;plot(time, real(oula_new));hold on;
%     plot(time , pos(:,4:6))
%     figure;plot(pos(:,4:6)-real(oula_new))
end
%% fuction else
function Ang3 = qtpoula(q)
%��Ԫ��תŷ����
x = atan2(2*(q(1)*q(2)+q(3)*q(4)),1 - 2*(q(2)^2+q(3)^2));
y = asin(2*(q(1)*q(3) - q(2)*q(4)));
z = atan2(2*(q(1)*q(4)+q(2)*q(3)),1 - 2*(q(3)^2+q(4)^2));
Ang3 = [x y z]';
end

function out = euler2Rotation(in) %%Rwb
x = in(1);
y = in(2);
z = in(3);
%ŷ����ת��ת����
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



