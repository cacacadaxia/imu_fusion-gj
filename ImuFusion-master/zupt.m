% =========================================================================
%
%                  IMU_kalman_filter
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 3��16��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.IMU��λdeg/s��[0,0,g]m/s^2,�Ƕȿ�����rad/s
%        2.ʹ�����ַ�����һ���ǿ����˲�������һ����IMU�򵥵���Ϣ�ں�
%        3.һ�ֻ���a zero velocity���˲��㷨���������ǹ۲��ٶ�Ϊ0m/s������IMU����
%           ���Ԥ��
%        4.�������� [https://github.com/xioTechnologies/Gait-Tracking-With-x-IMU]
%
%--------------------------------------------------------------------------

% This implements a zero velocity update EKF algorithm, which can be used
% for pedestrian tracking.
% By @TanTanTan
% ��ʵ�������ٶȸ���EKF�㷨�����������˸��١�
% Email: jinaqiaoqiao@gmail.com


clear; close all; 
% clc;

% --- Dataset parameters
dt = 1/256;                 %% �̶�ֵ
g = 9.81;
sigma_acc  = 0.1;           % accelerometer noise (m/s^2)
sigma_gyro = 0.1*pi/180;    % gyroscope noise (rad/s)
sigma_vel = 0.01;           % zero-velocity update measurement noise (m/s)


%% zero-velocity detector parameters
cova  = 0.01^2;
covw  = (0.1*pi/180)^2;
W     = 5;              % window size
gamma  = 0.3e5;

%% read data (Choose one)
% data = load('data/StaightLine.mat');
% data = load('data/StairsAndCorridor.mat');
data = load('data/SpiralStairs.mat');
imu_data = data.imu_data';
N = size(imu_data, 2);


%% Since we have got all the data, we can run the zero velocity detection
%  for the whole dataset.
%%���������Ѿ�������������ݣ���˿��Զ��������ݼ��������ٶȼ�⡣
iszv = zeros(1, N);
T = zeros(1, N-W+1);
for k = 1:N-W+1
    mean_a = mean(imu_data(1:3,k:k+W-1), 2);
    for l = k:k+W-1
        temp = imu_data(1:3,l) - g * mean_a / norm(mean_a);
        T(k) = T(k) + imu_data(4:6,l)'*imu_data(4:6,l)/covw + temp'*temp/cova;
    end
end
T = T./W;
for k = 1:size(T,2)
    if T(k) < gamma
        iszv(k:k+W-1) = ones(1,W);
    end
end


%%
%=========================================================================%
%==                             Initialization                           =%
%=========================================================================%
% We require that the IMU to be at static for first one second. So we can
% estimate the initial roll/pitch angle, as well as the biases.
%%����Ҫ��IMU�ڵ�һ���ӱ��־�̬����ˣ����ǿ��Թ��Ƴ�ʼ�Ĳ���/�������Լ�ƫ�
init_a = mean(imu_data(1:3,1:20),2);
init_a = init_a / norm(init_a);

init_psi =  0;
init_theta = -asin(init_a(1));
init_phi = atan2(init_a(2), init_a(3));

init_quat = angle2quat(init_psi, init_theta, init_phi); %%��ʼ�Ƕ�
%% test
% figure(5);plot(imu_data(1:3,1000:2000)');
% figure(6);plot(imu_data(4:6,1000:2000)');



% Estimate sensor bias. (���ֻ�Ƿǳ����Ե���ɱ궨����)
Rsw = quat2dcm(init_quat);
as  = Rsw * [0;0;g];
bias_a = mean(imu_data(1:3,1:500),2) - as;
bias_w = mean(imu_data(4:6,1:500),2);

% set the initial state vector
x = zeros(10,1);
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
    Q = Qc* dt^2;                 %%�趨Qֵ
    
    %% state propagation
    R_S_n = quat2dcm(quat');        %%guan to ben
    acc = R_S_n' * a - [0; 0;  g];  %%guan
%     Vq*quat - (acc+[0;0;g])
    x(1:3) = x(1:3) + x(4:6)* dt + 0.5*acc* dt^2;
    x(4:6) = x(4:6) + acc* dt;
    
    quat = (eye(4) + 0.5*Ow* dt)*quat;%%����Ĵ�������ٶȵĲ�һ��
    quat = quat/norm(quat);
    x(7:10) = quat;
    
    %% covariance propagation
    P = F*P*F' + Q;
    
    %% zero-velocity update
    if iszv(k) == 1
        K = (P*H')/(H*P*H'+R);
        y = -x(4:6);%%����û�������۲��ٶ�Ϊ0
        
        x = x + K*y;
        x(7:10) = x(7:10)/norm(x(7:10));
        
        P = (eye(10)-K*H)*P;
    end
    
    P = (P+P')/2;
    
    x_r(:,k) = x;
    x_oula(:,k) = qtpoula(quat);
    %% test for
    R_S_n = quat2dcm(quat'); %%
    acctp = R_S_n' * a - [0; 0;  g];%%guan
    xx(1:3) = xx(1:3) + xx(4:6)* dt + 0.5*acctp* dt^2;
    xx(4:6) = xx(4:6) + acctp* dt;
    
    
    Ksave(k,:,:) = K;
    xxsav = [xxsav , xx];
    
    %%
    Ksave2(:,k) = K*y;
    accsave(:,k) = acctp;
    quatsave(k,:) = quat';
end


%% View the results
% Plot the trajectory
        figure(1), hold on
        plot(x_r(1,:),x_r(2,:));
        plot(x_r(1,1),x_r(2,1),'ro');
        legend('Trajectory','Start point')
        axis equal
        grid on
        
        figure(2),
        plot(1:N,x_r(3,:));
        title('Height')
        grid on
    
%% �ǶȱȽ�
        % 
        % for i = 1:length(quatsave)
        %     oula = x_oula(:,i);
        % end
        % figure;plot(x_oula'/pi*180);
        % figure;plot(quatsave);

%% ֻ�ǶԽǶȽ��л��֣����ֻ��ַ�ʽ
init_quat_tp = init_quat;     %%��ʼ�Ƕȸ�ֵ
quat = init_quat_tp';
for k = 2:N
    w = imu_data(4:6, k-1);
    Ow = [0     -w(1)   -w(2)    -w(3);...
        w(1)   0       w(3)    -w(2);...
        w(2)  -w(3)    0        w(1);...
        w(3)   w(2)   -w(1)     0  ];
    Fc = 0.5*Ow;
    F = eye(4) + Fc*dt;%%����������ı���ֵķ���(���ھ���Ļ��ַ���)
    
    
    
    %%
    quat = F* quat;
    quat = quat/norm(quat);
    quatAll(k,:) = quat';
    oula_new(k,:) = qtpoula(quat)';
end
    figure;plot(quatsave);
    hold on;plot(quatAll);
    figure;plot(real(oula_new));hold on;
    plot(x_oula')

    
    figure;plot(oula_new-x_oula');


%% ���ýǶ�ֵ������ٶ�
accX = imu_data(1,:)';
accY = imu_data(2,:)';
accZ = imu_data(3,:)';
samplePeriod = dt;
% Compute accelerometer magnitude ������ٶȼƷ���
acc_mag = sqrt(accX.*accX + accY.*accY + accZ.*accZ);

% HP filter accelerometer data ��ͨ�˲���
filtCutOff = 0.001;  %%��ֹƵ���趨
[b, a] = butter(1, (2*filtCutOff)/(1/samplePeriod), 'high');
acc_magFilt = filtfilt(b, a, acc_mag);

% Compute absolute value
acc_magFilt = abs(acc_magFilt);

% LP filter accelerometer data
filtCutOff = 5;
[b, a] = butter(1, (2*filtCutOff)/(1/samplePeriod), 'low');
acc_magFilt = filtfilt(b, a, acc_magFilt);

% Threshold detection
stationary = acc_magFilt < 0.05*9.8;        %%��λ��һ��

%% �������˼���Ѿ�����˱Ƚ�׼ȷ����̬
accsave  = accsave';  %%3,n
% Integrate acceleration to yield velocity
vel = zeros(size(accsave));
for t = 2:length(vel)
    vel(t,:) = vel(t-1,:) + accsave(t,:) * 1/256;
    if(stationary(t) == 1)
        vel(t,:) = [0 0 0];     % force zero velocity when foot stationary
    end
end
% Compute integral drift during non-stationary periods
velDrift = zeros(size(vel));
stationaryStart = find([0; diff(stationary)] == -1);
stationaryEnd = find([0; diff(stationary)] == 1);
for i = 1:numel(stationaryEnd)
    driftRate = vel(stationaryEnd(i)-1, :) / (stationaryEnd(i) - stationaryStart(i));
    enum = 1:(stationaryEnd(i) - stationaryStart(i));
    drift = [enum'*driftRate(1) enum'*driftRate(2) enum'*driftRate(3)];
    velDrift(stationaryStart(i):stationaryEnd(i)-1, :) = drift;
end

% Remove integral drift
vel = vel - velDrift;
% Integrate velocity to yield position
pos = zeros(size(vel));
for t = 2:length(pos)
    pos(t,:) = pos(t-1,:) + vel(t,:) * 1/256;    % integrate velocity to yield position
end
time = 0:dt:dt*(length(pos)-1);
% Plot translational position
    % figure('Position', [9 39 900 600], 'NumberTitle', 'off', 'Name', 'Position');
    % hold on;
    % plot(time, pos(:,1), 'r');
    % plot(time, pos(:,2), 'g');
    % plot(time, pos(:,3), 'b');
    % title('Position');
    % xlabel('Time (s)');
    % ylabel('Position (m)');
    % legend('X', 'Y', 'Z');
    % hold off;

    figure;plot(pos(:,1),pos(:,2));grid on;axis equal;hold on;
    plot(x_r(1,:),x_r(2,:));plot(x_r(1,1),x_r(2,1),'ro');



%% another acc intel(һ�µľͲ�������)
        % vel_err = zeros(size(accsave));
        % for t = 2:length(vel_err)
        %     vel_err(t,:) = vel_err(t-1,:) + accsave(t,:) * 1/256;
        %     if(stationary(t) == 1)
        %         vel_err(t,:) = [0 0 0];     % force zero velocity when foot stationary
        %     end
        % end
        % vel_err = vel_err - velDrift;%%��һ������Ҫ��ֱ��Ӱ�쵽���
        % % Integrate velocity to yield position
        % pos_err = zeros(size(vel_err));
        % for t = 2:length(pos_err)
        %     pos_err(t,:) = pos_err(t-1,:) + vel_err(t,:) * 1/256;    % integrate velocity to yield position
        % end
        % 
        % figure;plot(pos_err);hold on;
        % plot(pos);
        % %%����Ľ�������һ�µ�

%% gyro�Ĵ��� (ֱ�ӶԽǶȽ��л��ִ����϶�ƫ��ܴ��)
%%x��ĳ�ֵƯ���Ѿ����˲�����
err = 0;
gyroX = imu_data(4,:) - err;
gyroY = imu_data(5,:);
gyroZ = imu_data(6,:);
gyro =  [gyroX ; gyroY ; gyroZ];
gyro = gyro';
angle = zeros(size(gyro));
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

%%ò�Ʋ�ıȽ϶�
figure;plot(angle(:,1));
hold on;
plot(x_oula(1,:)');
legend jifen kalman

%% Animation
%% The animation code was stolen from xioTechnologies.
%% https://github.com/xioTechnologies/Gait-Tracking-With-x-IMU
PLOT1 = 0
if PLOT1 ==1
    L = size(x_r,2);
    SamplePlotFreq = 4;
    Spin = 120;
    SixDofAnimation(x_r(1:3,:)', quat2dcm(quatinv(x_r(7:10,:)')), ...
        'SamplePlotFreq', SamplePlotFreq, 'Trail', 'All',...
        'Position', [9 39 1280 768],...
        'View', [(100:(Spin/(L-1)):(100+Spin))', 10*ones(L, 1)],...
        'AxisLength', 0.1, 'ShowArrowHead', false, ...
        'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)',...
        'ShowLegend', false,...
        'CreateVideo', false);
end


%% fuction else
function Ang3 = qtpoula(q)
%��Ԫ��תŷ����
x = atan2(2*(q(1)*q(2)+q(3)*q(4)),1 - 2*(q(2)^2+q(3)^2));
y = asin(2*(q(1)*q(3) - q(2)*q(4)));
z = atan2(2*(q(1)*q(4)+q(2)*q(3)),1 - 2*(q(3)^2+q(4)^2));
Ang3 = [x y z]';
end
