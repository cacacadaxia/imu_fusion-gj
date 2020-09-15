
% =========================================================================
%
%                  �������IMU
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 3��14��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.����ipython notebook�����ݣ���������ʹ���������ݼ����е�
%        2. ��λ��deg/s [0,0,-g]m/s^2
%        3.
%        4.
%
%--------------------------------------------------------------------------
close all;
clear all;
imu_data = dlmread('test.txt');
imu_data = imu_data';
N = size(imu_data, 2);
figure;plot(imu_data(1:3,:)');
figure;plot(imu_data(4:6,:)');














%% View the results
%% Plot the trajectory
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