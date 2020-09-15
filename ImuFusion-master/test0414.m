clear;close all;
% ----------------------------------------------------------------------
%  step 1 �����ź�
% ------------------------------------------------------------
sr = 1000;                       % �����ʣ�1s�ڷ��ŵ�ĸ���(������ʱ��)
t = 1:1/sr:10;                  % ��Ժ�����ʱ��
f_sig = [0:10:100 ];           % ʵ��
sig = 0;
sigmode = 'real';
switch sigmode
    case 'real'
        for i = 1:length(f_sig)
            sig = sig + sin(2*pi*f_sig(i)*t)*(1+1j);          % ģ��OFDM�ź�
        end
    case 'imag'
        for i = 1:length(f_sig)
            sig = sig + exp(1j*2*pi*f_sig(i)*t);
        end
end

f_B = max(f_sig)-min(f_sig);                        % Ƶ�ʣ�ͬʱҲ�Ǵ���B
% figure;plot(t,sig);hold on
% Ϊʲô�����ʺ��źŴ���û��ϵ����Ϊ����sin�ź�
% ���������ͨ�ĵ����źţ�QPSK����BPSK������ô�źŵĴ���B����sr/2�����Բ����Ļ���ֻ��Ҫ2*B=sr�Ϳ��Բ�����
% Ϊʲô����Ϊʵ�ź��븴�ź�֮�������
% ------------------------------------------------
% step 2 ���Բ�ֵ+����
% ------------------------------------------------
fs   = 180;
% fs   = 500;                     % ����Ƶ��
t_cy = 1:1/fs:10;               % ������ʱ��
sig_fj = interp1(t,sig,t_cy,'spline');      % �����ĺ���
% -----------------------------------------------
% fft ʵ���źŵ�Ƶ��
NN = length(sig);figure;plot((-NN/2+1:NN/2)/NN*sr,20*log10(abs(fftshift(fft(sig)))));axis([-inf,inf,-inf,inf]);xlabel('fs Hz');ylabel('������ dB')
% fft �������Ƶ��
NN = length(sig_fj);figure;plot((-NN/2+1:NN/2)/NN*fs,20*log10(abs(fftshift(fft(sig_fj)))));title('������');axis([-inf,inf,-inf,inf]);xlabel('fs Hz');ylabel('������ dB')

%% ��ԭʼ�ź����������ź��úܴ�Ĳ��������²������۲����ս��
fmax = 2000;
tmax = 1:1/fmax:10;
sig1 = interp1(t,sig,tmax,'spline');
sig2 = interp1(t_cy,sig_fj,tmax,'spline');


figure;plot(real(sig1));
figure;plot(real(sig2-sig1));
sig_err = sig2-sig1;
NN = length(sig_err);figure;plot((-NN/2+1:NN/2)/NN*fmax,20*log10(abs(fftshift(fft(sig_err)))));title('������');axis([-inf,inf,-inf,inf]);xlabel('fs Hz');ylabel('������ dB')

% ------------------------------------------------
% Ƶ�װ���