%%  ���Ȱ�ռ� ˮƽ��żԴ�ڵر�����Ӧ����
clc;clear all;close all;
%% ��������
u0 = 4*pi*1e-7;
fs = 1e7;
dt =1./fs;
Ns = 2e4;
Tob = (1:Ns)./fs;
rou = 100;% ������
I = 1;
ds = 1;
m = I.*ds;
xr = 100;
yr = 0;
zr = 0;
%%
[step_ex10,step_ey10,impulse_hx,impulse_hy,impulse_hz] = Horizontal_Electrical_Dipole_Response(u0,rou,m,xr,yr,zr,Tob);
%% ˮƽ�糡�Ľ�Ծ��Ӧ
figure;
plot(Tob,step_ex10,'r','linewidth',2);
hold on;
plot(Tob,step_ey10,'b','linewidth',2);
legend('ex','ey');
title('ˮƽ��ż��Դ��(0,100,0)�Ľ�Ծˮƽ�糡');
xlabel('Time/s');
ylabel('E/(V/m)');
%% �ų���������Ӧ hx
figure;
plot(Tob,u0.*impulse_hx,'r','linewidth',2);
% hold on;
% plot(Tob,impulse_hy,'b-.','linewidth',2);
% hold on;
% plot(Tob,impulse_hz,'c','linewidth',2);
% legend('hx','hy','hz');
title('ˮƽ��ż��Դ��(0,100,0)�����ų���������ӦBx');
xlabel('Time/s');
ylabel('B/(T)');
% �ų��Ľ�Ծ��Ӧ hx
%%  �ų���������Ӧ hy
figure;
% plot(Tob,impluse_hx,'r--','linewidth',2);
% hold on;
plot(Tob,u0.*impulse_hy,'b','linewidth',2);
% hold on;
% plot(Tob,impulse_hz,'c','linewidth',2);
% legend('hx','hy','hz');
grid on;
title('ˮƽ��ż��Դ��(0,100,0)�����ų���������ӦBy');
xlabel('Time/s');
ylabel('B/(T)');
%% �ų���������Ӧ hz
figure;
% plot(Tob,impluse_hx,'r--','linewidth',2);
% hold on;
% plot(Tob,impulse_hy,'b-.','linewidth',2);
% hold on;
plot(Tob,impulse_hz.*u0,'r','linewidth',2);
% legend('hx','hy','hz');
title('ˮƽ��ż��Դ��(0,100,0)�����ų���������ӦBz');
xlabel('Time/s');
ylabel('B/(T)');

%% 