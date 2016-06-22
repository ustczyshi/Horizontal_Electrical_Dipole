%%  ���Ȱ�ռ� ˮƽ��żԴ�ڵر�����Ӧ����
clc;clear all;close all;
%% ��������
u0 = 4*pi*1e-7;
fs = 1e5;
dt =1./fs;
Ns = 4e3;
Tob = (1:Ns)./fs;
rou = 100;% ������
I = 1;
ds = 1;
m = I.*ds;
xr = 100;
yr = 100;
zr = 0;
position = ['(' num2str(xr) ',' num2str(yr) ',' num2str(zr) ')' ];
%%
[step_ex01,step_ey01,impulse_hx,impulse_hy,impulse_hz] = Horizontal_Electrical_Dipole_Response(u0,rou,m,xr,yr,zr,Tob);

save('horizontal_electrical_dipole_impulse_jiexijie','step_ex01','step_ey01','impulse_hx','impulse_hy','impulse_hz');
%% ˮƽ�糡�Ľ�Ծ��Ӧ
%{
figure;
plot(Tob,step_ex01,'r','linewidth',2);
hold on;
plot(Tob,step_ey01,'b','linewidth',2);
legend('ex','ey');
title(['ˮƽ��ż��Դ��' position '�ĸ���Ծˮƽ�糡']);
xlabel('Time/s');
ylabel('E/(V/m)');
%% �ų���������Ӧ hx
figure;
plot(Tob,u0.*(impulse_hx),'r','linewidth',2);
% hold on;
% plot(Tob,impulse_hy,'b-.','linewidth',2);
% hold on;
% plot(Tob,impulse_hz,'c','linewidth',2);
% legend('hx','hy','hz');
grid on;
title([ 'ˮƽ��ż��Դ��' position '�����ų���������ӦBx']);
xlabel('Time/s');
ylabel('B/(T)');

%%  �ų���������Ӧ hy
figure;
% plot(Tob,impluse_hx,'r--','linewidth',2);
% hold on;
plot(Tob,u0.*impulse_hy,'b','linewidth',2);
% hold on;
% plot(Tob,impulse_hz,'c','linewidth',2);
% legend('hx','hy','hz');
grid on;
title([ 'ˮƽ��ż��Դ��' position '�����ų���������ӦBy']);
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
grid on;
title([ 'ˮƽ��ż��Դ��' position '�����ų���������ӦBz']);
xlabel('Time/s');
ylabel('B/(T)');
%}
%% 