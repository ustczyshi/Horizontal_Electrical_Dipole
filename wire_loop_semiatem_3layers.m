%--------------------------------------------------------------------------
%��������������TEM �뺽��TEM
% ����ˮƽ��żԴ�������������й۲�,��״��ز�������Ӧ����

%--------------------------------------------------------------------------
%%
format long;
clear all;clc;close all;
%%
u0 = 4*pi*1e-7;
load parameters.txt;
sigma1 = parameters(1,2);%��һ��ĵ絼��
rou = 1./sigma1;
%% �������������
x = 100; % �շ�ˮƽƫ�ƾ࣬��y��
y = 100;
L= 1; % �������³��ȣ���x��
I = 1; % �������
%%  �뺽���շ��߶Ȳ���
z =0;% �۲������ĸ߶ȣ���������Ϊ��ֵ
h =0;% Դ�����ĸ߶�

%% �����ʺ͹۲�ʱ�������
fs = 1e5;% ������
dt = 1./fs;
t = 1/fs:1/fs:4e-2;% ʱ������


%%
% [hz_01,hz_10,hz_1_impulse,hx_01,hx_10,hx_1_impulse,hy_01,hy_10,hy_1_impulse] = Calculate_Horizontal_Electrical_Dipole(I,L,h,x,y,z,t);
[hz_01,hz_10,hz_1_impulse,hx_01,hx_10,hx_1_impulse,hy_01,hy_10,hy_1_impulse,ex_01,ex_impulse,ey_01,ey_impulse] = Calculate_Horizontal_Electrical_Dipole(I,L,h,x,y,z,t);
save('horizontal_electrical_dipole_impulse_shuzhijie','hx_1_impulse','hy_1_impulse','hz_1_impulse');
%% z���Ծ
figure;
loglog(t.*10^3,u0.*abs(hz_01),'r','Linewidth',1);
hold on
loglog(t.*10^3,u0.*abs(hz_10),'b:','Linewidth',1);
grid on;
legend('��ֵ������Ծ��Ӧ','��ֵ�⸺��Ծ��Ӧ');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz'])
xlabel('Time/(ms)')
ylabel('Bz/(A/m)');
% x���Ծ
figure;
loglog(t.*10^3,u0.*abs(hx_01),'r','Linewidth',1);
hold on
loglog(t.*10^3,u0.*abs(hx_10),'b:','Linewidth',1);
grid on;
legend('��ֵ������Ծ��Ӧ','��ֵ�⸺��Ծ��Ӧ');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx'])
xlabel('Time/(ms)')
ylabel('Bx/(A/m)');
% y ���Ծ
figure;
loglog(t.*10^3,u0.*abs(hy_01),'r','Linewidth',1);
hold on
loglog(t.*10^3,u0.*abs(hy_10),'b:','Linewidth',1);
grid on;
legend('��ֵ������Ծ��Ӧ','��ֵ�⸺��Ծ��Ӧ');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By'])
xlabel('Time/(ms)')
ylabel('By/(A/m)');
%
%% z ������
load horizontal_electrical_dipole_impulse_jiexijie.mat;
figure;
plot(t(1:end).*10^3,real(ex_01),'r','Linewidth',2);
hold on
plot(t(1:end).*10^3,(step_ex01),'k:','Linewidth',2);
grid on;
legend('��ֵ��real(ex\_01)','������step\_ex01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ex step response'])
xlabel('Time/(ms)')
ylabel('Ex/(V/m)');
figure;
plot(t(1:end).*10^3,real(ey_01),'r','Linewidth',2);
hold on
plot(t(1:end).*10^3,(step_ey01),'k:','Linewidth',2);
grid on;
legend('��ֵ��real(ey\_01)','������step\_ey01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ey step response'])
xlabel('Time/(ms)')
ylabel('Ey/(V/m)');
% ������Ӧ
% figure;
% plot(t(1:end).*10^3,ex_impulse,'r','Linewidth',2);
% hold on
% plot(t(1:end-1).*10^3,diff(step_ex01)./dt,'k:','Linewidth',2);
% grid on;
% legend('��ֵ��ex\_01','������step\_ex01');
% title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ex step response'])
% xlabel('Time/(ms)')
% ylabel('Ex/(V/m)');
% figure;
% plot(t(1:end).*10^3,abs(ey_01),'r','Linewidth',2);
% hold on
% plot(t(1:end).*10^3,abs(step_ey01),'k:','Linewidth',2);
% grid on;
% legend('��ֵ��ey\_01','������step\_ey01');
% title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ey step response'])
% xlabel('Time/(ms)')
% ylabel('Ey/(V/m)');
%
%%
%{
figure;
plot(t(1:end).*10^3,u0.*(hz_1_impulse),'r','Linewidth',2);
hold on
plot(t(1:end-1).*10^3,u0.*(diff(hz_01)./dt),'b','Linewidth',2);
grid on;
plot(t(1:end).*10^3,u0.*(impulse_hz),'k:','Linewidth',2);
grid on;
legend('��ֵ��hz\_1\_impulse','��ֵ��diff(hz\_01)./dt','������impulse\_hz');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz impulse response'])
xlabel('Time/(ms)')
ylabel('Bz/(T)');

figure;
plot(t(1:end).*10^3,u0.*(hx_1_impulse),'r','Linewidth',2);
hold on
plot(t(1:end-1).*10^3,u0.*(diff(hx_01)./dt),'b','Linewidth',2);
grid on;
plot(t(1:end).*10^3,u0.*(impulse_hx),'k:','Linewidth',2);
grid on;
legend('��ֵ��hx\_1\_impulse','��ֵ��diff(hx\_01)./dt','������impulse\_hx');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx impulse response'])
xlabel('Time/(ms)')
ylabel('Bx/(T)');
figure;
plot(t(1:end).*10^3,u0.*(hy_1_impulse),'r','Linewidth',2);
hold on
plot(t(1:end-1).*10^3,u0.*(diff(hy_01)./dt),'b','Linewidth',2);
grid on;
plot(t(1:end).*10^3,u0.*(impulse_hy),'k:','Linewidth',2);
grid on;
legend('��ֵ��hy\_1\_impulse','��ֵ��diff(hy\_01)./dt','������impulse\_hy');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By impulse response'])
xlabel('Time/(ms)')
ylabel('By/(T)');
%}
%%

%% z ������
%{
figure;
plot(t(1:end).*10^3,u0.*hz_1_impulse,'r','Linewidth',1);
hold on
plot(t(1:end).*10^3,u0.*hz_1_impulse,'b','Linewidth',1);
grid on;
legend('��ֵ����������Ӧ','��ֵ�⸺������Ӧ');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz'])
xlabel('Time/(ms)')
ylabel('Bz/(A/m)');
% y ������
figure;
plot(t(1:end).*10^3,u0.*hy_1_impulse,'r','Linewidth',1);
hold on
plot(t(1:end).*10^3,-u0.*hy_1_impulse,'b','Linewidth',1);
grid on;
legend('��ֵ����������Ӧ','��ֵ�⸺������Ӧ');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By'])
xlabel('Time/(ms)')
ylabel('By/(A/m)');
% x ������
figure;
plot(t(1:end).*10^3,u0.*hx_1_impulse,'r','Linewidth',1);
hold on
plot(t(1:end).*10^3,-u0.*hx_1_impulse,'b','Linewidth',1);
grid on;
legend('��ֵ����������Ӧ','��ֵ�⸺������Ӧ');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx'])
xlabel('Time/(ms)')
ylabel('Bx/(A/m)');
%}
%% save data

