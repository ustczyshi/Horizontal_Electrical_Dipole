%--------------------------------------------------------------------------
%��������������TEM �뺽��TEM
% ����ˮƽ��żԴ�������������й۲�

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
%% �糡����
%  x�����
% E_x_impulse=zeros(1,length(t));% ��������Ӧʱ��
% E_x_step = E_x_impulse;%����Ծ
% E_x_impulse1=zeros(1,length(t));% ��������Ӧʱ��
% E_x_step1 = E_x_impulse1;%����Ծ
% %  y����
% E_y_impulse=zeros(1,length(t));% ��������Ӧʱ��
% E_y_step = E_y_impulse;%����Ծ
% E_y_impulse1=zeros(1,length(t));% ��������Ӧʱ��
% E_y_step1 =E_y_impulse1;%����Ծ
% %% �ų�����
% %  z�����
% h_z_impulse=zeros(1,length(t));% ��������Ӧʱ��
% h_z_step = h_z_impulse;%����Ծ
% h_z_impulse1=zeros(1,length(t));% ��������Ӧʱ��
% h_z_step1 = h_z_impulse1;%����Ծ
% %  x�����
% h_x_impulse=zeros(1,length(t));% ��������Ӧʱ��
% h_x_step = h_x_impulse;%����Ծ
% h_x_impulse1=zeros(1,length(t));% ��������Ӧʱ��
% h_x_step1 = h_x_impulse1;%����Ծ
% %  y����
% h_y_impulse=zeros(1,length(t));% ��������Ӧʱ��
% h_y_step = h_y_impulse;%����Ծ
% h_y_impulse1=zeros(1,length(t));% ��������Ӧʱ��
% h_y_step1 = h_y_impulse1;%����Ծ
%%
% [hz_01,hz_10,hz_1_impulse,hx_01,hx_10,hx_1_impulse,hy_01,hy_10,hy_1_impulse] = Calculate_Horizontal_Electrical_Dipole(I,L,h,x,y,z,t);
% [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_impulse,ey_01,ey_impulse] = Calculate_Horizontal_Electrical_Dipole(I,L,h,x,y,z,t)
[hz_01,hz_10,hz_1_impulse,hx_01,hx_10,hx_1_impulse,hy_01,hy_10,hy_1_impulse,ex_01,ex_impulse,ey_01,ey_impulse] = Calculate_Horizontal_Electrical_Dipole(I,L,h,x,y,z,t);
% save('horizontal_electrical_dipole_impulse_shuzhijie','hx_1_impulse','hy_1_impulse','hz_1_impulse');
%% z���Ծ
%{
%    ��ͼ��
hz_10 = hz_01(end) - hz_01;
figure;
plot(t.*10^3,u0.*hz_01,'r','Linewidth',1);
hold on
plot(t.*10^3,u0.*hz_10,'b','Linewidth',1);
grid on;
legend('��ֵ������Ծ��Ӧ','��ֵ�⸺��Ծ��Ӧ');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz'])
xlabel('Time/(ms)')
ylabel('Bz/(A/m)');
% x���Ծ
%��ͼ��
% hx_10 = hx_01(end) - hx_01;
figure;
plot(t.*10^3,u0.*(hx_01),'r','Linewidth',1);
hold on
plot(t.*10^3,-u0.*(hx_01),'b','Linewidth',1);
grid on;
legend('��ֵ������Ծ��Ӧ','��ֵ�⸺��Ծ��Ӧ');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx'])
xlabel('Time/(ms)')
ylabel('Bx/(A/m)');
% y ���Ծ
hy_10 = hy_01(end) - hy_01;
figure;
plot(t.*10^3,u0.*(hy_01),'r','Linewidth',1);
hold on
plot(t.*10^3,u0.*(hy_10),'b','Linewidth',1);
grid on;
legend('��ֵ������Ծ��Ӧ','��ֵ�⸺��Ծ��Ӧ');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By'])
xlabel('Time/(ms)')
ylabel('By/(A/m)');
%}
%% z ������
load horizontal_electrical_dipole_impulse_jiexijie.mat;
figure;
loglog(t(1:end).*10^3,abs(real(ex_01)),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,abs(step_ex01),'k:','Linewidth',2);
grid on;
legend('��ֵ��ex\_01','������step\_ex01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ex step response'])
xlabel('Time/(ms)')
ylabel('Ex/(V/m)');
figure;
loglog(t(1:end).*10^3,abs(real(ey_01)),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,abs(step_ey01),'k:','Linewidth',2);
grid on;
legend('��ֵ��ey\_01','������step\_ey01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ey step response'])
xlabel('Time/(ms)')
ylabel('Ey/(V/m)');
% ������
% ex_error = abs(abs(ex_01)-abs(step_ex01))./abs(step_ex01);
% ey_error = abs(abs(ey_01)-abs(step_ey01))./abs(step_ey01);
% figure;
% loglog(t.*1e3,ex_error.*100,'r','linewidth',2);
% hold on;
% loglog(t.*1e3,ey_error.*100,'k:','linewidth',2);
% grid on;
% legend('ex','ey');
% title('��ֵ��ͽ������������');
% xlabel('Time/(ms)')
% ylabel('error/(%)');
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
figure;
loglog(t(1:end).*10^3,u0.*abs(hz_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end-1).*10^3,u0.*abs(diff(hz_01)./dt),'b','Linewidth',2);
grid on;
loglog(t(1:end).*10^3,u0.*abs(impulse_hz),'k:','Linewidth',2);
grid on;
legend('��ֵ��hz\_1\_impulse','��ֵ��diff(hz\_01)./dt','������impulse\_hz');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz impulse response'])
xlabel('Time/(ms)')
ylabel('Bz/(T)');

figure;
loglog(t(1:end).*10^3,u0.*abs(hx_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end-1).*10^3,u0.*abs(diff(hx_01)./dt),'b','Linewidth',2);
grid on;
loglog(t(1:end).*10^3,u0.*abs(impulse_hx),'k:','Linewidth',2);
grid on;
legend('��ֵ��hx\_1\_impulse','��ֵ��diff(hx\_01)./dt','������impulse\_hx');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx impulse response'])
xlabel('Time/(ms)')
ylabel('Bx/(T)');
figure;
loglog(t(1:end).*10^3,u0.*abs(hy_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end-1).*10^3,u0.*abs(diff(hy_01)./dt),'b','Linewidth',2);
grid on;
loglog(t(1:end).*10^3,u0.*abs(impulse_hy),'k:','Linewidth',2);
grid on;
legend('��ֵ��hy\_1\_impulse','��ֵ��diff(hy\_01)./dt','������impulse\_hy');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By impulse response'])
xlabel('Time/(ms)')
ylabel('By/(T)');
% ������
hz_error = abs(abs(hz_1_impulse)-abs(impulse_hz))./abs(impulse_hz);
hx_error = abs(abs(hx_1_impulse)-abs(impulse_hx))./abs(impulse_hx);
hy_error = abs(abs(hy_1_impulse)-abs(impulse_hy))./abs(impulse_hy);
figure;
loglog(t.*1e3,hz_error.*100,'r','linewidth',2);
hold on;
loglog(t.*1e3,hx_error.*100,'b','linewidth',2);
hold on;
loglog(t.*1e3,hy_error.*100,'k:','linewidth',2);
grid on;
legend('hz','hx','hy');
title('��ֵ��ͽ������������');
xlabel('Time/(ms)')
ylabel('error/(%)');

%% save data

