%% �糡����ʹų�������м���,
% �м����Ƕ�Ӧ��Ծ��Ӧ���м���
% ��Ӧ�ĵ��е��м���C_e1~e11;
% ��Ӧ�ĵ��е��м���C_e2~e12;
% ��Ӧ�ĵ��е��м���C_h1~h11;
% ��Ӧ�ĵ��е��м���C_h2~h12;
function [ e11_1_impulse,e11_01_step,e12_1_impulse,e12_01_step,h11_1_impulse,h11_01_step,h12_1_impulse,h12_01_step] = calculate_temp(r1,r2,z,t,G_S,m2,J_1,delta_1)
u0 = 4*pi*1e-7;
% �糡������
% ��Ӧ�ĵ��е��м���C_e1~e11;
e11_1_impulse = zeros(1,length(t));
e11_01_step = zeros(1,length(t));
% ��Ӧ�ĵ��е��м���C_e2~e12;
e12_1_impulse = zeros(1,length(t));
e12_01_step = zeros(1,length(t));
% �ų�������
% ��Ӧ�ĵ��е��м���C_h1~h11;
h11_1_impulse = zeros(1,length(t));
h11_01_step = zeros(1,length(t));
% ��Ӧ�ĵ��е��м���C_h2~h12;
h12_1_impulse = zeros(1,length(t));
h12_01_step = zeros(1,length(t));

%����lambda������lambda��frequency��չ�ɶ�ά����-----��Դ�ֱ���򺺿˶��任
lambda_r1 = (1./r1) .*exp(delta_1); % ��μ���lambda���ɲ�����ĺ�����ƫ����ת��Ϊ���ֱ���lambda�����R1�ӵ���
lambda_r2 = (1./r2) .*exp(delta_1); % ��μ���lambda���ɲ�����ĺ�����ƫ����ת��Ϊ���ֱ���lambda,���R2�ӵ���

for ii=1:length(t)
     freq = (log(2)*1i/(t(ii)*2*pi))*m2;
     %-------------------------------------------------------
     %���ݵ��ƹ�ʽ�������-����ķ���ϵ��
      [lambda1_Array,frequency1_Array] = Array_trans(lambda_r1,freq); % ���R1�ӵ���
      [lambda2_Array,frequency2_Array] = Array_trans(lambda_r2,freq); % ���R2�ӵ���  
      r_TE1=calculate_r_TE(lambda_r1,freq); % ���R1�ӵ���
      r_TE2=calculate_r_TE(lambda_r2,freq); % ���R2�ӵ���  
      
      z1bar_1 = calculate_r_TM_Z1bar(lambda_r1,freq);%  ���R1�ӵ���  
      z1bar_2 = calculate_r_TM_Z1bar(lambda_r2,freq);%  ���R2�ӵ���  
      %%  ��������е��迹�� zo_bar
      w1 = 2.*pi.*frequency1_Array;%����
      w2 = 2.*pi.*frequency2_Array;%����
      
      z0_bar1 = 1i.*w1.*u0;
      z0_bar2 = 1i.*w2.*u0;
      %% ---------------------------------------------------���R1�ӵ�����м�ֵ
      % �糡�ӵ�����м�ֵ
      f_e11 = exp(lambda1_Array.*(z)).*(2.*z1bar_1-(1+r_TE1).*(z0_bar1)./lambda1_Array); 
      g_e11 = Fast_Hankel(r1,f_e11,J_1);%������
      e11_1_impulse(ii) = GS_Trans2(t(ii),g_e11,G_S); %��������Ӧʱ��
      e11_01_step(ii) = GS_Trans(t(ii),g_e11,freq,G_S);%����Ծ��Ӧʱ��
      % �ų��ӵ�����м�ֵ
      f_h11 = (1+r_TE1).*exp(lambda1_Array.*(z));% r_TM1 = 1;
      g_h11 = Fast_Hankel(r1,f_h11,J_1);%����Ծ��ӦƵ��f_h11��ͬһ�еĶ�Ӧ��ͬlambda�����˶��任�˲�ϵ����ͬƫ��������ֵ
      h11_1_impulse(ii) = GS_Trans2(t(ii),g_h11,G_S); %��������Ӧʱ��
      h11_01_step(ii) = GS_Trans(t(ii),g_h11,freq,G_S);%����Ծ��Ӧʱ��
      %% ---------------------------------------------------���R2�ӵ�����м�ֵ
      % �糡�ӵ�����м�ֵ
      f_e12 = exp(lambda2_Array.*(z)).*(2.*z1bar_2-(1+r_TE2).*(z0_bar2)./lambda2_Array); 
      g_e12 = Fast_Hankel(r2,f_e12,J_1);%������
      e12_1_impulse(ii) = GS_Trans2(t(ii),g_e12,G_S); %��������Ӧʱ��
      e12_01_step(ii) = GS_Trans(t(ii),g_e12,freq,G_S);%����Ծ��Ӧʱ��
      % �ų��ӵ�����м�ֵ
      f_h12 = (1+r_TE2).*exp(lambda2_Array.*(z));% r_TM2 = 1;
      g_h12 = Fast_Hankel(r2,f_h12,J_1);%����Ծ��ӦƵ��f_h11��ͬһ�еĶ�Ӧ��ͬlambda�����˶��任�˲�ϵ����ͬƫ��������ֵ
      h12_1_impulse(ii) = GS_Trans2(t(ii),g_h12,G_S); %��������Ӧʱ��
      h12_01_step(ii) = GS_Trans(t(ii),g_h12,freq,G_S);%����Ծ��Ӧʱ��
end
% temp_e1  = e11_01_step;
% temp_e2  = e12_01_step;
% temp_h1  = h11_01_step;
% temp_h2  = h12_01_step;


end
