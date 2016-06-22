function rtm_z1bar=calculate_r_TM_Z1bar(lambda,freq)
%% ����TM���ڷֲ����еķ���ϵ��,���ڵ���Դ�ַ���,���ڴ�����������������rtm=1
%�������
mu_0 = 4*pi*1e-7;
%��ȡ����
load parameters.txt;
sigma = parameters(:,2)';% �絼�� 1~n��Ԫ�ش���1~n��ĵ絼�� ������ת��Ϊ������
d =  parameters(:,3)'; % ���
N=length(sigma);
[lambda_Array,frequency_Array] = Array_trans(lambda,freq); % ������ת��Ϊ����
% forѭ���ж�μ��㲿�ֵ�������������м����ݣ���С������
lambda_2 = lambda_Array.^2; 
% kk = -1i * 2 * pi * frequency_Array * mu_0;
% u1_star =  (lambda_2 - kk * 1./rho(N)).^0.5;% 
kk = -1i * 2 * pi * frequency_Array * mu_0;%ʹ��ͬ��Ƶ�ʲ�ͬ��ͬһ��Ƶ����ͬ
% u1_star =  (lambda_2 - kk .* sigma(N)).^0.5;% �����ײ��u_n
%%   �������-------------------------------------20160511
z1_star =  (lambda_2 - kk .* sigma(N)).^0.5./sigma(N);
for k= N-1:-1:1 %% n��ֻ��Ҫ����n-1��
     u_n = (lambda_2 - kk.*sigma(k)).^0.5;% ��n���u_n��
    z_n = u_n./sigma(k);% ��n���z_n,�������------------20160511
%     z1_star = z_n .*( z1_star + z_n .* exp(-2*u_n * d(k)))./(z_n + z1_star .* exp(-2*u_n * d(k)));
    temp = tanh(u_n * d(k));
    z1_star = z_n .*( z1_star + z_n .* temp)./(z_n + z1_star .* temp);
end
rtm_z1bar = z1_star;

end
