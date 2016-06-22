function rte=calculate_r_TE(lambda,freq)
%% ����TE���ڷֲ����еķ���ϵ��
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
kk = -1i * 2 * pi * frequency_Array * mu_0;
u1_star =  (lambda_2 - kk .* sigma(N)).^0.5;% 
for k= N-1:-1:1 %% n��ֻ��Ҫ����n-1��
    u_n = (lambda_2 - kk.*sigma(k)).^0.5;% ��n���u_n��
%     u1_star = u_n .*( u1_star + u_n .* exp(-2*u_n * d(k)))./(u_n + u1_star .* exp(-2*u_n * d(k)));
    temp = tanh(u_n * d(k));
    u1_star = u_n .*( u1_star + u_n .* temp)./(u_n + u1_star .* temp);
end
rte = (lambda_Array-u1_star)./(lambda_Array+u1_star);

end
