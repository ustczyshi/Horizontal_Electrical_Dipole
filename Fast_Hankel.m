function y=Fast_Hankel(r,sum,J_zero)
%rΪ�շ��࣬
load parameters.txt;% ��ȡ����
m = parameters(1,1); % ��ȡ�ž�
y = m/r  *  sum * J_zero';% �þ���˷�������ֵ���֣�y�Ǹ���������ÿһ�ж�Ӧ��ͬ��Ƶ��
%./r���ɺ��˶��任�����е�ת������������
end