% function y = GS_Trans2(t,h_v)
% G_S=load ('G_S.txt')';% ��ȡGS�任ϵ��  G_SΪ������
function y = GS_Trans2(t,h_v,G_S)
y = log(2)./t .* G_S * h_v;%  ������Ӧ
end