% function y = GS_Trans(t,h_v,freq)
% G_S=load ('G_S.txt')';% G_SΪ������
function y = GS_Trans(t,h_v,freq,G_S)
y = log(2)./(-1i*2*pi*freq)./t .* G_S * h_v; % ��Ծ��Ӧ
end