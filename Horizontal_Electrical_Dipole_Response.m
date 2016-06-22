%% ˮƽ��ż��Դ�ڵ��棨h=z=0���������ĳ�
% step_ex10,step_ey10 ˮƽ�糡������Ծ��Ӧ
% impluse_hx,impulse_hy ˮƽ�ų���������Ӧ
% impulse_hz ��ֱ�ų���������Ӧ
function [step_ex01,step_ey01,impluse_hx,impulse_hy,impulse_hz] = Horizontal_Electrical_Dipole_Response(u0,rou,m,xr,yr,zr,Tob)
if zr ~= 0
    zr = 0;
end
r = (xr.^2+yr.^2).^0.5;
% alph = (u0./rou).^0.5.*r;
sita = (u0./rou./Tob).^(0.5)./2;
temp = sita.*r;
%% �糡�Ľ�Ծ��Ӧ��ϵ��
ex_coff = m.*rou./2./pi./r.^3;
ey_coff = m.*rou./2./pi./r.^5;
%% �ų���������Ӧ��ϵ��
hx_coff = -m.*xr.*yr./4./pi./r.^4./Tob.*exp(-temp.^2/2);
hy_coff = m.*sita.^2.*exp(-temp.^2/2).*rou./u0./pi./r^2;
hz_coff = m.*yr.*rou./2./pi./u0./r.^5;
% ��һ������Bessel����
I0 = besseli(0,(temp.^2)/2);
I1 = besseli(1,(temp.^2)/2);
% ex ey ��Ծ��Ӧ
step_ex01 =  ex_coff.*((3.*(xr./r).^2-1)-erf(temp)+2.*temp.*exp(-temp.^2)./(pi).^0.5);
step_ey01 = ey_coff.*3.*xr.*yr.*ones(1,length(Tob));
% ������Ӧ
impluse_hx = hx_coff.*((4+temp.^2).*I1-temp.^2.*I0);
impulse_hy = hy_coff.*(3.*I1+temp.^2.*(I1-I0)-xr.^2./r.^2.*((temp.^2+4).*I1-temp.^2.*I0));
impulse_hz = hz_coff .*(3.*erf(temp)-2./pi.^0.5.*temp.*(3+2.*temp.^2).*exp(-temp.^2));
% save data
save(['Horizontal_Electrical_Dipole_Response_rou' num2str(rou) '_xr' num2str(xr) '_yr' num2str(yr) '_zr' num2str(zr) '.mat'],...
    'step_ex01','step_ey01','impluse_hx','impulse_hy','impulse_hz');

end