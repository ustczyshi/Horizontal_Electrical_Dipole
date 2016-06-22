%%    ������ɢ��Ϊ���
% ������Ҫ�Ǹ����շ�ƫ�ƾ�ͷ���Դ�ߴ磬���������ɢ������С�߶�
% Դ��λ���ڣ�-L/2,L/2��֮�䣻
% �۲���ڣ�x,y,z��;
% �۲���Դ��ˮƽ��ֱ����Ϊ|y|
% �۲���Դ��ɢ�������ż���ĵ�ˮƽ����Ϊ
% min(|y|,((x-L/2)^2+y^2)^(0.5),((x+L/2)^2+y^2)^(0.5))
% |y| ��Ӧ�۲��x��(-L/2,L/2)֮��
% ((x-L/2)^2+y^2)^(0.5),((x+L/2)^2+y^2)^(0.5) ��Ӧ�۲���λ����|x|>L/2

%% ��ɢ���ı�׼
% ���շ����볬��Դ���ȵ�10���ǿ�����ΪԴΪ��ż��Դ
% ����ȡ΢Ԫdx�߶�Ϊ�۲��෢�䳤������̾����1/20Ϊ��׼
% dx = dmin /20;
% ��ɢ���γɵĵ�żԴ�ĸ���N = L/dx;
% r ���ظ���ż��Դ���ľ�۲��ľ���
function [N, dx,x_center,dmin,r] = Discrete_Source(x,y,z,L)
    r1 = ((x+L/2)^2+y^2)^(0.5);
    r2 = ((x-L/2)^2+y^2)^(0.5);
    %% ����۲�۲��೤����Դ����̾���
    if x >L/2
        dmin = r2;
    elseif x < -L/2
        dmin = r1;
    else
        dmin = abs(y);
    end
    zoomout = 20; 
if dmin >= zoomout*L
         r =  (x^2+y^2).^(0.5);
         dx = L;
         N = 1;
         x_center = 0;
         display('������Ϊ����ԴΪˮƽ��ż��Դ������ָ��С΢Ԫ');
         return;
     else
        % if the distance from the observation point to the source < 20 times the long of the source
    dx = dmin ./ zoomout;
end
%     if dx > 1
%         dx = 1;
    N =floor(ceil(L./dx)./2).*2+1 ;% the moment cell numbers
    dx = L/N; % real scale of the moment cell
    % ��Ϊ������Ϊʹ�м�ĵ�żԴ�е�����Դ�㴦��Ҳ��Ϊ�˶Գƿ���
%     x_p = ones(1,(N-1)./2);% locate + x axis
    x_p = dx.*(1:(N-1)/2); % the center position of the moment cell located + x axis direction   
    x_n = -x_p;% the center position of the moment cell located - x axis direction
    x_center = [x_n, 0 , x_p];% the center position of the moment cell
    r = ((x-x_center).^2+y.^2).^0.5;
    %% ���ƽ��ͼ���ֽ�Ч��
    figure;
    x_axis = (-L/2:dx:L/2);
    y_axis = zeros(1,length(x_axis));
    z_axis = zeros(1,length(x_axis));
    y_center = zeros(1,length(x_center));
    z_center = zeros(1,length(x_center));
    plot3(x_axis,y_axis,z_axis,'r','Linewidth',1);% ������Դ�ķ�Χ
    hold on;
    plot3(x_center,y_center,z_center,'bx','linewidth',1);% ������Դ�ֳɵĵ�ż��Դ������λ��
    hold on;
    plot3(x,y,z,'bo','linewidth',1);
    axis([-L L -2*L 2*L -10 10 ]);
    grid on;
    title('simulate model');
    xlabel('x-source scale/m');
    ylabel('y-offset/m');
    zlabel(' z-hight of the receiver/m');
%     end
    
end