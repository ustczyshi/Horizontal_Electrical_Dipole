%����ȡ����lambda��i��������lambda��Ƶ����չΪ��ά����
%ʹ��ͬ��Ƶ�ʲ�ͬ��ͬһ��Ƶ����ͬ����ͬ��lambda��ͬ��ͬһ��lambda��ͬ��
% lambda: ��������
% freq����������
function [lmd_array,fre_array] = Array_trans(lambda, freq)
if (size(lambda,1)>1) &&(size(lambda,2)>1)
    error('input lambda is not a vector');  % ����������򱨴�
else
    if (size(lambda,2)==1) %% �������Ϊ��������ת��Ϊ������
        lambda = lambda';
    end
end
    
if (size(freq,1)>1) &&(size(freq,2)>1)
    error('input freq is not a vector'); 
else
    if (size(freq,2)==1)
        freq = freq';
    end
end
col= length(lambda);
row = length(freq);
lmd_array = repmat(lambda,row,1);%��ͬ��lambda��ͬ��ͬһ��lambda��ͬ��
fre_array = repmat(freq',1,col);%ʹ��ͬ��Ƶ�ʲ�ͬ��ͬһ��Ƶ����ͬ��
% lmd_array = ones(row,1) * lambda;
% fre_array = freq' * ones(1,col);
end