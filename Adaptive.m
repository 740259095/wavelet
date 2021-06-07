function X = ada(x, wname, n)


% ��ʼ������ֵ
R       = 5;                                    
alpha   = 0.1;                                  
beta    = 0.3;
delta   = DELTA(x);                             
lambda2 = 4 * delta^2 * log(R);                 

[C, S] = wavedec2(x, n, wname);                 

% ��ȡÿ��ϵ�������д���
for i = n : -1 : 1
    cH = detcoef2('h', C, S, i);                
    cV = detcoef2('v', C, S, i);                
    cD = detcoef2('d', C, S, i);                
    
    dim = size(cH);
    % �ֱ������������ϵ��
    for j = 1 : dim(1)
        for k = 1 : dim(2)
            S_jk2 = energy(cH, j, k, R);
            cH(j, k) = shrink(cH(j, k), S_jk2, alpha, beta, lambda2);
            
            S_jk2 = energy(cV, j, k, R);
            cV(j, k) = shrink(cV(j, k), S_jk2, alpha, beta, lambda2);
            
            S_jk2 = energy(cD, j, k, R);
            cD(j, k) = shrink(cD(j, k), S_jk2, alpha, beta, lambda2);
        end
    end
    
  
	k     = size(S,1) - i;
	first = S(1,1)*S(1,2) + 3 * sum(S(2:k-1, 1).*S(2:k-1, 2)) + 1;  
	add   = S(k,1)*S(k,2);                                          
    
    C(first : first + add - 1) = reshape(cH, 1, add);
    C(first + add : first + 2*add - 1) = reshape(cV, 1, add);
    C(first + 2*add : first + 3*add - 1) = reshape(cD, 1, add);
end
% �ع�ͼ��
X = waverec2(C, S, wname);                      


function delta = DELTA(x)
%��������

[C, S] = wavedec2(x, 1, 'db1');                             % С���ֽ�
d = C( prod( S(1,:) ) + 2 * prod( S(2,:) ) + 1 : end);      % HH�Ӵ�ϵ��
delta = (median( abs(d) ) / 0.6745);                          % ����delta



function S_jk2 = energy(cM, j, k, R)
%   ����С��ϵ������������
%
dim = size(cM);

% �߽��ж�
row_min = (j-1 < fix(R/2)) * (1-j) + (j-1 >= fix(R/2)) * fix(-R/2);
row_max = (dim(1)-j < fix(R/2)) * (dim(1)-j) + (dim(1)-j >= fix(R/2)) * fix(R/2);
col_min = (k-1 < fix(R/2)) * (1-k) + (k-1 >= fix(R/2)) * fix(-R/2);
col_max = (dim(2)-k < fix(R/2)) * (dim(2)-k) + (dim(2)-k >= fix(R/2)) * fix(R/2);

s = 0;
for m = row_min : row_max
    for n = col_min : col_max
        s = cM(j + m, k + n)^2 + s;
    end
end
S_jk2 = s / R^2;



function d_jk = shrink(d, S_jk2, alpha, beta, lambda2)
%   ����С��ϵ��
%
if S_jk2 >= beta * lambda2
    d_jk = d * (1 - alpha * lambda2 / S_jk2);
else
    d_jk = 0;
end