function [fdr_pvals, rejected] = fdr_correction(pvals, alpha)
    % ����FDRУ����pֵ
    
    % Step 1: ��ԭʼpֵ��������
    sorted_pvals = sort(pvals);
    
    % Step 2: ����FDRУ������ֵ
    m = length(pvals);  % �ܼ�����
    threshold = (1:m) / m * alpha;
    
    % Step 3: ����FDRУ����pֵ
    fdr_pvals = NaN(size(pvals));
    for i = 1:m
        if sorted_pvals(i) <= threshold(i)
            fdr_pvals(pvals == sorted_pvals(i)) = sorted_pvals(i) * m / i;
        end
    end
    
    % Step 4: ����FDRУ����pֵ���м������
    rejected = fdr_pvals <= alpha;
end