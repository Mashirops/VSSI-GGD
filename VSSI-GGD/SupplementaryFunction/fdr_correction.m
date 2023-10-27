function [fdr_pvals, rejected] = fdr_correction(pvals, alpha)
    % 计算FDR校正的p值
    
    % Step 1: 对原始p值进行排序
    sorted_pvals = sort(pvals);
    
    % Step 2: 计算FDR校正的阈值
    m = length(pvals);  % 总假设数
    threshold = (1:m) / m * alpha;
    
    % Step 3: 计算FDR校正的p值
    fdr_pvals = NaN(size(pvals));
    for i = 1:m
        if sorted_pvals(i) <= threshold(i)
            fdr_pvals(pvals == sorted_pvals(i)) = sorted_pvals(i) * m / i;
        end
    end
    
    % Step 4: 根据FDR校正的p值进行假设检验
    rejected = fdr_pvals <= alpha;
end