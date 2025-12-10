function [X_selected, fit_selected,sim_selected_idx] = selectPopulation(X, fit, t, Max_iter, N1,dim)
    % 输入：
    % X - 种群矩阵 (N x dim)
    % fit - 适应度值 (N x 1)
    % t - 当前迭代次数
    % Max_iter - 最大迭代次数
    % N1 - 需要选择的个体数
    % 输出：
    % X_selected - 选出的个体矩阵
    % fit_selected - 对应的适应度值
    
    % 1. 按适应度排序选择前15个个体
    [fit_sorted, sort_idx] = sort(fit, 'ascend');
    X_top15 = X(sort_idx(1:15), :);
    fit_top15 = fit_sorted(1:15);
    
    % 2. 计算余弦相似度矩阵（相对于最优个体）
    x_best = X_top15(1,:); % 适应度最好的个体
    cosine_sim = zeros(15, 1);
    for i = 1:15
       sumx=0;sum1=0;sum2=0;
        for j=1:dim
            sumx=sumx+(x_best(1,j)* X_top15(i,j));
            sum1=sum1+(x_best(1,j)^2);
            sum2=sum2+( X_top15(i,j)^2);
        end
        cosine_sim(i,1)=sumx/(sqrt(sum1)*sqrt(sum2));
    end
    
    % 3. 按余弦相似度排序
    [~, sim_sort_idx] = sort(cosine_sim);
    
    % 4. 动态计算划分点
    num_from_sim =   N1-round(0.6*N1 - 0.2*N1*t/Max_iter);
    num_from_fit = round(0.6*N1 - 0.2*N1*t/Max_iter);
    
    % 边界检查
    num_from_sim = max(1, min(num_from_sim, 14)); % 至少选1个，最多选14个
    num_from_fit = max(1, min(num_from_fit, 14)); % 至少选1个，最多选14个
    
    % 5. 选择个体
    % 5.1 从相似度高的前num_from_sim个
    sim_selected_idx = sim_sort_idx(1:num_from_sim);
    
    % 5.2 从剩余的个体中选择适应度好的前num_from_fit个
    remaining_idx = setdiff(1:15, sim_selected_idx);
    [~, remain_fit_sort] = sort(fit_top15(remaining_idx), 'ascend');
    fit_selected_idx = remaining_idx(remain_fit_sort(1:num_from_fit));
    
    % 6. 合并选中的个体
    final_selected_idx = [sim_selected_idx; fit_selected_idx'];
    X_selected = X_top15(final_selected_idx, :);
    fit_selected = fit_top15(final_selected_idx);
    
    % 验证总数是否正确
    assert(length(final_selected_idx) == N1, '选择个体数量错误');
end