function nVOL = calculateHypervolumeDiversity1(pop, lb, ub)
% 计算超体积多样性指标 (公式30-31)
% 输入：当前种群、搜索空间边界、维度
% 输出：种群超体积和标准化指标

% 计算种群边界
pop_lb = min(pop,[],1);
pop_ub = max(pop,[],1);

% 计算极限超体积 (公式30)
V_lim = sqrt(prod(ub - lb));

% 计算种群超体积 (公式31)
V_pop = sqrt(prod((pop_ub - pop_lb)/2));

% 标准化指标
nVOL = sqrt(V_pop / V_lim);
end
