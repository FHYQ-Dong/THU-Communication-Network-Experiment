function x = mu_expand(y, x_max, mu)
    % MU_EXPAND mu律扩展函数
    %   x = MU_EXPAND(y, y_max, mu) 对输入信号y进行mu律扩展
    %
    % 输入：
    %   y     - 输入信号（向量）
    %   y_max - 输入信号的最大绝对值
    %   mu    - mu律扩展参数（通常为255）
    % 输出：
    %   x     - 扩展后的信号（列向量）

    if nargin < 3
        mu = 255; % Default mu value
    end
    y = y(:);
    x = sign(y) .* x_max ./ mu .* ((1+mu) .^ (abs(y) ./ x_max) - 1);
end
