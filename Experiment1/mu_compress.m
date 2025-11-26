function y = mu_compress(x, x_max, mu)
    % MU_COMPRESS mu律压缩函数
    %   y = MU_COMPRESS(x, x_max, mu) 对输入信号x进行mu律压缩
    %
    % 输入：
    %   x     - 输入信号（向量）
    %   x_max - 输入信号的最大绝对值
    %   mu    - mu律压缩参数（通常为255）
    % 输出：
    %   y     - 压缩后的信号（列向量）

    if nargin < 3
        mu = 255; % Default mu value
    end
    x = x(:);
    y = x_max .* log(1 + mu .* abs(x) ./ x_max) ./ log(1 + mu) .* sign(x);
end
