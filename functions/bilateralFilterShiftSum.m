function out = bilateralFilterShiftSum(src, windowSize, sigmaSpatial, sigmaRange)
% 整体错位累加双边滤波；src/out 为同尺寸 H×W 或 H×W×3 double。
% 幅值单位为 8 位灰度级；RGB 使用联合均方距离，边界为 replicate。

[H, W, ~] = size(src);
radius = (windowSize - 1) / 2;
padded = padarray(src, [radius radius], 'replicate', 'both');
[dx, dy] = meshgrid(-radius:radius);
spatialWeight = exp(-0.5 * ((dx ./ sigmaSpatial).^2 + (dy ./ sigmaSpatial).^2));
numerator = zeros(size(src));
denominator = zeros(H, W);
for u = 1:windowSize
    for v = 1:windowSize
        shifted = padded(u:u+H-1, v:v+W-1, :);
        rangeDistance = mean(((shifted - src) ./ sigmaRange).^2, 3);
        weight = spatialWeight(u, v) .* exp(-0.5 * rangeDistance);
        numerator = numerator + weight .* shifted;
        denominator = denominator + weight;
    end
end
out = numerator ./ denominator;
end
