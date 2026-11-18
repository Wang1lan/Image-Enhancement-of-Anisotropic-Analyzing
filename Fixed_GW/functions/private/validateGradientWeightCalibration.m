function [GLow, GHigh, Wmin] = validateGradientWeightCalibration(GLow, GHigh, Wmin, errorId)
% 两个参考系共用原始标定域，此处不做参数量化。
values = {GLow, GHigh, Wmin};
for k = 1:numel(values)
    v = values{k};
    assert(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v), ...
        errorId, 'GLow、GHigh、Wmin 必须是有限实数标量。');
end
GLow = double(GLow);
GHigh = double(GHigh);
Wmin = double(Wmin);
assert(GLow >= 0 && GLow < GHigh && GHigh <= 0.1 && Wmin >= 0.1 && Wmin <= 0.5, ...
    errorId, '要求 0 <= GLow < GHigh <= 0.1 且 0.1 <= Wmin <= 0.5。');
% 十进制门限读作 double 后，相减最多受两个端点的表示误差影响。
assert(GHigh - GLow >= 0.003 - 2*eps(max(GHigh, 0.003)), ...
    errorId, '原始高低阈值间隔必须至少为 0.003。');
end
