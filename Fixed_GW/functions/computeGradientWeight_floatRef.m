function [Wfloat, debugFloat] = computeGradientWeight_floatRef(Y12, GLow, GHigh, Wmin)
%COMPUTEGRADIENTWEIGHT_FLOATREF 同一 Y12 输入的原图梯度线性浮点参考。
%   [Wfloat,debugFloat] = computeGradientWeight_floatRef(Y12,GLow,GHigh,Wmin)
%   Y12 为至少 3x3 的二维 uint16，[0,4095]；参数与定点版本相同。
%   需要将项目原 functions 目录加入路径。输出 double，范围 [Wmin,1]。
%   浮点参考采用仓库 Sobel 卷积及 9x9 replicate 均值。

assert(isa(Y12, 'uint16') && isreal(Y12) && ismatrix(Y12) && ...
    all(size(Y12) >= 3) && all(Y12(:) <= uint16(4095)), ...
    'computeGradientWeight_floatRef:Input', ...
    'Y12 必须是至少 3x3、数值不超过4095的二维 uint16 图像。');
prepareGradientWeightFxParams(GLow, GHigh, Wmin);
opts = struct('inputMax', 4095, 'GLow', GLow, 'GHigh', GHigh, ...
    'Wmin', Wmin, 'denoise', struct('enable', false), ...
    'mapping', struct('type', 'linear'), 'debugLevel', 'off');
if nargout > 1
    opts.debugLevel = 'maps';
    [W0, reference] = computeGradientWeight(Y12, opts);
else
    W0 = computeGradientWeight(Y12, opts);
end
Wfloat = imboxfilt(W0, 9, 'Padding', 'replicate');
if nargout > 1
    debugFloat = struct('Gx', reference.Gx, 'Gy', reference.Gy, ...
        'G', reference.G, 'W0', W0, 'W_smooth', Wfloat);
end
end
