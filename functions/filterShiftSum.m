function out = filterShiftSum(img, kernel, padMode)
% 整体错位累加实现二维卷积；out = filterShiftSum(img,kernel,padMode)。
% img 为 H×W 或 H×W×3 double，kernel 为奇数方核；输出同尺寸 double。
% padMode 为 'replicate' 或 'symmetric'，仅填充空间维。

[H, W, ~] = size(img);
radius = (size(kernel, 1) - 1) / 2;
padded = padarray(img, [radius radius], padMode, 'both');
flippedKernel = rot90(kernel, 2);
out = zeros(size(img));
for u = 1:size(flippedKernel, 1)
    for v = 1:size(flippedKernel, 2)
        out = out + flippedKernel(u, v) .* ...
            padded(u:u+H-1, v:v+W-1, :);
    end
end
end
