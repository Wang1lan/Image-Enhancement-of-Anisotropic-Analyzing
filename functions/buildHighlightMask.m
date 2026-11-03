function [softMask, hardMask] = buildHighlightMask(img, maskParams)
% 构造单路亮斑保护；[softMask,hardMask] = buildHighlightMask(img,maskParams)。
% img 为 H×W×3 uint8；参数含 Te/Th/Tl、dilateSizes、numGaussian。
% 输出 H×W double softMask∈[0,1]，及由原图检测的 logical hardMask。

rgb = double(img);
Lmax = max(rgb, [], 3);
Lmin = min(rgb, [], 3);
hardMask = (Lmax >= maskParams.Te) | ...
    ((Lmax >= maskParams.Th) & ((Lmax - Lmin) <= maskParams.Tl));

core = hardMask;
[H, W] = size(core);
for dilateSize = maskParams.dilateSizes
    radius = (dilateSize - 1) / 2;
    padded = padarray(core, [radius radius], 'replicate', 'both');
    dilated = false(size(core));
    for u = 1:dilateSize
        for v = 1:dilateSize
            dilated = dilated | padded(u:u+H-1, v:v+W-1);
        end
    end
    core = dilated;
end

coreU8 = uint8(core) .* uint8(255);
blurU8 = coreU8;
g = [1 4 6 4 1];
gaussianKernel = g' * g;
for n = 1:maskParams.numGaussian
    weighted = filterShiftSum(double(blurU8), gaussianKernel, 'symmetric');
    % 每次平滑立即按非负整数“加128后右移8位”等价方式舍入。
    blurU8 = uint8(floor(weighted / 256 + 0.5));
end
softMask = double(max(coreU8, blurU8)) / 255;
end
