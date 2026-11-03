function noise = evaluateNoisePCA(img, imgEnh, roi, noiseParams)
% 在共享 ROI 和 patch 中心上估计 PCA 噪声；图像为同尺寸 RGB uint8。
% roi 为一基整数 [x y width height]；sigma 单位为 8 位灰度级。
% 返回 sigmaOrig、sigmaEnh、noiseGain、noiseGainDB、roiPos、validPatchNum。

    p = noiseParams.patchSize;
    stride = noiseParams.stride;
    tailNum = noiseParams.tailNum;
    assert(isscalar(p) && isfinite(p) && p >= 3 && p == round(p) && mod(p, 2) == 1, ...
        'patchSize 必须是不小于 3 的奇数。');
    assert(isscalar(stride) && isfinite(stride) && stride == round(stride) && ...
        stride >= 1 && stride < p, 'stride 必须为 1 到 patchSize-1 的整数。');
    assert(isscalar(tailNum) && isfinite(tailNum) && tailNum == round(tailNum) && ...
        tailNum >= 1 && tailNum <= p^2, 'tailNum 必须为 1 到 patchSize^2 的整数。');
    assert(isequal(size(img), size(imgEnh)), '原图与增强图尺寸必须一致。');
    assert(numel(roi) == 4 && all(isfinite(roi(:))) && ...
        all(roi(:) == round(roi(:))) && all(roi(:) > 0), ...
        'ROI 必须为 [x y width height] 正整数坐标。');
    roi = double(reshape(roi, 1, []));
    assert(roi(3) >= p && roi(4) >= p && ...
        roi(1) + roi(3) - 1 <= size(img, 2) && ...
        roi(2) + roi(4) - 1 <= size(img, 1), ...
        'ROI 必须位于图像内并能容纳完整 patch。');

    halfSize = floor(p / 2);
    rows = roi(2) + halfSize:stride:roi(2) + roi(4) - 1 - halfSize;
    columns = roi(1) + halfSize:stride:roi(1) + roi(3) - 1 - halfSize;
    [columnGrid, rowGrid] = meshgrid(columns, rows);
    patchPos = [rowGrid(:), columnGrid(:)];
    patchNum = size(patchPos, 1);
    assert(patchNum >= p^2 + 1, 'ROI 的 patch 数量不足，请选择更大的 ROI。');

    noise.sigmaOrig = estimateSigma(img, patchPos, p, tailNum);
    noise.sigmaEnh = estimateSigma(imgEnh, patchPos, p, tailNum);
    noise.noiseGain = NaN;
    noise.noiseGainDB = NaN;
    if noise.sigmaOrig > eps('double')
        noise.noiseGain = noise.sigmaEnh / noise.sigmaOrig;
        noise.noiseGainDB = 20 * log10(noise.noiseGain);
    end
    noise.roiPos = roi;
    noise.validPatchNum = patchNum;
end

function sigma = estimateSigma(img, patchPos, patchSize, tailNum)
    rgb = double(img);
    y = (77 * rgb(:, :, 1) + 150 * rgb(:, :, 2) + 29 * rgb(:, :, 3)) / 256;
    patchNum = size(patchPos, 1);
    halfSize = floor(patchSize / 2);
    x = zeros(patchSize^2, patchNum);
    for k = 1:patchNum
        row = patchPos(k, 1);
        column = patchPos(k, 2);
        patch = y(row-halfSize:row+halfSize, column-halfSize:column+halfSize);
        x(:, k) = patch(:);
    end

    % 按行跨 patch 去均值，保留仓库的协方差定义。
    x = x - mean(x, 2);
    covariance = x * x' / (patchNum - 1);
    covariance = (covariance + covariance') / 2;
    eigenvalues = sort(real(eig(covariance)), 'ascend');
    eigenvalues = max(eigenvalues, 0);
    sigma = sqrt(median(eigenvalues(1:tailNum)));
end
