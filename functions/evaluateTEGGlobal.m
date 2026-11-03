function teg = evaluateTEGGlobal(img, imgEnh, hardMask, evalParams)
% 计算非高光区域的梯度均值比；输入为同尺寸 RGB uint8 图像、logical 掩模。
% evalParams 指定高光距离和剔除边宽；输出为标量，无有效基准时为 NaN。

    textureROI = bwdist(hardMask) > evalParams.highlightDistance;
    b = evalParams.borderWidth;
    textureROI([1:b, end-b+1:end], :) = false;
    textureROI(:, [1:b, end-b+1:end]) = false;
    teg = NaN;
    if ~any(textureROI(:))
        return;
    end

    gradOrig = gradientMagnitude(img);
    gradEnh = gradientMagnitude(imgEnh);
    meanOrig = mean(gradOrig(textureROI));
    if meanOrig > eps('double')
        teg = mean(gradEnh(textureROI)) / meanOrig;
    end
end

function gradient = gradientMagnitude(img)
    rgb = double(img);
    y = (77 * rgb(:, :, 1) + 150 * rgb(:, :, 2) + 29 * rgb(:, :, 3)) / 256;
    sx = [-1 0 1; -2 0 2; -1 0 1];
    gx = filterShiftSum(y, sx, 'symmetric');
    gy = filterShiftSum(y, sx', 'symmetric');
    gradient = hypot(gx, gy);
end
