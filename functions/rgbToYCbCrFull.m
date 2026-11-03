function ycc = rgbToYCbCrFull(rgb)
% 全范围 RGB→YCbCr；输入输出为 H×W×3 double，不裁剪或舍入。
% RGB/Y 单位为 8 位灰度级，Cb/Cr 中性值为 128。

Kr = 0.299; Kb = 0.114; Kg = 1 - Kr - Kb;
R = rgb(:,:,1); G = rgb(:,:,2); B = rgb(:,:,3);
Y = Kr * R + Kg * G + Kb * B;
Cb = 128 + (B - Y) / (2 * (1 - Kb));
Cr = 128 + (R - Y) / (2 * (1 - Kr));
ycc = cat(3, Y, Cb, Cr);
end
