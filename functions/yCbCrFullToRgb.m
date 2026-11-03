function rgb = yCbCrFullToRgb(ycc)
% 全范围 YCbCr→RGB；输入输出为 H×W×3 double，不裁剪或舍入。
% 与 rgbToYCbCrFull 配套，Cb/Cr 中性值为 128。

Kr = 0.299; Kb = 0.114; Kg = 1 - Kr - Kb;
Y = ycc(:,:,1); Cb = ycc(:,:,2) - 128; Cr = ycc(:,:,3) - 128;
R = Y + 2 * (1 - Kr) * Cr;
B = Y + 2 * (1 - Kb) * Cb;
G = (Y - Kr * R - Kb * B) / Kg;
rgb = cat(3, R, G, B);
end
