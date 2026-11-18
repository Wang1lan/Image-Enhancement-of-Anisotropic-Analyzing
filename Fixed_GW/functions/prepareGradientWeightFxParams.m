function params = prepareGradientWeightFxParams(GLow, GHigh, Wmin)
%PREPAREGRADIENTWEIGHTFXPARAMS 原始浮点标定离线转换为四个整数寄存器。
%   params = prepareGradientWeightFxParams(GLow,GHigh,Wmin) 使用 Q17 斜率。
%   原始标定随参数保留供浮点验收。
[GLow,GHigh,Wmin] = validateGradientWeightCalibration(GLow,GHigh,Wmin, ...
    'prepareGradientWeightFxParams:Parameters');
w = uint16(round(4095*Wmin));                         % u12 [410,2048]
low = uint32(round(2096640*GLow));                    % u18, 32760*64
high = uint32(round(2096640*GHigh));                  % u18 <=209664
K_Q17 = gradientWeightSlopeWord(w, high-low);
params = loadGradientWeightFxRegs(w,low,high,K_Q17);
params.GLow = GLow;
params.GHigh = GHigh;
params.Wmin = Wmin;
end
