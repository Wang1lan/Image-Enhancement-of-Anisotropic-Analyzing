function filter = bandpassFilter_B(fRadius, sigmaHighFreq, sigmaLowFreq)
% bandpassFilter_B  方案B：零直流目标频响的双高斯差带通核。
% 用法：h = bandpassFilter_B(10, 0.5, 0.1);
% 输入：
%   fRadius       - 核半径，输出尺寸为 (2*fRadius+1)^2
%   sigmaHighFreq - 较宽频域高斯的尺度
%   sigmaLowFreq  - 较窄频域高斯的尺度（必须小于 sigmaHighFreq）
% 输出：
%   filter        - 实数二维带通核，sum(filter(:))≈0，正系数之和≈1
% 依赖：fsamp2（与原仓库 bandpassFilter.m 相同）
%
% 与方案A的区别：
%   直接设 H(r)=exp(-r^2/sigmaHighFreq^2)-exp(-r^2/sigmaLowFreq^2)，
%   理论直流响应 H(0)=0；核生成后仅均匀修正有限尺寸造成的残余DC，
%   不再将中心系数减1，因此不人为引入全频段的常量频响项。
%   原方案A中的 scaleLowFreq 在此固定为1，不再需要作为输入。

if nargin < 1 || isempty(fRadius)
    fRadius = 10;
end
if nargin < 2 || isempty(sigmaHighFreq)
    sigmaHighFreq = 0.5;
end
if nargin < 3 || isempty(sigmaLowFreq)
    sigmaLowFreq = 0.1;
end

validateattributes(fRadius, {'numeric'}, ...
    {'scalar','real','finite','integer','>=',1}, mfilename, 'fRadius');
validateattributes(sigmaHighFreq, {'numeric'}, ...
    {'scalar','real','finite','positive'}, mfilename, 'sigmaHighFreq');
validateattributes(sigmaLowFreq, {'numeric'}, ...
    {'scalar','real','finite','positive'}, mfilename, 'sigmaLowFreq');
assert(sigmaHighFreq > sigmaLowFreq, ...
    'sigmaHighFreq 必须大于 sigmaLowFreq。');

sz = 2*fRadius + 1;
sampleScale = 32;                         % 与原版的频率采样密度一致
freq = linspace(-1, 1, sampleScale*sz);
[fx, fy] = meshgrid(freq, freq);
r2 = fx.^2 + fy.^2;

% 零直流带通目标：两高斯在 r=0 处的幅度相同，差值严格为零。
Htarget = exp(-r2/sigmaHighFreq^2) - exp(-r2/sigmaLowFreq^2);

% 有限尺寸 FIR 近似；保持与方案A一致的 fsamp2 构造方法。
filter = real(fsamp2(freq, freq, Htarget, [sz, sz]));

% 仅消除有限核的残余直流响应，不执行“中心减1”。
filter = filter - sum(filter(:))/numel(filter);

% 与方案A一样，正系数之和归一化为1，便于比较细节层的幅度。
posSum = sum(filter(filter > 0));
assert(isfinite(posSum) && posSum > eps, ...
    '带通核正系数之和无效，请调整 sigma 或核尺寸。');
filter = filter / posSum;
end
