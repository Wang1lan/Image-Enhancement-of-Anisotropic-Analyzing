function tests = test_gradientWeight_fx
%TEST_GRADIENTWEIGHT_FX 定点权重、独立整数黄金模型及浮点精度验收。
%   在 Fixed_GW 目录运行 runtests('tests/test_gradientWeight_fx.m')。
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.oldPath = path;
fixedDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(fixedDir, 'functions'));
addpath(fullfile(fileparts(fixedDir), 'functions'));
testCase.TestData.params = prepareGradientWeightFxParams( ...
    0.0024894783978878089, 0.0071887207526713029, 0.3);
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testDefaultParameters(testCase)
p = testCase.TestData.params;
verifyEqual(testCase, p.Tlow_Q4, uint32(1305));
verifyEqual(testCase, p.Thigh_Q4, uint32(3768));
verifyEqual(testCase, p.D_Q4, uint32(2463));
verifyEqual(testCase, p.Wmin12, uint16(1229));
verifyEqual(testCase, p.K_Q20, uint32(1220146));
verifyEqual(testCase, double(p.fracBitsGradient), 4);
verifyEqual(testCase, double(p.fracBitsSlope), 20);
verifyEqual(testCase, double(p.inputMax), 4095);
end

function testParameterQuantization(testCase)
ui = [0 1 0; 0 1 1; 0.1 0.9 0.5; 0 1/524160 0];
for k = 1:size(ui, 1)
    p = prepareGradientWeightFxParams(ui(k,1), ui(k,2), ui(k,3));
    low = round(ui(k,1)*524160);
    high = round(ui(k,2)*524160);
    minimum = round(ui(k,3)*4095);
    slope = round((4095-minimum)*1048576/(high-low));
    verifyEqual(testCase, p.Tlow_Q4, uint32(low));
    verifyEqual(testCase, p.Thigh_Q4, uint32(high));
    verifyEqual(testCase, p.Wmin12, uint16(minimum));
    verifyEqual(testCase, p.K_Q20, uint32(slope));
    verifyEqual(testCase, [p.GLow p.GHigh p.Wmin], ui(k,:));
end
p = prepareGradientWeightFxParams(0, 1/524160, 0);
verifyEqual(testCase, p.D_Q4, uint32(1));
verifyEqual(testCase, p.K_Q20, uint32(4293918720));
end

function testInvalidParameters(testCase)
invalid = { ...
    {[], .1, .3}, {[.01 .02], .1, .3}, {'0', .1, .3}, ...
    {NaN, .1, .3}, {0, Inf, .3}, {complex(0,1), .1, .3}, ...
    {-eps, .1, .3}, {0, 1+eps, .3}, {.1, .1, .3}, ...
    {.2, .1, .3}, {0, .1, -eps}, {0, .1, 1+eps}, ...
    {0, .1, []}, {0, .1, NaN}, {0, .1, Inf}, ...
    {0, 1e-10, .3}};
for k = 1:numel(invalid)
    args = invalid{k};
    verifyError(testCase, @() prepareGradientWeightFxParams(args{:}), ...
        'prepareGradientWeightFxParams:Parameters');
end
end

function testInvalidY12(testCase)
p = testCase.TestData.params;
invalid = {zeros(3), zeros(3,'uint8'), zeros(3,'int16'), ...
    zeros(2,3,'uint16'), zeros(3,3,2,'uint16'), ...
    zeros(0,0,'uint16'), repmat(uint16(4096),3,3)};
for k = 1:numel(invalid)
    image = invalid{k};
    verifyError(testCase, @() computeGradientWeight_fx(image,p), ...
        'computeGradientWeight_fx:Input');
end
end

function testConstantImages(testCase)
p = testCase.TestData.params;
for value = [0 73 2048 4095]
    image = repmat(uint16(value), 7, 11);
    [weight, dbg] = computeGradientWeight_fx(image, p);
    verifyEqual(testCase, dbg.Gx, zeros(size(image),'int32'));
    verifyEqual(testCase, dbg.Gy, zeros(size(image),'int32'));
    verifyEqual(testCase, dbg.G, zeros(size(image),'uint16'));
    verifyEqual(testCase, weight, repmat(p.Wmin12,size(image)));
    verifyAgainstGolden(testCase, image, p);
end
end

function testDirectionalStepsAndSigns(testCase)
p = testCase.TestData.params;
horizontal = zeros(11,15,'uint16');
horizontal(:,8:end) = uint16(4095);
[~, dx] = computeGradientWeight_fx(horizontal, p);
verifyEqual(testCase, dx.Gx(:,7:8), repmat(int32(-16380),11,2));
verifyEqual(testCase, dx.Gy, zeros(size(horizontal),'int32'));
verifyEqual(testCase, dx.G(:,7:8), repmat(uint16(16380),11,2));
vertical = horizontal';
[~, dy] = computeGradientWeight_fx(vertical, p);
verifyEqual(testCase, dy.Gy(7:8,:), repmat(int32(-16380),2,11));
verifyEqual(testCase, dy.Gx, zeros(size(vertical),'int32'));
verifyAgainstGolden(testCase, horizontal, p);
verifyAgainstGolden(testCase, vertical, p);
end

function testSyntheticImagesAgainstGolden(testCase)
images = syntheticImages;
for k = 1:numel(images)
    verifyAgainstGolden(testCase, images{k}, testCase.TestData.params);
end
end

function testSmallAndNarrowImages(testCase)
shapes = [3 3; 3 13; 14 3; 4 5; 5 4; 9 9];
for k = 1:size(shapes,1)
    count = prod(shapes(k,:));
    image = uint16(reshape(mod((0:count-1)*887,4096),shapes(k,:)));
    verifyAgainstGolden(testCase, image, testCase.TestData.params);
end
end

function testWeightEndpoints(testCase)
image = uint16(reshape(mod((0:142)*137,4096),11,13));
for minimum = [0 1]
    p = prepareGradientWeightFxParams(.0024894783978878089, ...
        .0071887207526713029, minimum);
    [weight, dbg] = computeGradientWeight_fx(image, p);
    verifyAgainstGolden(testCase, image, p);
    if minimum == 1
        verifyEqual(testCase, weight, repmat(uint16(4095),size(image)));
        verifyEqual(testCase, dbg.W0_12, weight);
        verifyEqual(testCase, p.K_Q20, uint32(0));
    end
end
end

function testExtremeThresholdConfigurations(testCase)
ui = [0 1 0; 0 1/524160 0; 1-1/524160 1 0; ...
    0 1 1; .2 .6 .9999];
image = uint16(reshape(mod((0:98)*911,4096),9,11));
for k = 1:size(ui,1)
    p = prepareGradientWeightFxParams(ui(k,1),ui(k,2),ui(k,3));
    verifyAgainstGolden(testCase, image, p);
end
end

function testThresholdBranchesAndClipping(testCase)
p = prepareGradientWeightFxParams(10/32760,14/32760,.3);
gradient = uint16([0 9 10 11 12 13 14 15 24570]);
[actual, q4] = mapGradientWeight_fx(gradient,p);
verifyEqual(testCase, q4, uint32(gradient)*uint32(16));
verifyEqual(testCase, actual, uint16([1229 1229 1229 1946 2662 3379 4095 4095 4095]));
verifyEqual(testCase, actual, uint16(goldenMap(int64(gradient),p)));
p = testCase.TestData.params;
gradient = uint16([80 81 82 158 234 235 236 237]);
verifyEqual(testCase, mapGradientWeight_fx(gradient,p), ...
    uint16(goldenMap(int64(gradient),p)));
end

function testExhaustiveGradientMapping(testCase)
p = testCase.TestData.params;
gradient = uint16(0:24570);
[actual, q4] = mapGradientWeight_fx(gradient,p);
verifyEqual(testCase, actual, uint16(goldenMap(int64(gradient),p)));
verifyEqual(testCase, q4, uint32(gradient)*uint32(16));
normalized = double(gradient)/32760;
reference = p.Wmin+(1-p.Wmin)*min(max( ...
    (normalized-p.GLow)/(p.GHigh-p.GLow),0),1);
error = abs(double(actual)-round(reference*4095));
verifyLessThanOrEqual(testCase, max(error), 1);
verifyGreaterThanOrEqual(testCase, diff(int32(actual)), int32(0));
verifyEqual(testCase, actual([1 end]), [p.Wmin12 uint16(4095)]);
end

function testExhaustiveDivision81(testCase)
sum81 = uint32(0:331695);
actual = divide81Round_fx(sum81);
expected = uint16(idivide(uint64(sum81)+uint64(40),uint64(81),'floor'));
verifyEqual(testCase, actual, expected);
verifyClass(testCase, actual, 'uint16');
verifyEqual(testCase, actual([1 end]), uint16([0 4095]));
product = (uint64(sum81)+uint64(40)).*uint64(828505); % u38
verifyEqual(testCase, max(product), uint64(274844106175));
verifyLessThanOrEqual(testCase, product, bitshift(uint64(1),38)-uint64(1));
end

function testSobelBinaryVertices(testCase)
% 每个3x3块的中心保留独立邻域，一帧覆盖512个二值极点。
image = zeros(3,3*512,'uint16');
expectedX = zeros(1,512,'int32');
expectedY = zeros(1,512,'int32');
kernelX = int64([1 0 -1; 2 0 -2; 1 0 -1]);
for vertex = 0:511
    patch = reshape(bitget(uint16(vertex),1:9),3,3).*uint16(4095);
    image(:,3*vertex+(1:3)) = patch;
    x = int64(0);
    y = int64(0);
    for row = 1:3
        for col = 1:3
            x = x+int64(patch(row,col))*kernelX(row,col); % s15
            y = y+int64(patch(row,col))*kernelX(col,row); % s15
        end
    end
    expectedX(vertex+1) = int32(x);
    expectedY(vertex+1) = int32(y);
end
[~, dbg] = computeGradientWeight_fx(image,testCase.TestData.params);
centers = 2:3:size(image,2);
verifyEqual(testCase, dbg.Gx(2,centers), expectedX);
verifyEqual(testCase, dbg.Gy(2,centers), expectedY);
verifyEqual(testCase, dbg.G(2,centers), ...
    uint16(abs(expectedX)+abs(expectedY)));
verifyEqual(testCase, max(dbg.G(2,centers)), uint16(24570));
verifyEqual(testCase, max(expectedX), int32(16380));
verifyEqual(testCase, min(expectedX), int32(-16380));
verifyGreaterThanOrEqual(testCase, expectedX, int32(-16384));
verifyLessThanOrEqual(testCase, expectedX, int32(16383));
verifyGreaterThanOrEqual(testCase, expectedY, int32(-16384));
verifyLessThanOrEqual(testCase, expectedY, int32(16383));
end

function testFloatStageAgreement(testCase)
p = testCase.TestData.params;
images = syntheticImages;
for k = 1:numel(images)
    [weight, fixed] = computeGradientWeight_fx(images{k},p);
    [reference, floating] = computeGradientWeight_floatRef( ...
        images{k},p.GLow,p.GHigh,p.Wmin);
    verifyEqual(testCase, fixed.Gx, int32(round(floating.Gx*32760)));
    verifyEqual(testCase, fixed.Gy, int32(round(floating.Gy*32760)));
    verifyEqual(testCase, fixed.G, uint16(round(floating.G*32760)));
    initialError = abs(double(fixed.W0_12)-round(floating.W0*4095));
    finalError = abs(double(weight)-round(reference*4095));
    verifyLessThanOrEqual(testCase, max(initialError(:)),1);
    verifyLessThanOrEqual(testCase, max(finalError(:)),2);
    verifyEqual(testCase, reference, floating.W_smooth);
    verifyGreaterThanOrEqual(testCase, reference, p.Wmin-1e-12);
    verifyLessThanOrEqual(testCase, reference, 1+1e-12);
end
end

function testOutputAndIntermediateRanges(testCase)
p = testCase.TestData.params;
images = syntheticImages;
seenMiddle = false;
for k = 1:numel(images)
    [weight, d] = computeGradientWeight_fx(images{k},p);
    verifyClass(testCase, weight,'uint16');
    verifySize(testCase, weight,size(images{k}));
    verifyClass(testCase, d.Gx,'int32');
    verifyClass(testCase, d.Gy,'int32');
    verifyClass(testCase, d.G,'uint16');
    verifyClass(testCase, d.G_Q4,'uint32');
    verifyClass(testCase, d.W0_12,'uint16');
    verifyClass(testCase, d.S9,'uint32');
    verifyClass(testCase, d.S81,'uint32');
    verifyGreaterThanOrEqual(testCase, d.Gx,int32(-16380));
    verifyLessThanOrEqual(testCase, d.Gx,int32(16380));
    verifyGreaterThanOrEqual(testCase, d.Gy,int32(-16380));
    verifyLessThanOrEqual(testCase, d.Gy,int32(16380));
    verifyLessThanOrEqual(testCase, abs(d.Gx),int32(16380)); % u14
    verifyLessThanOrEqual(testCase, abs(d.Gy),int32(16380)); % u14
    verifyLessThanOrEqual(testCase, d.G,uint16(24570));      % u15
    verifyLessThanOrEqual(testCase, d.G_Q4,uint32(393120)); % u19
    verifyLessThanOrEqual(testCase, d.S9,uint32(36855));    % u16
    verifyLessThanOrEqual(testCase, d.S81,uint32(331695));  % u19
    verifyGreaterThanOrEqual(testCase, d.W0_12,p.Wmin12);
    verifyLessThanOrEqual(testCase, d.W0_12,uint16(4095));  % u12
    verifyGreaterThanOrEqual(testCase, weight,p.Wmin12);
    verifyLessThanOrEqual(testCase, weight,uint16(4095));   % u12
    middle = d.G_Q4>p.Tlow_Q4 & d.G_Q4<p.Thigh_Q4;
    seenMiddle = seenMiddle || any(middle(:));
    delta = uint64(d.G_Q4(middle))-uint64(p.Tlow_Q4);       % u19
    product = delta.*uint64(p.K_Q20);                     % u33 reserved
    roundedNumerator = product+uint64(524288);            % u33 reserved
    verifyTrue(testCase, all(roundedNumerator(:)<=uint64(4294705088)));
    verifyEqual(testCase, d.W12,weight);
end
verifyTrue(testCase, seenMiddle, '图案组必须覆盖线性映射中间段。');
end

function testOutputOnlyEntry(testCase)
image = uint16(reshape(0:62,7,9));
p = testCase.TestData.params;
actual = computeGradientWeight_fx(image,p);
[expected, ~] = computeGradientWeight_fx(image,p);
verifyEqual(testCase, actual,expected);
end

function images = syntheticImages
impulse = zeros(11,13,'uint16');
impulse(6,7) = uint16(4095);
[cols,rows] = meshgrid(0:12,0:10);
checker = uint16(mod(rows+cols,2)*4095);
ramp = uint16(rows*149+cols*173);
stream = RandStream('mt19937ar','Seed',20261112);
random = uint16(randi(stream,[0 4095],11,13));
lowContrast = uint16(randi(stream,[0 75],11,13)+1000);
images = {impulse,checker,ramp,random,lowContrast};
end

function verifyAgainstGolden(testCase,image,params)
[weight, actual] = computeGradientWeight_fx(image,params);
expected = goldenModel(image,params);
names = fieldnames(expected);
for k = 1:numel(names)
    verifyEqual(testCase, actual.(names{k}),expected.(names{k}));
end
verifyEqual(testCase, weight,expected.W12);
end

function result = goldenModel(image,params)
% 独立逐点黄金模型，所有数据运算用int64，按定义直接求81项和。
[height,width] = size(image);
kernelX = int64([1 0 -1; 2 0 -2; 1 0 -1]);
gx = zeros(height,width,'int64');
gy = zeros(height,width,'int64');
for row = 1:height
    for col = 1:width
        for dy = -1:1
            inputRow = max(1,min(height,row+dy));
            for dx = -1:1
                inputCol = max(1,min(width,col+dx));
                pixel = int64(image(inputRow,inputCol)); % u12
                termX = pixel*kernelX(dy+2,dx+2);         % s14
                termY = pixel*kernelX(dx+2,dy+2);         % s14
                gx(row,col) = gx(row,col)+termX;          % s15
                gy(row,col) = gy(row,col)+termY;          % s15
            end
        end
    end
end
gradient = abs(gx)+abs(gy);                               % u15
q4 = gradient*int64(16);                                 % u19
initial = goldenMap(gradient,params);                    % u12
sum9 = zeros(height,width,'int64');
sum81 = zeros(height,width,'int64');
for row = 1:height
    for col = 1:width
        for dx = -4:4
            inputCol = max(1,min(width,col+dx));
            sum9(row,col) = sum9(row,col)+initial(row,inputCol); % u16
        end
        for dy = -4:4
            inputRow = max(1,min(height,row+dy));
            for dx = -4:4
                inputCol = max(1,min(width,col+dx));
                sum81(row,col) = sum81(row,col)+initial(inputRow,inputCol); % u19
            end
        end
    end
end
weight = idivide(sum81+int64(40),int64(81),'floor');       % u12
weight = min(max(weight,int64(params.Wmin12)),int64(4095)); % u12
result = struct('Gx',int32(gx),'Gy',int32(gy), ...
    'G',uint16(gradient),'G_Q4',uint32(q4), ...
    'W0_12',uint16(initial),'S9',uint32(sum9), ...
    'S81',uint32(sum81),'W12',uint16(weight));
end

function weight = goldenMap(gradient,params)
% 独立精确整数映射，逐点比较和idivide，不复用正式映射/移位函数。
q4 = gradient*int64(16);                                 % u19
minimum = int64(params.Wmin12);                          % u12
weight = repmat(minimum,size(gradient));                 % u12
for k = 1:numel(gradient)
    if q4(k)>=int64(params.Thigh_Q4)
        weight(k) = int64(4095);                         % u12
    elseif q4(k)>int64(params.Tlow_Q4)
        delta = q4(k)-int64(params.Tlow_Q4);              % u19
        product = delta*int64(params.K_Q20);             % u33 reserved
        numerator = product+int64(524288);               % u33 reserved
        increment = idivide(numerator,int64(1048576),'floor'); % u13
        estimate = minimum+increment;                   % u13
        weight(k) = min(max(estimate,minimum),int64(4095)); % u12
    end
end
end
