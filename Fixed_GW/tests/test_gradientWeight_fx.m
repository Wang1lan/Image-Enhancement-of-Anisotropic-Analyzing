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

function testComenColorConversion(testCase)
% 每行一个 RGB8 色块，顺序为黑、白、中灰、红、绿、蓝。
rgb8 = uint8([0 0 0; 255 255 255; 128 128 128; 255 0 0; 0 255 0; 0 0 255]);
patches = reshape(double(rgb8)/255,6,1,3);
ycc = comen_rgb2ycbcr_ft(patches);
expected = [0 0 0; 1 0 0; 128/255 0 0; ...
    .2126 -.1146 .5000; .7152 -.3854 -.4542; .0722 .5000 -.0458];
verifyEqual(testCase,reshape(ycc,6,3),expected,'AbsTol',2*eps);
verifyEqual(testCase,uint16(round(4095*ycc(:,:,1))),uint16([0;4095;2056;871;2929;296]));
rgb = comen_ycbcr2rgb_ft(ycc);
verifyEqual(testCase,uint8(round(255*reshape(rgb,6,3))),rgb8);
neutral = zeros(2,1,3);
neutral(:,1,1) = [-.1;1.1];
verifyEqual(testCase,reshape(comen_ycbcr2rgb_ft(neutral),2,3),[zeros(1,3);(1-2^-13)*ones(1,3)]);
end

function testQ17RegisterVectors(testCase)
ui = [.0024894783978878089 .0071887207526713029; ...
      .0049333333333332879 .0090784313725489965];
words = [5220 15072 38130; 10343 19034 43223];
for k = 1:2
    p = prepareGradientWeightFxParams(ui(k,1),ui(k,2),.3);
    regs = loadGradientWeightFxRegs(1229,words(k,1),words(k,2),words(k,3));
    verifyEqual(testCase,rmfield(p,{'GLow','GHigh','Wmin'}),regs);
    verifyEqual(testCase,[p.GLow p.GHigh p.Wmin],[ui(k,:) .3]);
    verifyEqual(testCase,regs.Wmin12,uint16(1229));
    verifyEqual(testCase,regs.D_Q6,uint32(words(k,2)-words(k,1)));
    verifyEqual(testCase,regs.K_Q17,uint32(words(k,3)));
end
p = testCase.TestData.params;
verifyEqual(testCase, p.Format, "Q17");
verifyEqual(testCase, [p.fracBitsGradient p.fracBitsSlope p.inputMax], [6 17 4095]);
gradient = uint16([0 81 82 102 235 236 4096 24570]);
expected = uint16([1229 1229 1237 1610 4086 4095 4095 4095]);
verifyEqual(testCase,mapGradientWeight_fx(gradient,p),expected);
end

function testParameterQuantization(testCase)
low = 2000.53/2096640;
p = prepareGradientWeightFxParams(low,low+.003,.1);
verifyEqual(testCase,[p.Tlow_Q6 p.Thigh_Q6 p.D_Q6 p.K_Q17],uint32([2001 8290 6289 76801]));
verifyEqual(testCase,p.Wmin12,uint16(410));
p = prepareGradientWeightFxParams(0,.1,.5);
verifyEqual(testCase,p.Wmin12,uint16(2048));
verifyEqual(testCase,[p.Tlow_Q6 p.Thigh_Q6 p.K_Q17],uint32([0 209664 1280]));
p = prepareGradientWeightFxParams(.097,.1,.3);
verifyEqual(testCase,p.Tlow_Q6,uint32(203374));
p = prepareGradientWeightFxParams(.081,.084,.3);
verifyGreaterThanOrEqual(testCase,p.D_Q6,uint32(6289));
end

function testInvalidRegisters(testCase)
cases = {{1228.9,5220,15072,38130}, {409,5220,15072,38130}, ...
    {2049,5220,15072,38130}, {1229,-1,15072,38130}, {1229,5220,209665,38130}, ...
    {1229,NaN,15072,38130}, {1229,5220,6000,38130}, {1229,15072,5220,38130}, ...
    {1229,5220,15072,76802}, {1229,5220,15072,1279}};
for k = 1:numel(cases)
    args = cases{k};
    verifyError(testCase,@()loadGradientWeightFxRegs(args{:}),'loadGradientWeightFxRegs:Range');
end
mismatch = {{1229,5220,15072,38129}, {1229,5220,15072,12958}};
for k = 1:numel(mismatch)
    args = mismatch{k};
    verifyError(testCase,@()loadGradientWeightFxRegs(args{:}),'loadGradientWeightFxRegs:Mismatch');
end
end

function testInvalidParameters(testCase)
invalid = { ...
    {[], .1, .3}, {[.01 .02], .1, .3}, {'0', .1, .3}, ...
    {NaN, .1, .3}, {0, Inf, .3}, {complex(0,1), .1, .3}, ...
    {-eps, .1, .3}, {0, 1+eps, .3}, {.1, .1, .3}, ...
    {.2, .1, .3}, {0, .1, -eps}, {0, .1, 1+eps}, ...
    {0, .1, []}, {0, .1, NaN}, {0, .1, Inf}, ...
    {0, 1e-10, .3}, {.081, .084-1e-12, .3}, {0,.10000001,.3}, ...
    {0,.01,.1-1e-10}, {0,.01,.5+1e-10}};
for k = 1:numel(invalid)
    args = invalid{k};
    verifyError(testCase, @() prepareGradientWeightFxParams(args{:}), ...
        'prepareGradientWeightFxParams:Parameters');
end
end

function testFloatReferenceRejectsInvalidCalibration(testCase)
image = zeros(3,5,'uint16');
verifyError(testCase,@()computeGradientWeight_floatRef(image,.08,.083-1e-12,.3), ...
    'computeGradientWeight_floatRef:Parameters');
verifyError(testCase,@()computeGradientWeight_floatRef(image,.07,.09,.05), ...
    'computeGradientWeight_floatRef:Parameters');
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
    p = testCase.TestData.params;
    verifyAgainstGolden(testCase,images{k},p);
end
end

function testSmallAndNarrowImages(testCase)
shapes = [3 3; 3 13; 14 3; 4 5; 5 4; 9 9];
for k = 1:size(shapes,1)
    count = prod(shapes(k,:));
    image = uint16(reshape(mod((0:count-1)*887,4096),shapes(k,:)));
    p = testCase.TestData.params;
    verifyAgainstGolden(testCase,image,p);
    [w,df] = computeGradientWeight_floatRef(image,p.GLow,p.GHigh,p.Wmin);
    verifyEqual(testCase,w,imfilter(df.W0,ones(7)/49,'replicate','same'),'AbsTol',1e-12);
end
end

function testFrameConfigurations(testCase)
image = zeros(9,13,'uint16');
image(5,7) = 51;  % G=102 使两组 Q17 标定产生不同 W0。
q = testCase.TestData.params;
other = prepareGradientWeightFxParams(0,.1,.5);
[~,first] = computeGradientWeight_fx(image,q);
[~,second] = computeGradientWeight_fx(image,other);
verifyTrue(testCase,any(first.W0_12(:) ~= second.W0_12(:)));
for config = {q,other,q}
    verifyAgainstGolden(testCase,image,config{1});
end
verifyError(testCase,@()loadGradientWeightFxRegs(1229,5220,15072,43223), ...
    'loadGradientWeightFxRegs:Mismatch');
[~,last] = computeGradientWeight_fx(image,q);
verifyEqual(testCase,last,first);
end

function testThresholdBranchesAndClipping(testCase)
p = prepareGradientWeightFxParams(10/32760,150/32760,.3);
g = uint16([0 9 10 11 149 150 151 4096 24570]);
[actual,x] = mapGradientWeight_fx(g,p);
verifyEqual(testCase,x,uint32(g).*uint32(64));
verifyEqual(testCase,x([3 6]),[p.Tlow_Q6 p.Thigh_Q6]);
verifyEqual(testCase,actual([1:3 6:9]),uint16([1229 1229 1229 4095 4095 4095 4095]));
verifyEqual(testCase,actual,uint16(goldenMap(int64(g),p)));
p = testCase.TestData.params;
middleX = uint64(235*64)-uint64(p.Tlow_Q6);
maxDefaultNumerator = middleX.*uint64(p.K_Q17)+uint64(65536);
verifyEqual(testCase,maxDefaultNumerator,uint64(374502136));
verifyGreaterThan(testCase,maxDefaultNumerator,bitshift(uint64(1),28)-uint64(1));
end

function testExhaustiveGradientMapping(testCase)
low = 2000.53/2096640;
ui = [.0024894783978878089 .0071887207526713029 .3; ...
    .0049333333333332879 .0090784313725489965 .3; 0 .1 .1; 0 .1 .5; ...
    .097 .1 .5; low low+.003 .1; .00179 .00479 .1; 0 .09 .25; ...
    .0007 .004 .1; 10/32760 150/32760 .3; 0 .02 .3];
stream = RandStream('mt19937ar','Seed',20261118);
lows = .097*rand(stream,128,1);
ui = [ui; lows, lows+.003+(.097-lows).*rand(stream,128,1), .1+.4*rand(stream,128,1)];
g = uint16(0:24570);
for k = 1:size(ui,1)
    low = ui(k,1); high = ui(k,2); minimum = ui(k,3);
    reference = round(4095*(minimum+(1-minimum)*min(max( ...
        (double(g)/32760-low)/(high-low),0),1)));
    p = prepareGradientWeightFxParams(low,high,minimum);
    expectedK = round((4095-double(p.Wmin12))*2^17/double(p.D_Q6));
    verifyEqual(testCase,p.K_Q17,uint32(expectedK));
    [actual,x] = mapGradientWeight_fx(g,p);
    verifyEqual(testCase,actual,uint16(goldenMap(int64(g),p)));
    verifyTrue(testCase,all(diff(int32(actual))>=0));
    verifyLessThanOrEqual(testCase,max(abs(double(actual)-reference)),1);
    verifyTrue(testCase,all(actual(x<=p.Tlow_Q6)==p.Wmin12));
    verifyTrue(testCase,all(actual(x>=p.Thigh_Q6)==4095));
end
end

function testNarrowFormatRegressions(testCase)
ui = [.00179 .00479 .1; 0 .09 .25; .0007 .004 .1];
g = [60 2800 130];
floatCodes = [460 3940 4060];
for k = 1:3
    low = ui(k,1); high = ui(k,2); w = ui(k,3);
    verifyEqual(testCase,round(4095*(w+(1-w)*(g(k)/32760-low)/(high-low))),floatCodes(k));
    p = prepareGradientWeightFxParams(low,high,w);
    verifyLessThanOrEqual(testCase,abs(double(mapGradientWeight_fx(uint16(g(k)),p))-floatCodes(k)),1);
end
% 过窄格式重现同一组输入的两码偏差。
for bits = [4 5]
    ends = round(32760*2^bits*ui(1,1:2));
    verifyEqual(testCase,410+round((60*2^bits-ends(1))*3685/diff(ends)),462);
end
for k = 2:3
    shifts = [16 12]; expected = [3942 4058];
    ends = round(2096640*ui(k,1:2));
    w = round(4095*ui(k,3));
    s = shifts(k-1);
    c = round((4095-w)*2^s/diff(ends));
    verifyEqual(testCase,w+round((64*g(k)-ends(1))*c/2^s),expected(k-1));
end
end

function testExhaustiveDivision49(testCase)
totals = uint32(0:200655);
actual = divide49Round_fx(totals);
x = uint64(totals)+uint64(24);
exact = idivide(x,uint64(49),'floor');
product = x.*uint64(2675);
verifyEqual(testCase,actual,uint16(idivide(product,uint64(131072),'floor')));
verifyEqual(testCase,max(product),uint64(536816325));
delta = int32(actual)-int32(exact);
verifyTrue(testCase,all(delta==0 | delta==1));
verifyEqual(testCase,nnz(delta==1),7466);
verifyLessThanOrEqual(testCase,max(actual),uint16(4095));
verifyEqual(testCase,actual([1 26 43684 end]),uint16([0 1 892 4095]));
verifyEqual(testCase,exact(43684),uint64(891));          % S49=43683
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
    verifyEqual(testCase,reference,imfilter(floating.W0,ones(7)/49,'replicate','same'),'AbsTol',1e-12);
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
    verifyClass(testCase, d.G_Q6,'uint32');
    verifyClass(testCase, d.W0_12,'uint16');
    verifyClass(testCase, d.S7,'uint32');
    verifyClass(testCase, d.S49,'uint32');
    verifyGreaterThanOrEqual(testCase, d.Gx,int32(-16380));
    verifyLessThanOrEqual(testCase, d.Gx,int32(16380));
    verifyGreaterThanOrEqual(testCase, d.Gy,int32(-16380));
    verifyLessThanOrEqual(testCase, d.Gy,int32(16380));
    verifyLessThanOrEqual(testCase, abs(d.Gx),int32(16380)); % u14
    verifyLessThanOrEqual(testCase, abs(d.Gy),int32(16380)); % u14
    verifyLessThanOrEqual(testCase, d.G,uint16(24570));      % u15
    verifyLessThanOrEqual(testCase, d.G_Q6,uint32(1572480)); % u21
    verifyLessThanOrEqual(testCase, d.S7,uint32(28665));    % u15
    verifyLessThanOrEqual(testCase, d.S49,uint32(200655));  % u18
    verifyGreaterThanOrEqual(testCase, d.W0_12,p.Wmin12);
    verifyLessThanOrEqual(testCase, d.W0_12,uint16(4095));  % u12
    verifyGreaterThanOrEqual(testCase, weight,p.Wmin12);
    verifyLessThanOrEqual(testCase, weight,uint16(4095));   % u12
    middle = d.G_Q6>p.Tlow_Q6 & d.G_Q6<p.Thigh_Q6;
    seenMiddle = seenMiddle || any(middle(:));
    delta = uint64(d.G_Q6(middle))-uint64(p.Tlow_Q6); % u18
    numerator = delta.*uint64(p.K_Q17)+uint64(65536);     % u29
    verifyTrue(testCase,all(numerator(:)<=uint64(483170688)));
    verifyEqual(testCase, d.W12,weight);
end
verifyTrue(testCase, seenMiddle, '图案组必须覆盖线性映射中间段。');
end

function testOutlierCoordinatesAndDecomposition(testCase)
p = testCase.TestData.params;
image = zeros(9,11,'uint16'); image(5,6) = 51;
[~,fx] = computeGradientWeight_fx(image,p);
[~,df] = computeGradientWeight_floatRef(image,p.GLow,p.GHigh,p.Wmin);
idx = find(fx.G==102,1);
fx.W0_12(idx) = fx.W0_12(idx)+uint16(4);
fx.W12(1) = fx.W12(1)+uint16(5);
rows = gradientWeightAccuracyOutliers(fx,df,p);
verifyEqual(testCase,sort(rows.Stage),["W0";"W7"]);
mid = rows(rows.Stage=="W0",:);
verifyEqual(testCase,sub2ind(size(image),mid.PixelRow,mid.PixelColumn),idx);
verifyEqual(testCase,mid.Branch,"middle");
verifyEqual(testCase,[mid.Delta mid.C mid.s],[1308 38130 17]);
verifyEqual(testCase,mid.S49,fx.S49(idx));
verifyEqual(testCase,mid.EndpointError+mid.SlopeError+mid.PixelRoundError, ...
    double(fx.W0_12(idx))-4095*df.W0(idx),'AbsTol',1e-12);
edge = rows(rows.Stage=="W7",:);
verifyEqual(testCase,[edge.PixelRow edge.PixelColumn],[1 1]);
verifyEqual(testCase,edge.Branch,"low");
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
reference = computeGradientWeight_floatRef(image,params.GLow,params.GHigh,params.Wmin);
errorCodes = double(weight)-round(4095*reference);
verifyLessThanOrEqual(testCase,max(abs(errorCodes),[],'all'),2);
end

function result = goldenModel(image,params)
% 独立逐点黄金模型，所有数据运算用int64，按定义直接求49项和。
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
q6 = gradient*int64(64);                                % u21
initial = goldenMap(gradient,params);
sum7 = zeros(height,width,'int64');
sum49 = zeros(height,width,'int64');
for row = 1:height
    for col = 1:width
        for dx = -3:3
            inputCol = max(1,min(width,col+dx));
            sum7(row,col) = sum7(row,col)+initial(row,inputCol);
        end
        for dy = -3:3
            inputRow = max(1,min(height,row+dy));
            for dx = -3:3
                inputCol = max(1,min(width,col+dx));
                sum49(row,col) = sum49(row,col)+initial(inputRow,inputCol);
            end
        end
    end
end
weight = idivide((uint64(sum49)+uint64(24)).*uint64(2675),uint64(131072),'floor');
result = struct('Gx',int32(gx),'Gy',int32(gy),'G',uint16(gradient),'G_Q6',uint32(q6), ...
    'W0_12',uint16(initial),'S7',uint32(sum7),'S49',uint32(sum49),'W12',uint16(weight));
end

function weight = goldenMap(gradient,params)
% 用 Q17 寄存器值和整数除法对照像素端乘加、移位及饱和。
coefficient = int64(params.K_Q17);
x = gradient.*int64(64);
minimum = int64(params.Wmin12);
delta = max(x-int64(params.Tlow_Q6),int64(0));
increment = idivide(delta.*coefficient+int64(65536),int64(131072),'floor');
weight = min(minimum+increment,int64(4095));
weight(x>=int64(params.Thigh_Q6)) = int64(4095);
end
