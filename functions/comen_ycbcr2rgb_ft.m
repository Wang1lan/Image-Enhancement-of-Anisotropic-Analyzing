function ImgRGB = comen_ycbcr2rgb_ft(ImgYcbcr)

ImgYcbcr = double(ImgYcbcr);

para = [1.0000   0.0000    1.5748
        1.0000   -0.1873   -0.4681
        1.0000    1.8556    0.0000];
    
    
Y = ImgYcbcr(:,:,1);
Cb = ImgYcbcr(:,:,2);
Cr = ImgYcbcr(:,:,3);

r = para(1,1) * Y + para(1,2) * Cb + para(1,3) * Cr;
g = para(2,1) * Y + para(2,2) * Cb + para(2,3) * Cr;
b = para(3,1) * Y + para(3,2) * Cb + para(3,3) * Cr;

ImgRGB(:,:,1) = r;
ImgRGB(:,:,2) = g;
ImgRGB(:,:,3) = b;

ImgRGB(ImgRGB < 0) = 0;
ImgRGB(ImgRGB > 1 - 2^(-13)) = 1 - 2^(-13);