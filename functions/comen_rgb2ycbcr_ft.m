function ImgYcbcr = comen_rgb2ycbcr_ft(imgin)

imgin = double(imgin);
[M,N,D] = size(imgin);
ImgYcbcr = zeros(M,N,D);

para = [ 0.2126   0.7152   0.0722
        -0.1146  -0.3854   0.5000
         0.5000  -0.4542  -0.0458];
    
R = imgin(:,:,1);
G = imgin(:,:,2);
B = imgin(:,:,3);

Y  = para(1,1) * R + para(1,2) * G + para(1,3) * B;
Cb = para(2,1) * R + para(2,2) * G + para(2,3) * B;
Cr = para(3,1) * R + para(3,2) * G + para(3,3) * B;

ImgYcbcr(:,:,1) = Y;
ImgYcbcr(:,:,2) = Cb;
ImgYcbcr(:,:,3) = Cr;