clc;
close all;
clear all;

%% Load and show images
% Read initial image <-> template pair
img=imread('BaboonImage.png');
%tmp=imread('BaboonTemplate.png');
tmp=imread('BaboonTemplateIntensityChange.png');


% Plot both of these images
figure;imshow(img);title('Image','Fontsize',14);
figure;imshow(tmp);title('Template','Fontsize',14);

transform = 'euclidean';


% parameters for ECC and Lucas-Kanade 
par = [];
par.levels =    2;
par.iterations = 30;
par.transform = transform;

%% Lucas-Kanade algorithm
[LKWarp]=iat_LucasKanade(img,tmp,par);

% Compute the warped image and visualize the error
[wimageLK, supportLK] = iat_inverse_warping(img, LKWarp, par.transform, 1:size(tmp,2),1:size(tmp,1));

%plot the warped image
figure;imshow(uint8(wimageLK)); title('Warped image by Lucas-Kanade', 'Fontsize', 14);

% draw mosaic
LKMosaic = iat_mosaic(tmp,img,[LKWarp; 0 0 1]);
figure;imshow(uint8(LKMosaic));title('Mosaic after Lucas-Kanade','Fontsize',14);


%% ECC algorithm
[ECCWarp]=iat_ecc(img,tmp,par);

% Compute the warped image and visualize the error
[wimageECC, supportECC] = iat_inverse_warping(img, ECCWarp, par.transform, 1:size(tmp,2),1:size(tmp,1));

%plot the warped image
figure;imshow(uint8(wimageECC)); title('Warped image by ECC', 'Fontsize', 14);

% draw mosaic
ECCMosaic = iat_mosaic(tmp,img,[ECCWarp; 0 0 1]);
figure;imshow(uint8(ECCMosaic));title('Mosaic after ECC','Fontsize',14);

