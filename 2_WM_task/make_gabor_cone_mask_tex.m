function [imgout1,imgout2,imgout3,imgout4] = make_gabor_cone_mask_tex(ppd,res,g_r)

% ppd      =  pixels per d.v.a
% res      =  resolution of display [horizontal, vertical]
% g_r      =  radius (d.v.a.) of gabors

dva_res = res./2./ppd;   % screen resolution expressed as d.v.a.

% make transparent mask for 4 possible quadrants of centrally-presented gabor (i.e. memoranda)
[x, y] = meshgrid(linspace(-dva_res(1), dva_res(1), res(1)), linspace(-dva_res(2), dva_res(2), res(2)));  % create matrix of X and Y d.v.a. values for each pixel
r = sqrt(y.^2 + x.^2);  % calculate polar radius coordinate for each pixel

mask1 = ones(size(r));
mask1(r<g_r & x<0 & y<0) = 0;   % mask for TOP-LEFT quadrant

mask2 = ones(size(r));
mask2(r<g_r & x>0 & y<0) = 0;   % mask for TOP-RIGHT quadrant

mask3 = ones(size(r));
mask3(r<g_r & x<0 & y>0) = 0;   % mask for BOTTOM-LEFT quadrant

mask4 = ones(size(r));
mask4(r<g_r & x>0 & y>0) = 0;   % mask for BOTTOM-RIGHT quadrant

unigrey = ones(size(r)).*127;  % initializing pure grey image
alphaout = ones(size(r)).*255;   % intializing matrix of transparency values
alphaout(mask1==0) = 0;  % setting everything within this gabor quadrant to be completely transparent (0 = fully transparent, 255 = full opaque)
imgout1 = cat(3,unigrey,unigrey,unigrey,alphaout);

alphaout = ones(size(r)).*255;   % intializing matrix of transparency values
alphaout(mask2==0) = 0;  % setting everything within this gabor quadrant to be completely transparent (0 = fully transparent, 255 = full opaque)
imgout2 = cat(3,unigrey,unigrey,unigrey,alphaout);

alphaout = ones(size(r)).*255;   % intializing matrix of transparency values
alphaout(mask3==0) = 0;  % setting everything within this gabor quadrant to be completely transparent (0 = fully transparent, 255 = full opaque)
imgout3 = cat(3,unigrey,unigrey,unigrey,alphaout);

alphaout = ones(size(r)).*255;   % intializing matrix of transparency values
alphaout(mask4==0) = 0;  % setting everything within this gabor quadrant to be completely transparent (0 = fully transparent, 255 = full opaque)
imgout4 = cat(3,unigrey,unigrey,unigrey,alphaout);