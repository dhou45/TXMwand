sample_idx = 'R1_025';
% plot mask for 95% rsquare
	img_test = fit_gauss_adjrsquare_3d_resize4;
	sizeR = size(img_test);
	noise_level = 0.95;
	img_test_bw = img_test;
	img_test_bw(img_test_bw>noise_level)=1;img_test_bw(img_test_bw<=noise_level)=0;
	img_test_bw(~isfinite(img_test_bw))=0;
	fig = figure(9); orthosliceViewer(img_test_bw);
	print(gcf,['noise_map_resize4_', sample_idx, '_',num2str(noise_level), '.png'],'-dpng', '-r300');pause(0.5)


% segmentation, use peak hight as mask
	sizeR = size(fit_gauss_intensity_3d_resize4);
	img_test = fit_gauss_intensity_3d_resize4;
	figure; orthosliceViewer(img_test); colormap jet; pause(0.5)

	normalizedImage = uint8(255*mat2gray(img_test));
	Binary_Mn = imbinarize(normalizedImage);
	figure; volshow(Binary_Mn); pause(0.5)

	SE=strel('cube',1);
	Binary_Mn_close=imclose(Binary_Mn,SE);
	Binary_Mn_fill = imfill(Binary_Mn_close,'holes');
	CC = bwconncomp(Binary_Mn_fill, 6);
	numPixels = cellfun(@numel,CC.PixelIdxList);
	[~,idx] = max(numPixels);
	filtered_vol = false(size(Binary_Mn_fill));
	filtered_vol(CC.PixelIdxList{idx}) = true;
	figure; volshow(filtered_vol); pause(0.5)
	figure; orthosliceViewer(img_test .* filtered_vol); colormap jet; pause(0.5)

	% erode x pixel, then choose largest, then dilate back (change erodeSize to smooth the mask, but no inner gaps)
	erodeSize = 5;
	dilateSize = erodeSize; % consider if dilate or not.

	SE=strel('cube',erodeSize);
	filtered_vol_erode = imerode(filtered_vol,SE);
	CC = bwconncomp(filtered_vol_erode, 6);
	numPixels = cellfun(@numel,CC.PixelIdxList);
	[~,idx] = max(numPixels);
	filtered_vol_mask = false(size(filtered_vol_erode));
	filtered_vol_mask(CC.PixelIdxList{idx}) = true;
	SE=strel('cube',dilateSize);
	filtered_vol_erode_dilate = imdilate(filtered_vol_mask,SE);
	figure; orthosliceViewer(img_test .* filtered_vol_erode_dilate); colormap jet; pause(0.5)

	figure; volshow(filtered_vol_erode_dilate); pause(0.5)
	print(gcf,['Seg_BW_after_resize4_', sample_idx, '.png'],'-dpng', '-r300');


% apply mask to soc mapping
	close all
	% img_test = fit_gauss_center_3d_resize4 .* filtered_vol;
	soc_masked_R1_025 = fit_gauss_center_3d_resize4 .* filtered_vol_erode_dilate;
	tomo_masked_R1_025 = fit_gauss_intensity_3d_resize4 .* filtered_vol_erode_dilate;
