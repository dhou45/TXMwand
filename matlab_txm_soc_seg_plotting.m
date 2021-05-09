%SoC map plotting, no segmentation
	close all; clc; clear;
	fileList = dir('fit_gauss_intensity*resize4*.mat');fileList.name
	load(fileList.name);
	fileList = dir('fit_gauss_center*resize4*.mat');fileList.name
	load(fileList.name);

sample_idx = 'R2_205_202104';

% plot mask for 95% rsquare
	img_test = fit_gauss_adjrsquare_3d_resize4;
	sizeR = size(img_test);
	noise_level = 0.95;
	% figure; imshow(img_test); colormap gray;
	img_test_bw = img_test;
	img_test_bw(img_test_bw>noise_level)=1;img_test_bw(img_test_bw<=noise_level)=0;
	img_test_bw(~isfinite(img_test_bw))=0;
	fig = figure(9); orthosliceViewer(img_test_bw);

	print(gcf,['noise_map_resize4_', sample_idx, '_',num2str(noise_level), '.png'],'-dpng', '-r300');pause(0.5)

	img_test_bw(~isfinite(img_test_bw))=0;
	noise_map_mask = img_test_bw;

% check soc range
	img_test = fit_gauss_center_3d_resize4;
	% Ni
	eng_min_data = min(img_test,[],'all')	
	eng_max_data = max(img_test,[],'all')
	fig = figure(1); histogram(img_test, 200, 'BinLimits', [eng_min_data+0.0001, eng_max_data-0.0001], 'Normalization', 'probability');

	% Ni	
	eng_min = 8.3460; 
	eng_max = 8.3520; 

	% % Mn
	% eng_min = 6.5540; % 7th
	% eng_max = 6.5630; % 15th

	img_test(img_test < max(eng_min, eng_min_data)+0.0001)= 0; img_test(img_test > min(eng_max, eng_max_data)-0.0001)= 0;
	fig = figure(2); fig = orthosliceViewer(img_test); colormap jet; caxis([eng_min, eng_max]); set(fig,'CrosshairEnable','off'); colorbar('Position', [0.7853 0.1656 0.0258 0.2436]); 
	fig = figure(3); histogram(img_test, 200, 'BinLimits', [max(eng_min, eng_min_data)+0.0001, min(eng_max, eng_max_data)-0.0001], 'Normalization', 'probability');
	title(strcat('mean:  ', num2str(mean(img_test(~(img_test == 0)),'all')),'    median:  ',num2str(median(img_test(~(img_test == 0)),'all'))))

	print('-f2',['whiteline_pos_resize4_', sample_idx, '.png'],'-dpng', '-r300');pause(0.5)
	print('-f3',['whiteline_posdist_resize4_', sample_idx, '.png'],'-dpng', '-r300');pause(0.5)


	%whiteline height in gray scale
	fig = figure(4); orthosliceViewer(fit_gauss_intensity_3d_resize4); colormap gray; pause(0.5)
	print('-f4',['whiteline_height_resize4_', sample_idx, '.png'],'-dpng', '-r300');
	% fig = figure(5); orthosliceViewer(fit_gauss_intensity_3d_resize4); colormap jet; pause(0.5)
	% print('-f5','whiteline_height_jet_resize4_', sample_idx, '.png','-dpng', '-r300');


% whiteline with noise mask
	% img_test = fit_gauss_center_3d_resize4 .* noise_map_mask;
	% img_test(img_test < eng_min+0.0001)= 0; img_test(img_test > eng_max-0.0001)= 0;

	% eng_median = median(img_test(img_test ~=0), 'all');	
	% fig = figure(5); fig = orthosliceViewer(img_test); colormap jet; caxis([eng_min, eng_max]); set(fig,'CrosshairEnable','off'); colorbar('Position', [0.7853 0.1656 0.0258 0.2436]); title(strcat('median: ', num2str(eng_median))); pause(0.5)
	% print('-f5',['whiteline_pos_noise_masked_resize4_', sample_idx, '.png'],'-dpng', '-r300');

	% fig = figure(3); histogram(img_test, 200, 'BinLimits', [eng_min+0.0001, eng_max-0.0001], 'Normalization', 'probability');
	% title(strcat('mean:  ', num2str(mean(img_test(~(img_test == 0)),'all')),'    median:  ',num2str(median(img_test(~(img_test == 0)),'all'))))
	% print('-f3',['whiteline_posdist_noise_masked_resize4_', sample_idx, '.png'],'-dpng', '-r300');pause(0.5)


% SoC map segmentation and plotting
% need 2nd round seg (erode then dilate)
	
		% close all; 	clc; clear;
		% % load fitted results
		% fileList = dir('fit_gauss_intensity*resize4*.mat');fileList.name
		% load(fileList.name);
		% fileList = dir('fit_gauss_center*resize4*.mat');fileList.name
		% load(fileList.name);

	close all; 
	% segmentation, use peak hight as mask
	sizeR = size(fit_gauss_intensity_3d_resize4);
	img_test = fit_gauss_intensity_3d_resize4;
	% img_test(img_delta<0.25*int_max) = 0;
	figure; orthosliceViewer(img_test); colormap jet; pause(0.5)
	% print(gcf,'Seg_before_Gravel1_025_resize4.png','-dpng', '-r300');

	normalizedImage = uint8(255*mat2gray(img_test));
	Binary_Mn = imbinarize(normalizedImage);
	figure; volshow(Binary_Mn); pause(0.5)
	% print(gcf,'Seg_BW_before_Gravel1_025_resize4.png','-dpng', '-r300');

	% figure; orthosliceViewer(Binary_Mn);
	% figure; sliceViewer(Binary_Mn);
	SE=strel('cube',1);
	Binary_Mn_close=imclose(Binary_Mn,SE);
	Binary_Mn_fill = imfill(Binary_Mn_close,'holes');
	% figure; volshow(Binary_Mn_fill);
	% choose the largest particle
	CC = bwconncomp(Binary_Mn_fill, 6);
	numPixels = cellfun(@numel,CC.PixelIdxList);
	[~,idx] = max(numPixels);
	filtered_vol = false(size(Binary_Mn_fill));
	filtered_vol(CC.PixelIdxList{idx}) = true;
	figure; volshow(filtered_vol); pause(0.5)
	figure; orthosliceViewer(img_test .* filtered_vol); colormap jet; pause(0.5)
	% print(gcf,'Seg_after_Gravel1_025_resize4.png','-dpng', '-r300');
	% figure; sliceViewer(filtered_vol);


	% erode x pixel, then choose largest, then dilate back (change erodeSize to smooth the mask, but no inner gaps)
	erodeSize = 5;
	dilateSize = erodeSize; % consider if dilate or not.

	SE=strel('cube',erodeSize);
	filtered_vol_erode = imerode(filtered_vol,SE);
	% figure; volshow(filtered_vol_erode); pause(0.5)
	CC = bwconncomp(filtered_vol_erode, 6);
	numPixels = cellfun(@numel,CC.PixelIdxList);
	[~,idx] = max(numPixels);
	filtered_vol_mask = false(size(filtered_vol_erode));
	filtered_vol_mask(CC.PixelIdxList{idx}) = true;
	% figure; volshow(filtered_vol_mask);
	SE=strel('cube',dilateSize);
	filtered_vol_erode_dilate = imdilate(filtered_vol_mask,SE);
	figure; orthosliceViewer(img_test .* filtered_vol_erode_dilate); colormap jet; pause(0.5)
	% print(gcf,'Seg_after_masked_resize4_', sample_idx, '.png','-dpng', '-r300');

	figure; volshow(filtered_vol_erode_dilate); pause(0.5)
	% print(gcf,['Seg_BW_after_resize4_', sample_idx, '.png'],'-dpng', '-r300');


% apply mask to soc mapping
	close all
	% img_test = fit_gauss_center_3d_resize4 .* filtered_vol;
	soc_masked = fit_gauss_center_3d_resize4 .* filtered_vol_erode_dilate;
	tomo_masked = fit_gauss_intensity_3d_resize4 .* filtered_vol_erode_dilate;
	save(['fit_gauss_masked_resize4_',sample_idx, '.mat'], 'soc_masked', 'tomo_masked')

	%whiteline height, masked in gray scale
	% fig = figure(4); orthosliceViewer(tomo_masked); colormap gray; pause(0.5)
	% print('-f4',['whiteline_height_masked_resize4_', sample_idx, '.png'],'-dpng', '-r300');

	%whiteline height, no mask in gray scale
	fig = figure(4); orthosliceViewer(fit_gauss_intensity_3d_resize4); colormap gray; pause(0.5)
	print('-f4',['whiteline_height_resize4_', sample_idx, '.png'],'-dpng', '-r300');
	% fig = figure(5); orthosliceViewer(fit_gauss_intensity_3d_resize4); colormap jet; pause(0.5)
	% print('-f5','whiteline_height_jet_resize4_', sample_idx, '.png','-dpng', '-r300');
	

	img_test = soc_masked;
	% Ni
	eng_max_data = max(img_test,[],'all')
	eng_min_data = min(img_test,[],'all')

	eng_min = 8.3460; 
	eng_max = 8.3520; 
	% fig = figure(1); histogram(img_test, 200, 'BinLimits', [eng_min_data+0.0001, eng_max_data-0.0001], 'Normalization', 'probability'); 

	img_test(img_test < max(eng_min, eng_min_data)+0.0001)= 0; img_test(img_test > min(eng_max, eng_max_data)-0.0001)= 0;
	fig = figure(2); fig = orthosliceViewer(img_test); colormap jet; caxis([eng_min, eng_max]); set(fig,'CrosshairEnable','off'); colorbar('Position', [0.7853 0.1656 0.0258 0.2436]); 
	fig = figure(3); histogram(img_test, 200, 'BinLimits', [max(eng_min, eng_min_data)+0.0001, min(eng_max, eng_max_data)-0.0001], 'Normalization', 'probability');
	title(strcat('mean:  ', num2str(mean(img_test(~(img_test == 0)),'all')),'    median:  ',num2str(median(img_test(~(img_test == 0)),'all'))))
	print('-f2',['whiteline_pos_masked_resize4_', sample_idx, '.png'],'-dpng', '-r300');pause(0.5)
	print('-f3',['whiteline_posdist_masked_resize4_', sample_idx, '.png'],'-dpng', '-r300');pause(0.5)

	% 3ev, smaller eng range, median around 8.3515
	eng_median = median(img_test(img_test ~=0), 'all');	
	eng_min = eng_median-0.0015;
	eng_max = eng_median+0.0015;
	fig = figure(5); fig = orthosliceViewer(img_test); colormap jet; caxis([eng_min, eng_max]); set(fig,'CrosshairEnable','off'); colorbar('Position', [0.7853 0.1656 0.0258 0.2436]); title(strcat('median: ', num2str(eng_median))); pause(0.5)
	print('-f5',['whiteline_pos_masked_resize4_', sample_idx, '_gap3eV.png'],'-dpng', '-r300');


% hist fit
	img_test = soc_masked;
	img_mu = mean(img_test(img_test ~=0));
	img_sigma = std(img_test(img_test ~=0));
	% mode(img_test(img_test ~=0))
	% var(img_test(img_test ~=0))
	eng_min = 8.3495; 
	eng_max = 8.353; 
	eng_max_data = max(img_test,[],'all')
	eng_min_data = min(img_test,[],'all')
	img_test(img_test < max(eng_min, eng_min_data)+0.0001)= 0; img_test(img_test > min(eng_max, eng_max_data)-0.0001)= 0;
	[counts1,edges1] = histcounts(img_test, 200, 'BinLimits', [eng_min, eng_max], 'Normalization', 'probability');
	positions1 = edges1(1:end-1) + diff(edges1) / 2;
	counts = [positions1(counts1>0);counts1(counts1>0)]; counts = counts'; 
	figure; scatter(counts(:,1), counts(:,2), 'go');  
	csvwrite(['soc_hist_',sample_idx,'.csv'],counts);

	% one gauss profile fit
		f = fit(counts(:,1),counts(:,2),'gauss1')
		hold on; plot(f, 'r-'); title([num2str(f.b1)]); hold off
		print(gcf,['hist_fit_gauss1_', sample_idx, '.png'],'-dpng', '-r300');


	% two gauss fit
		figure; scatter(counts(:,1), counts(:,2), 'go');  
		options = fitoptions('gauss2');
		options.Lower = [0 8.35050 0 0 -Inf 0];
		options.Upper = [inf 8.35175 inf inf Inf inf];
		f = fit(counts(:,1),counts(:,2),'gauss2', options)
		hold on; plot(f, 'r-'); title([num2str(f.b1), '    ', num2str(f.b2)]); hold off
		x = counts(:,1);
		peak2 = f.a2*exp(-((x-f.b2)/f.c2).^2);
		peak1 = f.a1*exp(-((x-f.b1)/f.c1).^2);
		hold on; plot(x,peak1, x, peak2)
		print(gcf,['hist_fit_gauss2_', sample_idx, '.png'],'-dpng', '-r300');




% +/- 2ev, smaller eng range
	eng_median = median(img_test(img_test ~=0), 'all');
	eng_min = eng_median-0.002; 
	eng_max = eng_median+0.002;
	fig = figure(6); fig = orthosliceViewer(img_test); colormap jet; caxis([eng_min, eng_max]); set(fig,'CrosshairEnable','off'); colorbar('Position', [0.7853 0.1656 0.0258 0.2436]); title(strcat('median: ', num2str(eng_median))); pause(0.5)
	print('-f6',['whiteline_pos_masked_resize4_', sample_idx, '_gap4eV.png'],'-dpng', '-r300');
