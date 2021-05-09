sample_idx = 'BA';

	
	% segmentation, use peak hight as mask
	sizeR = size(fit_gauss_intensity_3d_resize2);
	img_test = fit_gauss_intensity_3d_resize2;
	% img_test(img_delta<0.25*int_max) = 0;
	figure; orthosliceViewer(img_test); colormap jet; pause(0.5)
	% print(gcf,'Seg_before_Gravel1_025_resize2.png','-dpng', '-r300');

	normalizedImage = uint8(255*mat2gray(img_test));
	Binary_Mn = imbinarize(normalizedImage);
	figure; volshow(Binary_Mn); pause(0.5)
	% print(gcf,'Seg_BW_before_Gravel1_025_resize2.png','-dpng', '-r300');

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
	% print(gcf,'Seg_after_Gravel1_025_resize2.png','-dpng', '-r300');
	% figure; sliceViewer(filtered_vol);


	% erode x pixel, then choose largest, then dilate back (change erodeSize to smooth the mask, but no inner gaps)
	erodeSize = 5;
	dilateSize = erodeSize;

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
	% print(gcf,'Seg_after_masked_resize2_', sample_idx, '.png','-dpng', '-r300');

	figure; volshow(filtered_vol_erode_dilate); pause(0.5)



% apply mask to soc mapping
	close all
	% img_test = fit_gauss_center_3d_resize2 .* filtered_vol;
	soc_masked = fit_gauss_center_3d_resize2 .* filtered_vol_erode_dilate;
	tomo_masked = fit_gauss_intensity_3d_resize2 .* filtered_vol_erode_dilate;

	%whiteline height in gray scale
	fig = figure(4); orthosliceViewer(tomo_masked); colormap gray; pause(0.5)
	print('-f4',['whiteline_height_masked_resize2_', sample_idx, '.png'],'-dpng', '-r300');
	save(['fit_gauss_masked_resize2_',sample_idx, '.mat'], 'soc_masked', 'tomo_masked')

	img_test = soc_masked;
	img_mu = mean(img_test(img_test ~=0));
	img_sigma = std(img_test(img_test ~=0));
	% mode(img_test(img_test ~=0))
	% var(img_test(img_test ~=0))
	
	img_test = soc_masked;
	eng_min = 8.350; 
	eng_max = 8.3535; 
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


% upper 2.2% distribution, 1sigma, full size
	sample_idx = 'BBMFCD70';
	img_pairs = soc_masked;
	img_pairs(img_pairs<(img_mu+2*img_sigma))= 0; img_pairs(img_pairs>(img_mu+3*img_sigma))=0;
	img_pairs(img_pairs ~=0)=1;
	figure; volshow(img_pairs)
	save(['img_pairs_',sample_idx, '.mat'], 'img_pairs')

	% soc pdf plots, full size
	sample_idx = 'BBMFCD70';
	sizeR = size(img_pairs);
	X = []; Y = []; Z = [];
	% poolobj = parpool;
	time_before = datetime('now'); 
	for ii = 1: sizeR(1)
	    if rem(ii, 10) == 1    % (ii/5) == round(ii/5)
	    	time_after = datetime('now'); 
			disp([num2str(ii/sizeR(1)*100),'%  uptime:  ' ,datestr(time_after-time_before, 'HH:MM:SS'),'  total est:  ' ,datestr((time_after-time_before)/(ii/sizeR(1)), 'HH:MM:SS')])
		end
		parfor jj = 1:sizeR(2)
			for kk = 1:sizeR(3)
				if img_pairs(ii, jj, kk) == 1
					X = [X, ii]; Y = [Y, jj]; Z = [Z, kk];
				end
			end
		end
	end
	% delete(poolobj);
	img_clusters = [X; Y; Z]';
	% save('soc_clusters_positions.mat', 'img_clusters');

	distances = pdist2(img_clusters,img_clusters);
	% sum (distances, 'all')/max(distances,[],'all')
	% figure; histogram(distances, 200, 'BinLimits', [0.1, max(distances,[],'all')], 'Normalization', 'probability');
	norm_distances = distances./max(distances,[],'all');
	figure; histogram(norm_distances, 200, 'BinLimits', [0, max(norm_distances,[],'all')], 'Normalization', 'probability');
	[counts1,edges1] = histcounts(norm_distances, 200, 'BinLimits', [0, max(norm_distances,[],'all')], 'Normalization', 'probability');
	positions1 = edges1(1:end-1) + diff(edges1) / 2;
	counts = [positions1(counts1>0);counts1(counts1>0)]; counts = counts'; 
	hold on; scatter(counts(:,1), counts(:,2), 'go'); 
	print(gcf,['pdf_norm_', sample_idx, '.png'],'-dpng', '-r300');
	csvwrite(['cluster_counts_norm_',sample_idx,'.csv'],counts);

	close all


	% plot absolute distances
	figure; histogram(distances, 200, 'BinLimits', [0.1, max(distances,[],'all')], 'Normalization', 'probability');
	[counts1,edges1] = histcounts(distances, 200, 'BinLimits', [0, max(distances,[],'all')], 'Normalization', 'probability');
	positions1 = edges1(1:end-1) + diff(edges1) / 2;
	counts = [positions1(counts1>0);counts1(counts1>0)]; counts = counts'; 
	hold on; scatter(counts(:,1), counts(:,2), 'go'); 
	print(gcf,['pdf_abs_', sample_idx, '.png'],'-dpng', '-r300');
	csvwrite(['cluster_counts_abs_',sample_idx,'.csv'],counts);



% soc pdf plots, resize2 again.
	sample_idx = 'BBMFCD70_resize2';
	img_pairs_resize2_temp = imresize3(img_pairs, 0.5);
	img_pairs_resize2 = imbinarize(img_pairs_resize2_temp);
	sizeR = size(img_pairs_resize2);

	figure; orthosliceViewer(img_pairs);
	figure; orthosliceViewer(img_pairs_resize2);
	% poolobj = parpool;

	X = []; Y = []; Z = [];
	time_before = datetime('now'); 
	for ii = 1: sizeR(1)
	    if rem(ii, 10) == 1    % (ii/5) == round(ii/5)
	    	time_after = datetime('now'); 
			disp([num2str(ii/sizeR(1)*100),'%  uptime:  ' ,datestr(time_after-time_before, 'HH:MM:SS'),'  total est:  ' ,datestr((time_after-time_before)/(ii/sizeR(1)), 'HH:MM:SS')])
		end
		parfor jj = 1:sizeR(2)
			for kk = 1:sizeR(3)
				if img_pairs_resize2(ii, jj, kk) == 1
					X = [X, ii]; Y = [Y, jj]; Z = [Z, kk];
				end
			end
		end
	end
	% delete(poolobj);
	img_clusters = [X; Y; Z]';
	% save(['soc_clusters_positions_',sample_idx,'.mat'], 'img_clusters');

	distances = pdist2(img_clusters,img_clusters);
	% sum (distances, 'all')/max(distances,[],'all')
	norm_distances = distances./max(distances,[],'all');
	figure; histogram(norm_distances, 200, 'BinLimits', [0, max(norm_distances,[],'all')], 'Normalization', 'probability');
	[counts1,edges1] = histcounts(norm_distances, 200, 'BinLimits', [0, max(norm_distances,[],'all')], 'Normalization', 'probability');
	positions1 = edges1(1:end-1) + diff(edges1) / 2;
	counts = [positions1(counts1>0);counts1(counts1>0)]; counts = counts'; 
	hold on; scatter(counts(:,1), counts(:,2), 'go'); 
	print(gcf,['pdf_norm_', sample_idx, '.png'],'-dpng', '-r300');
	csvwrite(['cluster_counts_norm_',sample_idx,'.csv'],counts);



% upper 15.9% distribution, 1sigma, and resize2
	sample_idx = 'BACD70_1sigma';

	img_test = soc_masked;
	img_mu = mean(img_test(img_test ~=0));
	img_sigma = std(img_test(img_test ~=0));
	img_pairs = soc_masked;
	img_pairs(img_pairs<(img_mu+1*img_sigma))= 0; img_pairs(img_pairs>(img_mu+3*img_sigma))=0;
	img_pairs(img_pairs ~=0)=1;
	figure; volshow(img_pairs)
	save(['img_pairs_',sample_idx, '.mat'], 'img_pairs')

	% soc pdf plots, resize2 again.
	sample_idx = 'BACD70_1sigma_resize2';
	img_pairs_resize2_temp = imresize3(img_pairs, 0.5);
	img_pairs_resize2 = imbinarize(img_pairs_resize2_temp);
	sizeR = size(img_pairs_resize2);

	figure; orthosliceViewer(img_pairs);
	figure; orthosliceViewer(img_pairs_resize2);
	% poolobj = parpool;

	X = []; Y = []; Z = [];
	time_before = datetime('now'); 
	for ii = 1: sizeR(1)
	    if rem(ii, 10) == 1    % (ii/5) == round(ii/5)
	    	time_after = datetime('now'); 
			disp([num2str(ii/sizeR(1)*100),'%  uptime:  ' ,datestr(time_after-time_before, 'HH:MM:SS'),'  total est:  ' ,datestr((time_after-time_before)/(ii/sizeR(1)), 'HH:MM:SS')])
		end
		parfor jj = 1:sizeR(2)
			for kk = 1:sizeR(3)
				if img_pairs_resize2(ii, jj, kk) == 1
					X = [X, ii]; Y = [Y, jj]; Z = [Z, kk];
				end
			end
		end
	end
	% delete(poolobj);
	img_clusters = [X; Y; Z]';
	% save(['soc_clusters_positions_',sample_idx,'.mat'], 'img_clusters');
	distances = pdist2(img_clusters,img_clusters);
	% sum (distances, 'all')/max(distances,[],'all')
	norm_distances = distances./max(distances,[],'all');
	figure; histogram(norm_distances, 200, 'BinLimits', [0, max(norm_distances,[],'all')], 'Normalization', 'probability');
	[counts1,edges1] = histcounts(norm_distances, 200, 'BinLimits', [0, max(norm_distances,[],'all')], 'Normalization', 'probability');
	positions1 = edges1(1:end-1) + diff(edges1) / 2;
	counts = [positions1(counts1>0);counts1(counts1>0)]; counts = counts'; 
	hold on; scatter(counts(:,1), counts(:,2), 'go'); 
	print(gcf,['pdf_norm_', sample_idx, '.png'],'-dpng', '-r300');
	csvwrite(['cluster_counts_norm_',sample_idx,'.csv'],counts);

	close all


% lower 2.2% distribution, 2sigma, full size
	sample_idx = 'BAMFCD80';
	img_pairs = soc_masked;
	img_mu = mean(img_pairs(img_pairs ~=0));
	img_sigma = std(img_pairs(img_pairs ~=0));
	img_pairs(img_pairs<(img_mu-3*img_sigma))= 0; img_pairs(img_pairs>(img_mu-2*img_sigma))=0;
	img_pairs(img_pairs ~=0)=1;
	% fig = figure; volshow(img_pairs);
	% print(gcf,['soc_low_1sigma_3D_', sample_idx, '.png'],'-dpng', '-r300'); pause(0.5)
	figure; orthosliceViewer(img_pairs); title(['-2sigma', num2str(img_mu-2*img_sigma)])
	print(gcf,['soc_low_2sigma_2D_', sample_idx, '.png'],'-dpng', '-r300'); pause(0.5)

	save(['soc_low_2sigma_',sample_idx, '.mat'], 'img_pairs')

% upper 2.2% distribution, 2sigma, full size
	img_pairs = soc_masked;
	img_mu = mean(img_pairs(img_pairs ~=0));
	img_sigma = std(img_pairs(img_pairs ~=0));
	img_pairs(img_pairs<(img_mu+2*img_sigma))= 0; img_pairs(img_pairs>(img_mu+3*img_sigma))=0;
	img_pairs(img_pairs ~=0)=1;
	% fig = figure; volshow(img_pairs);
	% print(gcf,['soc_high_1sigma_3D_', sample_idx, '.png'],'-dpng', '-r300'); pause(0.5)
	figure; orthosliceViewer(img_pairs); title(['2sigma: ', num2str(img_mu+2*img_sigma)])
	print(gcf,['soc_high_2sigma_2D_', sample_idx, '.png'],'-dpng', '-r300'); pause(0.5)
	save(['soc_high_2sigma_',sample_idx, '.mat'], 'img_pairs')




% lower 15% distribution, 1sigma, full size
clc; clear all; close all;

	sample_idx = 'BBMFCD70';
	img_pairs = soc_masked;
	img_mu = mean(img_pairs(img_pairs ~=0));
	img_sigma = std(img_pairs(img_pairs ~=0));
	img_pairs(img_pairs<(img_mu-3*img_sigma))= 0; img_pairs(img_pairs>(img_mu-1*img_sigma))=0;
	img_pairs(img_pairs ~=0)=1;
	% fig = figure; volshow(img_pairs);
	% print(gcf,['soc_low_1sigma_3D_', sample_idx, '.png'],'-dpng', '-r300'); pause(0.5)
	figure; orthosliceViewer(img_pairs); title(['-1sigma', num2str(img_mu-1*img_sigma)])
	print(gcf,['soc_low_1sigma_2D_', sample_idx, '.png'],'-dpng', '-r300'); pause(0.5)

	save(['soc_low_1sigma_',sample_idx, '.mat'], 'img_pairs')

% upper 15% distribution, 1sigma, full size
	img_pairs = soc_masked;
	img_mu = mean(img_pairs(img_pairs ~=0));
	img_sigma = std(img_pairs(img_pairs ~=0));
	img_pairs(img_pairs<(img_mu+1*img_sigma))= 0; img_pairs(img_pairs>(img_mu+3*img_sigma))=0;
	img_pairs(img_pairs ~=0)=1;
	% fig = figure; volshow(img_pairs);
	% print(gcf,['soc_high_1sigma_3D_', sample_idx, '.png'],'-dpng', '-r300'); pause(0.5)
	figure; orthosliceViewer(img_pairs); title(['1sigma: ', num2str(img_mu+1*img_sigma)])
	print(gcf,['soc_high_1sigma_2D_', sample_idx, '.png'],'-dpng', '-r300'); pause(0.5)

	save(['soc_high_1sigma_',sample_idx, '.mat'], 'img_pairs')


