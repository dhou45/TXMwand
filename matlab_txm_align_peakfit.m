close all; clc; clear;

fileList = dir('x_eng*.mat'); fileList.name
load(fileList.name);
fileList = dir('crop_reg*resize2*.mat'); fileList.name
load(fileList.name);

sample_idx = '2021AA4';
sizeR = size(crop_reg_xanes_tomo_resize2)
% delete(myCluster.Jobs);
% myCluster = parcluster('local')

%STEP 1, align based on whole image
	fixedVolume = crop_reg_xanes_tomo_resize2(:,:,:,sizeR(4));
	movingRegisteredVolume = zeros([sizeR(1:3),(sizeR(4)-1)]);
	[optimizer,metric] = imregconfig('multimodal');
	Rfixed  = imref3d(size(fixedVolume));
	%OnePlusOne, finer Reg
	time_before = datetime('now')
	optimizer = registration.optimizer.OnePlusOneEvolutionary;
	metric = registration.metric.MattesMutualInformation;
	optimizer.GrowthFactor = 1.01;
	optimizer.Epsilon = 1.0e-6;
	optimizer.InitialRadius = 3.0e-3;
	optimizer.MaximumIterations = 500;
	optimizer
	parfor ii = 1:sizeR(4)-1
		disp(ii)
		% helperVolumeRegistration(fixedVolume,movingVolume);
		Rmoving_temp = imref3d(size(crop_reg_xanes_tomo_resize2(:,:,:,ii)));
		% optimizer.InitialRadius = 0.004;
		movingRegisteredVolume(:,:,:,ii) = imregister(crop_reg_xanes_tomo_resize2(:,:,:,ii),Rmoving_temp, fixedVolume,Rfixed, 'rigid', optimizer, metric);
	end
	time_after = datetime('now')
	disp(time_after-time_before)


%STEP 2, save matlab aligned xanes tomo
	%check the alignment after matlab reg
	figure; sliceViewer(squeeze(movingRegisteredVolume(:,:,round(sizeR(3)/2),:)));
	% figure; sliceViewer(squeeze(movingRegisteredVolume(:,round(sizeR(2)/2),:,:))); %check image from another axies

	matlab_aligned_crop_reg_xanes_tomo_resize2 = zeros(sizeR);
	for i = 1:(sizeR(4)-1)
		matlab_aligned_crop_reg_xanes_tomo_resize2(:,:,:,i) = movingRegisteredVolume(:,:,:,i);
	end
	matlab_aligned_crop_reg_xanes_tomo_resize2(:,:,:,sizeR(4)) = fixedVolume;
	save(['matlab_aligned_crop_reg_xanes_tomo_resize2_',sample_idx,'.mat'], 'matlab_aligned_crop_reg_xanes_tomo_resize2','eng_whiteline','-v7.3')
	close all


%STEP 3, whiteline fitting after matlab alignment
	% choose voxel for xanes profile
	sizeR = size(matlab_aligned_crop_reg_xanes_tomo_resize2)
	figure(1); imshow(squeeze(matlab_aligned_crop_reg_xanes_tomo_resize2(:,:,round(sizeR(3)/2),1)), [])

	Eng_start_index = 7; Eng_end_index = 15;
	jj = 90; % x in matlab plot
	ii= 75; % y in matlab plot
	gaussian_factor = 5; 

	% check pixel xanes spectra
	kk = round(sizeR(3)/2);
	Spectra1 = squeeze(matlab_aligned_crop_reg_xanes_tomo_resize2(ii,jj,kk,:));
	Energy = eng_whiteline/1000; %Energy must in keV!!!
	% figure; scatter(Energy, Spectra1); hold on; plot(Energy, Spectra1); hold off
	Spectra1_smooth = smoothdata(Spectra1,'gaussian',gaussian_factor);
	figure(2); scatter(Energy, Spectra1, 'k'); hold on; scatter(Energy(Eng_start_index: Eng_end_index), Spectra1_smooth(Eng_start_index: Eng_end_index), 'r'); plot(Energy, Spectra1_smooth, 'r'); title("gaussian smooth"); hold off; pause(0.5);

	%check fitted spectra
	Spectra1_crop = double(Spectra1_smooth(Eng_start_index:Eng_end_index));
	Eng_crop = double(Energy(Eng_start_index:Eng_end_index));
	k = 1:length(Eng_crop);
	ki = linspace(1,length(Spectra1_crop),100);
	Eng_interp = interp1(k,Eng_crop,ki,'linear');
	Spectra1_interp = interp1(k,Spectra1_crop,ki,'linear');
	% figure; scatter(Eng_interp, Spectra1_interp); title('interpolation to 100 points')
	xData = Eng_interp.';
	yData = Spectra1_interp.';
	ft = fittype( {'(sin(x-pi))', '((x-10)^2)', '1'}, 'independent', 'x', 'dependent', 'y', 'coefficients', {'a', 'b', 'c'} );
	% Fit linear combination model to data.
	[fitresult, gof] = fit( xData, yData, ft );
	figure(3); plot( fitresult, xData, yData ); title("combined linear fit"); pause(0.5);
	a = fitresult.a; b = fitresult.b; c = fitresult.c; 
	yData_fit = a*((sin(xData-pi))) + b*(((xData-10).^2)) + c;
	[peak_max,peak_index] = max(yData_fit);
	peak_position = xData(peak_index)

	
%STEP 4, batch whiteline fitting
saveas(figure(2),['spectra_matlab_aligned_crop_reg_xanes_resize2_',sample_idx, '.png']);
saveas(figure(3),['spectra_fit_matlab_aligned_crop_reg_xanes_resize2',sample_idx, '.png']); close all;

time_before = datetime('now') 
fit_gauss_center_3d_resize2 = zeros([sizeR(1), sizeR(2), sizeR(3)]);
fit_gauss_adjrsquare_3d_resize2 = zeros([sizeR(1), sizeR(2), sizeR(3)]);
fit_gauss_intensity_3d_resize2 = zeros([sizeR(1), sizeR(2), sizeR(3)]);
for ii = 1 : sizeR(1)
    if rem(ii, 10) == 1    % (ii/5) == round(ii/5)
    	time_after = datetime('now'); 
		disp([num2str(ii/sizeR(1)*100),'%  uptime:  ' ,datestr(time_after-time_before, 'HH:MM:SS'),'  total est:  ' ,datestr((time_after-time_before)/(ii/sizeR(1)), 'HH:MM:SS')])
	end
	for jj = 1 : sizeR(2)
		parfor kk = 1 : sizeR(3)
			Spectra1 = squeeze(matlab_aligned_crop_reg_xanes_tomo_resize2(ii,jj,kk,:));
			Spectra1_smooth = smoothdata(Spectra1,'gaussian',gaussian_factor);
			% figure; scatter(Energy, Spectra1); hold on; plot(Energy, Spectra1); title("gaussian smooth")
			Spectra1_crop = double(Spectra1_smooth(Eng_start_index:Eng_end_index));
			Eng_crop = double(Energy(Eng_start_index:Eng_end_index));
			k = 1:length(Eng_crop);
			ki = linspace(1,length(Spectra1_crop),100);
			Eng_interp = interp1(k,Eng_crop,ki,'linear');
			Spectra1_interp = interp1(k,Spectra1_crop,ki,'linear');
			% figure; scatter(Eng_interp, Spectra1_interp); title('interpolation to 100 points')
			xData = Eng_interp.';
			yData = Spectra1_interp.';
			ft = fittype( {'(sin(x-pi))', '((x-10)^2)', '1'}, 'independent', 'x', 'dependent', 'y', 'coefficients', {'a', 'b', 'c'} );
			% Fit linear combination model to data.
			[fitresult, gof] = fit( xData, yData, ft );
			a = fitresult.a; b = fitresult.b; c = fitresult.c; 
			yData_fit = a*((sin(xData-pi))) + b*(((xData-10).^2)) + c;
			[peak_max,peak_index] = max(yData_fit);
			peak_position = xData(peak_index);
			% Fit gauss model to data.
			% ft = fittype( 'gauss1' );
			% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
			% opts.Display = 'Off';
			% % Fit model to data.
			% [fitresult, gof] = fit( xData, yData, ft, opts );
			% % Plot fit with data.
			% figure( 'Name', 'gaussian fit 1' );
			% h = plot( fitresult, xData, yData );
			fit_gauss_center_3d_resize2(ii,jj,kk) = peak_position;
			fit_gauss_intensity_3d_resize2(ii,jj,kk) = peak_max;
			fit_gauss_adjrsquare_3d_resize2(ii,jj,kk) = gof.adjrsquare;
		end
	end
end

save(strcat('fit_gauss_3d_resize2_matlab_aligned_',sample_idx,'.mat'),'fit_gauss_center_3d_resize2', 'fit_gauss_adjrsquare_3d_resize2','fit_gauss_intensity_3d_resize2');
time_after = datetime('now')
disp(time_after-time_before)


% resize4 and fit if you want


