% read raw h5 energy list
	work_path = uigetdir(); 
	[work_path,filename,ext] = fileparts(filename);

	work_path = '/run/media/VTLinlab/Lin_Lab_10/MR_3DTXM/Mn_Pristine/';
	cd(work_path);
	files=dir(['fly_scan_id*.h5']);

	[work_path,filename,ext] = fileparts([files(1).folder, '/',files(1).name]);
	ind_s_file = filename;
	ind_s = str2num(ind_s_file(end-4:end))

	numberofE = length(files); % number of energy points
	x_eng = zeros(numberofE, 1);
	proj_num = zeros(numberofE, 1);
	for ii=1:numberofE
		disp(ii/numberofE*100)
	    fname=files(ii).name;
	    x_eng(ii) = h5read([fname],'/X_eng');
	    temp = size(h5read([fname],'/img_tomo'));
	    proj_num(ii) = temp(3);
	end

	ind_e = ind_s+length(x_eng)-1;
	ind = ind_s:ind_e;
	info_tomo = [ind', x_eng, proj_num];
	csvwrite('info_tomo.csv', info_tomo);
	figure; scatter(ind_s:ind_e, x_eng)
	% save(strcat('x_eng','.mat'),'x_eng');
	print(gcf, 'eng_list','-dpng','-r600');

	large_combined_h5_file = dir([work_path,'/',filename,ext]);
	fname = large_combined_h5_file.name
	angles = h5read([work_path,'/', fname],'/angle');
	figure; scatter(1:length(angles), angles)
	csvwrite('angles_ind_1.csv', angles);
	print(gcf, 'angles_ind_1','-dpng','-r600');


% plot angle vs proj#
	[filename, work_path] = uigetfile('*.h5');
	cd(work_path);
	large_combined_h5_file = dir([work_path,filename]);
	fname = large_combined_h5_file.name;
	angles = h5read([work_path,fname],'/angle');
	figure; scatter(1:length(angles), angles)
	csvwrite('angles_ind_1.csv', angles);


clear; clc

% read, resize and save mat for aligned zoomed regions in h5 file from TXMgui, only for whiteline recons
	[filename, work_path] = uigetfile('*.h5');
	% fullfilename = fullfile(work_path, filename)
	cd(work_path);
	large_combined_h5_file = dir([work_path,filename]);
	fname = large_combined_h5_file.name;
	x_eng = h5read([work_path,fname],'/registration_results/reg_results/eng_list');
	figure; scatter(1:length(x_eng), x_eng)

	crop_reg_xanes_tomo= single(h5read([work_path,fname],'/registration_results/reg_results/registered_xanes3D'));
	%check alignment after TXMgui
	sizeR = size(crop_reg_xanes_tomo)
	figure; sliceViewer(squeeze(crop_reg_xanes_tomo(:,:,round(sizeR(3)/2),:)));

	
	% resize by factor 2
	resize_factor = 2;
	crop_reg_xanes_tomo_resize2 = zeros([size(imresize3(crop_reg_xanes_tomo(:,:,:,1), 1/resize_factor)), sizeR(4)]);
	for zz = 1 : sizeR(4)
		crop_reg_xanes_tomo_resize2(:,:,:,zz) = imresize3(crop_reg_xanes_tomo(:,:,:,zz), 1/resize_factor);
	end
	% larger than 2GB, save as mat v73
	save(strcat('crop_reg_xanes_tomo_resize2_202104','.mat'),'crop_reg_xanes_tomo_resize2','x_eng','-v7.3');
	clear crop_reg_xanes_tomo;

% resize by factor 4
	resize_factor = 4;
	crop_reg_xanes_tomo_resize4 = zeros([size(imresize3(crop_reg_xanes_tomo(:,:,:,1), 1/resize_factor)), sizeR(4)]);
	for zz = 1 : sizeR(4)
		crop_reg_xanes_tomo_resize4(:,:,:,zz) = imresize3(crop_reg_xanes_tomo(:,:,:,zz), 1/resize_factor);
	end
	% larger than 2GB, save as mat v73
	save(strcat('crop_reg_xanes_tomo_resize4_202104','.mat'),'crop_reg_xanes_tomo_resize4','-v7.3');
