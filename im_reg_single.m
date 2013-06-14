function im_reg_single(filename, line_range, islowbaseline)
%
% Image registration for frames in a single file
% Can be compiled and run on the cluster
% filename, the image data file, can contain full path.
%
% NX -2012

if isdeployed
    line_range = eval(line_range);
    islowbaseline = eval(islowbaseline);
end

if ~exist(filename, 'file')
    error('Data file not exist!')
end

[pthstr,temp_name,ext] = fileparts(filename);
if isempty(pthstr)
    out_path = 'reg_single/';
else
    out_path = [pthstr '/reg_single/'];
end
if ~exist(out_path, 'dir')
    mkdir(out_path);
end

I = strfind(temp_name, 'tif');
out_name = [temp_name(1:end-3) 'greenberg_' temp_name(end-2:end)];
dest = [out_path filesep out_name '.tif'];

info = imfinfo(filename);
if isfield(info(1), 'ImageDescription')
    im_descr = info(1).ImageDescription;
else
    im_descr = '';
end

% Read raw data to im_s_0
% im_s_0 = imread_multi(filename,'g');
[im_s_0, header] = load_scim_data(filename);
% dft registration within the same file.
if islowbaseline == 1
    im_tg_dft = get_active_frame_mean(im_s_0); 
else
    im_tg_dft = mean(im_s_0(:,:,end-9:end),3);
end
im_s_dft = dft_reg(im_s_0, im_tg_dft, line_range);
% disp(['Using subregion of [' num2str(line_range) '] for dft registration!']);
if islowbaseline == 1
    im_t = [];
    disp('Low Baseline! Using preceding frames as TARGET frame! ');
else
    % Use mean of the last 10 frames as the target image for Greenberg Reg.
    im_t = mean(im_s_dft(:,:,end-9:end),3);
    disp('Using the mean of last 10 frame as TARGET frame! ');
end
% im_t = mean(im_s_dft(:,:,:),3);

% Options settings for greenberg registration
opt.norm_int = 0;
opt.settings.move_thresh = 0.06; % .075 % .06 paper
opt.settings.max_iter = 50 ;% 50; % 120 paper too
	
opt.settings.scanlinesperparameter = 2; % 1 (2 paper) % NX 
	%settings.scanlinesperparameter = 4; % 1 (2 paper)
opt.settings.corr_thresh = 0.55; % 0.9; % .75 % .85 paper ... .9 works 4 me
    
opt.settings.pregauss = .75 ; %0.75
opt.settings.haltcorr = 0.99;% .95 ; .99 paper
opt.settings.dampcorr = 0.8; % .8 paper as well

disp(['...............Start ''Greenberg_Reg'' for file ' temp_name ' .....................................']);
warning('off')
[im_c dx_r dy_r E] = imreg_greenberg_nxmod(im_s_dft, im_t, opt);
disp(['===============Done ''Greenberg_Reg'' for file ' temp_name ' ===========================']);

imwrite(uint16(im_c(:,:,1)), dest, 'tif', 'Compression', 'none', 'Description',im_descr, 'WriteMode', 'overwrite');
for f=2:size(im_c,3)
    imwrite(uint16(im_c(:,:,f)), dest, 'tif', 'Compression', 'none', 'WriteMode', 'append');
end

function im_t = get_active_frame_mean(im_s)
    % z-score
    frame_mean_fluo = squeeze(mean(mean(im_s,1), 2));
    zsc = (frame_mean_fluo - mean(frame_mean_fluo))/std(frame_mean_fluo);
    % Find the peak in z-score, and average the 10 frames around the peak
    % as the target frame.
    [~, inds_max] = max(zsc);
    if inds_max<=5
        ind1 = 1;
    else
        ind1 = inds_max-5;
    end
    if inds_max+5 > length(zsc)
        ind2 = length(zsc);
    else
        ind2 = inds_max+5;
    end
    im_t = mean(im_s(:,:, ind1 : ind2),3);
