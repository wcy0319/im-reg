function dft_reg_dir(src_dir, save_path, main_fname, targetImage, varargin)
% shift = batch_dft_reg(source_filenames, targetImage, padding_flag,save_path)
%
% Perform dft whole frame registration for the a whole directory. All image
%   files will be registered to the same target image.
%
% Optional input: main_fname, register imaging files only containing this
%                 main_fname. Leave [] would take all files in the dir.
%                 targetImage, If leave [], use default target image. 
%
% varargin{1}, 1x2 vector specify the range of rows in the image used for
%               wlhole frame registration.
% varargin{2}, save_to_image, default, 1. If 1, save results to tiff files.
%                   
% - NX 2012-04-10
% 
% Update, use ScanImage 4 data, with Tiff class, and int16 data format.
% - NX 2013-06-12
%

if isempty(varargin) || isempty(varargin{1})
    line_range = [];
else
    line_range = varargin{1};
end

if length(varargin) > 1
    save_to_image = varargin{2};
else
    save_to_image = 1;
end 

if  nargin < 2 || isempty(save_path)
    save_path = [src_dir filesep 'dft_reg'];
end

if ~isdir(save_path)
    mkdir(save_path);
end

if nargin < 3 || isempty(main_fname)
    main_fname = '*';
end

if nargin < 4 || isempty(targetImage)
    default_target = 1;
else
    default_target = 0;
end


fprintf('Output path: %s\n', save_path);

datafiles = dir([src_dir filesep main_fname '*.tif']);
fprintf('total data files %d\n', length(datafiles));

% DETERMINE THE TARGET IMAGE
if default_target == 1
    %%% mean image of trial No 10 (arbitrary) as target image
    if length(datafiles) > 10
        n = 10;
    else
        n = length(datafiles);
    end
    [im_temp, ~] = load_scim_data([src_dir filesep datafiles(n).name]);
    im_temp_reg = dft_reg(im_temp, [], line_range);
    im_tg = mean(im_temp_reg, 3);
else
    im_tg = targetImage;
end

if matlabpool('size') == 0
    matlabpool open 8
end

parfor i = 1:length(datafiles)
    %%%% Determine source file names and output file names and directory
    %%%% and data file header info.
    a =  datafiles(i).name;
    % Check whether the file name fits the trial data file
    % "main_fname_000.tif"
    if ~strcmp(main_fname, '*')
        if length(a) == length(main_fname) + 7
            fname = [src_dir filesep a];
        end
    else
        fname = [src_dir filesep a];
    end
    file_basename = a(1: end-7);
    
    im_info = imfinfo(fname);
    if isfield(im_info(1),'ImageDescription')
        im_describ = im_info(1).ImageDescription; % to be put back to the header
    else
        im_describ = '';
    end
    
    % Change channel number if necessary, the out put should be single channel
    % data
    if ~isempty(strfind(im_describ,'numberOfChannelsSave=2'))
        im_describ = strrep(im_describ, 'numberOfChannelsSave=2','numberOfChannelsSave=1');
    end
    if ~isempty(strfind(im_describ, 'saveDuringAcquisition=1'))&& ~isempty(strfind(im_describ, 'numberOfChannelsAcquire=2'))
        im_describ = strrep(im_describ, 'numberOfChannelsAcquire=2','numberOfChannelsAcquire=1');
    end
    if ~isempty(strfind(im_describ, 'channelsSave = 2'))
        im_describ = strrep(im_describ, 'channelsSave = 1');
    end
    
    save_name = [save_path filesep  file_basename 'dftReg_' fname(end-6:end-4) '.tif'];
    
    % %%% READ SOURCE IMAGE DATA
    [im_s, header] = load_scim_data(fname);
    imTagStruct = get_tiff_tag_to_struct(fname);
    
    % DFT REGISTRATION
    im_dft_reg = dft_reg(im_s, im_tg);
    
    if save_to_image == 1
        write_data_to_tiff(save_name, im_dft_reg, imTagStruct);
%         imwrite(uint16(im_dft_reg(:,:,1)), save_name, 'tif', 'Compression', 'none', 'Description',im_describ, 'WriteMode', 'overwrite');
%         for f=2:size(im_dft_reg,3)
%             imwrite(uint16(im_dft_reg(:,:,f)), save_name, 'tif', 'Compression', 'none', 'WriteMode', 'append');
%         end
        fprintf('Results saved to %s\n', save_name);
    end
%     dft_reg_trial = struct([]);
%     dft_reg_trial(1).filename = [file_basename 'dftReg_' fname(end-6:end-4)];
%     dft_reg_trial(1).dft_im_data = im_dft_reg;
%     parsave_dft_dir(dft_reg_trial(1).filename, dft_reg_trial);
end


% save([save_path 'im_dft_reg_results'], 'im_dft_reg_results','-v7.3');

