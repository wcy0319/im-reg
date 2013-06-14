function batch_greenberg_reg2(src_path, main_fname, out_path, varargin)

% Run image registration on multiple trials in the same folder, specified
% by the main_fname. Save registered image files to "out_path".
% Should first start "matlabpool open 8" to open up multiple cores for
% parallel processing.

% wholeFrame_reg_flag = varargin{1}; An option that skips whole frame dft registration
% 

if ~isempty(varargin)
    wholeFrame_reg_flag = varargin{1};
else
    wholeFrame_reg_flag = 1;
end

if wholeFrame_reg_flag == 1
    datafiles = dir([src_path filesep main_fname '*.tif']);
    fprintf('total data files %d\n', length(datafiles));
    if ~exist('out_path','var')
        out_path = [src_path filesep main_fname '_reg'];
    elseif ~isdir(out_path)
        mkdir(out_path);
    end
    
    % out_path = [out_path filesep main_fname '_reg'];
    
    fprintf('Output path: %s\n', out_path);
    fnames = {};
    for i = 1:length(datafiles)
        fname =  datafiles(i).name;
        % Check whether the file name fits the trial data file
        % "main_fname_000.tif"
        if ~strcmp(main_fname, '*')
            if length(fname) == length(main_fname) + 7
                fnames{i} = [src_path filesep fname];
            end
        else
            fnames{i} = [src_path filesep fname];
        end
    end
    
    % Take mean image of trial No 10 (arbitrary) as target image for the whole session
    im_temp = imread_multi(fnames{10},'g');
    
    im_tg = dft_reg(im_temp);
    
    im_tg_mean = mean(im_tg,3);
    
    % dft_shift = {}; im_descr = {};
    dft_dir = [src_path filesep 'dft_reg'];
    
    dft_shift_all = batch_dft_reg(im_tg_mean, fnames, 0, dft_dir);
    dft_files = dir([dft_dir filesep main_fname, '*.tif']);
    
else
    dft_dir = src_path;
    dft_files = dir([dft_dir filesep main_fname, '*.tif']);
    fprintf('Skip whole frame reg, start greenberg reg ....\nFirst file: %s\n', dft_files(1).name);
    if isempty(findstr(dft_files(1).name, 'dftReg'))
        warning('The images seems not registered by dft, still proceed?')
        ans = input('Y / N :');
        if strcmpi(ans, 'N')
            return
        end
    end
    dft_shift_all = [];
end

if matlabpool('size') == 0
    matlabpool open 8
end

parfor i=1:length(dft_files)
    % if isdir(datafiles(i).name)
      %  continue;
    % end
    filename = [dft_dir filesep dft_files(i).name];
    info = imfinfo(filename);
    if isfield(info(1), 'ImageDescription')
        im_descr = info(1).ImageDescription;
    else
        im_descr = '';
    end
    
    im_s = imread_multi(filename,'g');
    
    % Use as target the mean of late frames of the trial, which is more
    % still, and not illuminated by light stimulation.
    im_t = mean(im_s(:,:,89:end),3); %mean(im_s,3);
    
    disp(['Start ''greenberg_reg'' for file ' dft_files(i).name ' ...']);
    [im_c dx_r dy_r E] = imreg_greenberg(im_s, im_t, []);
    
    [pthstr,temp_name,ext] = fileparts(filename);
    out_name = [temp_name(1:end-3) 'greenberg_' temp_name(end-2:end)];
    dest = [out_path filesep out_name];
    
    imwrite(uint16(im_c(:,:,1)), [dest '.tif'], 'tif', 'Compression', 'none', 'Description',im_descr, 'WriteMode', 'overwrite');
    for f=2:size(im_c,3)
        imwrite(uint16(im_c(:,:,f)), [dest '.tif'], 'tif', 'Compression', 'none', 'WriteMode', 'append');
    end
    if ~isempty(dft_shift_all)
        dft_shift = dft_shift_all(:,:,i); % [];%
    else 
        dft_shift = [];
    end
    parsave([dest '_reginfo'], dx_r,dy_r,E,dft_shift);
    fprintf('Reg data saved to %s\n',out_name);
end
if matlabpool('size') > 0
    matlabpool close
end

function parsave(fn, dx_r, dy_r, E, dft_shift)
    save(fn, 'dx_r','dy_r','E','dft_shift');


function [im_dft, shift] = dft_reg(im_s)
im_tg = mean(im_s(:,:,1:10),3);
for i=1:size(im_s,3);
    output(:,i) = dftregistration(fft2(double(im_tg)),fft2(double(im_s(:,:,i))),1);
end
shift = output(3:4,:);
padding = [0 0 0 0];
im_dft = ImageTranslation_nx(im_s,shift,padding,0);
        
