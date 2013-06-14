function batch_greenberg_reg3(src_path, main_fname, out_path, varargin)

% Run image registration on multiple trials in the same folder, specified
% by the main_fname. Save registered image files to "out_path".
% Should first start "matlabpool open 8" to open up multiple cores for
% parallel processing.

%
% varargin{1}, 1x2 vector specify the range of rows in the image used for
%               wlhole frame registration.
% varargin{2}, imageSize, default [128 512]. This is used to determine
%               whether to skip image files.
% varargin{3}, minNumFrames, default 2. Need at least 2 frames for registration.


if isempty(varargin) || isempty(varargin{1})
    line_range = [];
else
    line_range = varargin{1};
end

if length(varargin)<2 || isempty(varargin{2})
    imageSize = [128 512];
else
    imageSize = varargin{2};
end
if length(varargin)<3 || isempty(varargin{3})
    minNumFrames = 2;
else
    minNumFrames = varargin{3};
end

if matlabpool('size') == 0
    matlabpool open 8
end

% Get data files, and setup output directory
datafiles = dir([src_path filesep main_fname '*.tif']);
fprintf('total data files %d\n', length(datafiles));
if ~exist('out_path','var')
    out_path = [src_path filesep main_fname '_reg'];
elseif ~isdir(out_path)
    mkdir(out_path);
end
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

parfor i=1:length(datafiles)
    % Get the full file name with path
    a =  datafiles(i).name;
    if ~strcmp(main_fname, '*')
        % Check whether the file name fits the trial data file
        % "main_fname_000.tif"
        if length(a) == length(main_fname) + 7
            filename = [src_path filesep a];
        end
    else
        filename = [src_path filesep a];
    end
    
    info = imfinfo(filename);
    % Check if file is a valid data file
    if info(1).Height == imageSize(1) && length(info) >= minNumFrames
        isvalidfile = 1;
    else
        isvalidfile = 0;
    end
    if isvalidfile == 1
        if isfield(info(1), 'ImageDescription')
            im_descr = info(1).ImageDescription;
        else
            im_descr = '';
        end
        
        % Read raw data to im_s_0
        im_s_0 = imread_multi(filename,'g');
        
        % dft registration within the same file.
        im_s_dft = dft_reg(im_s_0, [], line_range);
        
        % Use mean of the last 10 frames as the target image for Greenberg Reg.
        im_t = mean(im_s_dft(:,:,end-9:end),3);
        
        
        [pthstr,temp_name,ext] = fileparts(filename);
        out_name = [temp_name(1:end-3) 'greenberg_' temp_name(end-2:end)];
        dest = [out_path filesep out_name '.tif'];
        
        disp(['...............Start ''Greenberg_Reg'' for file ' temp_name ' .....................................']);
        [im_c dx_r dy_r E] = imreg_greenberg(im_s_dft, im_t, []);
        
        
        imwrite(uint16(im_c(:,:,1)), dest, 'tif', 'Compression', 'none', 'Description',im_descr, 'WriteMode', 'overwrite');
        for f=2:size(im_c,3)
            imwrite(uint16(im_c(:,:,f)), dest, 'tif', 'Compression', 'none', 'WriteMode', 'append');
        end
        
        fprintf('Reg data saved to %s\n',out_name);
    end
end
if matlabpool('size') > 0
    matlabpool close
end


% function parsave(fn, dx_r, dy_r, E, dft_shift)
%     save(fn, 'dx_r','dy_r','E','dft_shift');


    
    
% function [im_dft, shift] = dft_reg(im_s, im_tg, line_range)
% 
% % line_range, [startLine endLine], Only use a subportion of the image, specified by the range of
% %               lines (rows of the image matrix), to do registration and
% %               comput the shift. Th shift is then used to register the
% %               whole frame.
% 
% if nargin < 2 || isempty(im_tg)
%     % By default, use mean of the last 10 frames as the target image.
%     im_tg = mean(im_s(:,:,end-9:end),3);
% end
% if nargin < 3 || isempty(line_range)
%     row_num = 1 : size(im_tg,1);
% else
%     row_num = line_range(1) : line_range(2);
% end
% 
% for i=1:size(im_s,3);
%     output(:,i) = dftregistration(fft2(double(im_tg(row_num,:))),fft2(double(im_s(row_num, :,i))),1);
% end
% shift = output(3:4,:);
% padding = [0 0 0 0];
% im_dft = ImageTranslation_nx(im_s,shift,padding,0);
