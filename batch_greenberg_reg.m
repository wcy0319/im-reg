function batch_greenberg_reg(src_path, main_fname, out_path)

datafiles = dir([src_path filesep basename '*.tif']);
if ~exist(out_path,'var')
    out_path = src_path;
elseif ~isdir(out_path)
	mkdir(out_path);
end
% fnames = {};
% for i = 1:length(datafiles)
% 	if ~isempty(strfind(datafiles(i).name,'tif')) && strfind(datafiles(i).name, main_fname) == 1
% 		fnames = [fnames datafiles(i).name];
% 	end
% end
dft_shift = {};
for i=1:length(datafiles)
    % if isdir(datafiles(i).name)
      %  continue;
    % end
    filename = [src_path filesep datafiles(i).name];
    info = imfinfo(filename);
    if isfield(info(1), 'ImageDescription')
        im_descr = info(1).ImageDescription;
    else
        im_descr = '';
    end
    
    im_s = imread_multi(filename,'g');
    im_s_o = im_s;
    % Do dft translational registration before proceed
    [im_s, dft_shift{i}] = dft_reg(im_s_o);
    disp('Whole frame registered with ''dft_reg''');
    
    im_t = mean(im_s,3);
    disp(['Start ''greenberg_reg'' for file ' datafiles(i).name ' ...']);
    [im_c dx_r dy_r E] = imreg_greenberg(im_s, im_t, []);
    
    [pthstr,temp_name,ext] = fileparts(filename);
    out_name = [datafiles(i).name(1:end-3) 'greenberg_' datafiles(i).name(end-2:end)];
    dest = [out_path filesep out_name];
    
    imwrite(uint16(im_c(:,:,1)), [dest '.tif'], 'tif', 'Compression', 'none', 'Description',im_descr, 'WriteMode', 'overwrite');
    for f=2:size(im_c,3)
        imwrite(uint16(im_c(:,:,f)), [dest '.tif'], 'tif', 'Compression', 'none', 'WriteMode', 'append');
    end
    save([dest '_reginfo'], 'dx_r','dy_r','E','dft_shift');
    
end

function [im_dft, shift] = dft_reg(im_s)
im_tg = mean(im_s(:,:,1:10),3);
for i=1:size(im_s,3);
    output(:,i) = dftregistration(fft2(double(im_tg)),fft2(double(im_s(:,:,i))),1);
end
shift = output(3:4,:);
padding = [0 0 0 0];
im_dft = ImageTranslation_nx(im_s,shift,padding,0);
        
