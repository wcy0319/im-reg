%
% [im_c dx_r dy_r E] = imreg_greenberg(im_s, im_t, opt)
%
% Given two images, finds optimal x/y offset using the gradient descent
%  method from Greenberg et al., 2009:
%     "Automated correction of fast motion artifacts for two-photon
%      imaging of awake animals," D.S. Greenberg & J.N.D. Kerr, Journal of
%      Neuroscience Methods, 2009.
%      http://dx.doi.org/10.1016/j.jneumeth.2008.08.020
%  That algorithm is implemented in extern_greenberg_motioncorrect.m; see
%  it for details.  This is basically a wrapper that returns im_c,
%  the corrected image -- im_s is source image, im_t is 
%  target image (i.e., shifts source so it fits with target). Also returns
%  displacement in x (dx_r) and y (dy_r) that needs to be applied to source image
%  to match target.  dx_r > 0 and dy_r > 0 imply right and down movement, resp.
%  E is the error for each frame, measured as the absolute value of im_s-im-c.
%  Since this is a line-by-line correction algorithm, dx_r will be n x t where n
%  is the number of lines, t is the number of trials.  Same for dy.
%
% Note that im_s can be a stack; in this case, so will im_c.
%
% opt - structure containing parameters ; use structure to allow you to
%       vary the options without many function variables.  Right now blank.
%       norm_int: if 1, it will normalize so that inter-frame intensity for
%                 each *line* is same, using median as norm. factor.  
%                 Applies this to target and source. Default 1. This
%                 is useful if there are large luminance variations, and it
%                 basically assures that dark regions don't break the algo.
%
function [im_c dx_r dy_r E] = imreg_greenberg(im_s, im_t, opt)

  % --- opt check:
  if (length(opt) == 0) % defaults
	  opt.norm_int = 1;
	else % sanity checks - user does not have to pass all opts, so default unassigned ones
	  if (isfield(opt,'norm_int') == 0)
		  opt.norm_int = 1;
		end
	end

	% --- prelimes
	S_im_c = size(im_s);
	im_s_p = double(im_s);
	im_t_p = im_t;
    if length(S_im_c) < 3
        S_im_c(3) = 1;
    end
    %{
  % --- normalize line intensity?
  if (opt.norm_int)
	  disp(['imreg_greenberg::Normalizing image line intensities']);
		msf = median(reshape(im_s_p,[],1));
        % source image normalize
        for f=1:S_im_c(3)
            intensities = mean(im_s_p(:,:,f),2);
            sf = msf./intensities;
            for i=1:length(sf)
                im_s_p(i,:,f) = im_s_p(i,:,f)*sf(i);
            end
        end

        % target image normalize
        intensities = mean(im_t_p,2);
        sf = msf./intensities;
        for i=1:length(sf)
            im_t_p(i,:) = im_t_p(i,:)*sf(i);
        end
        
        disp(['imreg_greenberg::Done normalizingimage line intensities']);
	end
%}
	% --- run the algorithm: these is the setting structure it takes
	settings.move_thresh = 0.06; % .075 % .06 paper
	settings.max_iter = 50 ;% 50; % 120 paper too
	
    settings.scanlinesperparameter = 2; % 1 (2 paper) % NX 
	%settings.scanlinesperparameter = 4; % 1 (2 paper)
	settings.corr_thresh = 0.75; % 0.9; % .75 % .85 paper ... .9 works 4 me
    
    settings.pregauss = .75 ; %0.75
	settings.haltcorr = 0.99;% .95 ; .99 paper
	settings.dampcorr = 0.8; % .8 paper as well

    disp(sprintf('Relevant settings: scanlinesperparameter = %d, corr_thresh = %.2f',settings.scanlinesperparameter, settings.corr_thresh));
  % established deviations
	
	% experiments
%	settings.scanlinesperparameter = 4; % 1 (2 paper)
%	settings.move_thresh = 0.5; % .075 % .06 paper
%	settings.dampcorr = 0.7; % .8 paper as well
% settings.haltcorr = 0.85;% .95 ; .99 paper
%	settings.corr_thresh = 0.75; % .75 % .85 paper ... .9 works 4 me
	%settings.haltcorr = 0.95;% .95 ; .99 paper
    subrec = [];
  [p, iter_used, corr, failed, settings, xpixelposition, ypixelposition] = ...
    extern_greenberg_motioncorrect(im_s_p,im_t_p,subrec,settings,[]);

	% --- construct the corrected image -- p consists of [dx dy] -- based on size of p 
	nc_p = size(p,2);

  p2im_c = S_im_c(1)/((nc_p/2)-1); % how many img columns does each p entry span
	if (p2im_c < 1) ; dpsl = round(1/p2im_c) ; else ; dpsl =1 ; end % how many divisions per single scan line?
	if (p2im_c == 0) ; p2im_c = 1; end

  if (length(S_im_c) == 2) 
	  nframes = 1;
	else
	  nframes = S_im_c(3);
	end
  
	% --- loop over all frames and determine dx, dy
  for f=1:nframes
	  if (failed(f)) 
		  disp(['Warning: imreg_greenberg could not fit frame ' num2str(f)]);
		  continue ;
		end
	  dx = p(f,2:nc_p/2); % drop first term -- alternative would be to average
		                   % since you have n+1 pts for n (image) rows
		dy = p(f, (nc_p/2)+2:nc_p); % again drop the first of n+1 to get n entries

		% post processing -- discontinuities not allowed
    % fix by setting all dy to median
		fs = round(length(dx)/5); % filter size
		sdx = zeros(1,length(dx) + 2*fs);
		sdx = [dx(5)*ones(1,fs) dx dx(length(dx)-5)*ones(1,fs)]; % dont take edge values bc those are often strange
		sdx = medfilt1(sdx,fs); 
		dx = sdx(fs+1:length(sdx)-fs);

		fs = round(length(dy)/10); % filter size
		sdy = zeros(1,length(dy) + 2*fs);
		sdy = [dy(5)*ones(1,fs) dy dy(length(dy)-5)*ones(1,fs)]; % dont take edge values bc those are often strange
		sdy = medfilt1(sdy,fs); 
		dy = sdy(fs+1:length(sdy)-fs);

		dx_r(f,:) = dx;
		dy_r(f,:) = dy;
  end

	% --- call imreg_wrapup and get your final image
	wrap_opt.err_meth = 3; % correlation based
	wrap_opt.debug = 0;
  [im_c E] = imreg_wrapup (im_s, im_t, dx_r, dy_r, [], wrap_opt);

