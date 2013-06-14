% 
% S Peron Oct. 2009
%
% [im_c E] = imreg_wrapup (im_s, im_t, dx, dy, dtheta, opt)
%
% This is the common final-step script for the imreg_xxx methods.  It will return
%  and construct the corrected image from the source (which can be single or 
%  multi frame) given displacement vectors dx, dy (which can be one-per-frame
%  or many-per-frame -- i.e., rigid or "semi-rigid", where sub-segments are
%  rigid but the image itself is basically elastic). Note that pixels in im_c that
%  have no value due to image shift will be assigned -1.  *IT IS ASSUMED THAT YOUR
%  im_s SOURCE IMAGE HAS NO NEGATIVE VALUES* In addition, rotation thru angle
%  dtheta will be applied.  im_t is the target image, and is passed if you want to
%  have error metrics be measurable...
%
% In all cases (dx, dy, dtheta), specifying 1 value will mean that it is applied to
%  *all* frames; specifying the same number of values as frames will do the obvious
%  thing.  If you want to have more dx/dy (but NOT dtheta) than frames -- i.e., displace
%  frame chunks -- this will work provided that your number of frames < number of
%  divisions (i.e., this algo is kinda dumb and assumes that you have more frame subsegments
%  than frames; to change this, you will have to change the code - should be easy).
%
% For multi-lines per frame, var(f,:) implies variable for frame f (var is dx, dy).
%
% Rotation is always applied after translation, which may matter if you combine it 
%  with sub-frame translation.  imrotate is used, with bicubic interpolation.  The
%  center of rotation is the image center - this is non negotiable!
%
% E is a per-frame error metric. This can be used for post-processing correction, 
%  as well as optimization.
%
% opt specifies any additional options and is a structure (* denotes default; pass as
%   [] to use defaults)
%   opt.err_meth:  1: raw difference
%                  2: difference only for matched regiond (i.e., pixels where im_c
%                    dne -1)
%                 *3: normalized cross-correlation
%   opt.debug: 0: nothing - defaukt
%              1: show wait bars
%              2: show movie comparing target img to src
%   opt.post_proc: post-processing technique(s) to apply ; set flags to true to use.
%                  Consists of a vector with 1 if you want it used; 0 if not.  Element
%                  number v. function: (*: only available for cases where you 
%                  have a dx/dy for mutliple lines/line blocks per frame)
%                  1: perform post-hoc median filter* (size = 1/10th frame height)
%                  2: post-hoc edge correction* (corrects top and bottom)
%                  3: adaptive correction based on line-by-line correlation* -- 
%                     lines with correlations below median will have dx/dy set
%                     to proximal lines with acceptable correlations
%                  4: best-of*: if enabled, will try #3 above, 'rigid' displacement
%                     with median corrected displacement (post #3) and median non-corr
%                     displacement -- and keep the one with the lowest error/highest
%                     corr
%   opt.lE: for opt.post_proc 3 to work, you must pass line-by-line errors so that
%           adaptive rescaling of dx & dy can be performed.  Same size as dx, dy.
%
function [im_c E] = imreg_wrapup (im_s, im_t, dx, dy, dtheta, opt)
  % prelims:
  sim = size(im_s);
  sdx = size(dx);
	if (length(sdx) == 1) ; sdx = [1 sdx]; end
	if (length(sim) > 2) 
	  nframes = sim(3);
	else
	  nframes = 1;
	end

	% is this rigid or piecewise?
	piecewise = 1; % 0: rigid ; 1: piecewise
	if ((sdx(2) == 1 & nframes > 1) | (sdx(2) == 1 & nframes == 1))
	  piecewise = 0;
	end

  % opt check:
  if (length(opt) == 0) % defaults
	  opt.err_meth = 3;
		opt.debug = 0;
		opt.post_proc = [0 0 0 0];
		opt.lE = [];
	else % sanity checks
	  if (isfield(opt,'err_meth') == 0)
		  opt.err_meth = 3;
		end
	  if (isfield(opt,'debug') == 0)
		  opt.debug = 0;
		end
		if (isfield(opt,'post_proc') == 0)
		  opt.post_proc = [0 0 0 0];
		elseif (~piecewise)
		  opt.post_proc([1 2 3 4]) = 0; % disable piecewise-only post-processing steps
		end
		if (isfield(opt,'lE') == 0)
		  opt.lE = [];
		end
	end
	if (opt.post_proc(4)) ; opt.post_proc(3) = 1; end % can't have one without the other

	if (length(im_t) == 0) % no target image? then no error
	  opt.err_meth = 0; 
	end

  % waitbar
	if (opt.debug  >= 1) ; wb = waitbar(0, 'imreg-wrapup::post-processing'); end

  % --- im_c construction
	im_c = -1*ones(sim);

	% based on size of dx, do your im_c construction -- assume dx and dy are same
	%  size
	if ( ~piecewise)
	     
    % SINGLE dx or dy? cheap but effective ;-)
		if (length(sdx) == 1 & sdx(1) == 1)
		  dx = dx*ones(1,nframes);
		  dy = dy*ones(1,nframes);
		end

    % loop over frames
		for f=1:nframes
            im_c(:,:,f) = -im_c(:,:,f).*mean(mean(im_s(:,:,f)));  % added by NX, assign a mean value to blank area for rigid correction.
			if (opt.debug  >= 1) ;  waitbar(f/nframes,wb, 'imreg-wrapup::alignment'); end
			if (dy(f) > 0)
				y1_s = 1;
				y2_s = sim(1)-dy(f);
				y1_c = dy(f)+1;
				y2_c = sim(1);
			else
				y1_s = -1*dy(f) + 1;
				y2_s = sim(1);
				y1_c = 1;
				y2_c = sim(1)+dy(f);
			end
			if (dx(f) > 0)
				x1_s = 1;
				x2_s = sim(2)-dx(f);
				x1_c = dx(f)+1;
				x2_c = sim(2);
			else
				x1_s = -1*dx(f) + 1;
				x2_s = sim(2);
				x1_c = 1;
				x2_c = sim(2)+dx(f);
			end
			im_c(y1_c:y2_c,x1_c:x2_c,f) = im_s(y1_s:y2_s,x1_s:x2_s,f);
        end

	else % here we assume the only alternative -- piecewise displacement 

    % preliminary variable calculation
		dpc = sim(1)/sdx(2); % how many img columns does each dx entry span
		if (dpc < 1) ; dpsl = round(1/dpc) ; else ; dpsl =1 ; end % how many divisions per single scan line?
		if (dpc == 0) ; dpc = 1; end

    % --- keep originals around
    dx_o = dx;
		dy_o = dy;

    % --- median filter? if so, first 
    if (opt.post_proc(1))
			if (opt.debug  >= 1) ;  waitbar(0.25,wb, 'imreg-wrapup::median filtering'); end
			mfs = round(sdx(2)/10); % 1/10th of height
			dxt = medfilt1(reshape(dx',[],1),mfs);
			dyt = medfilt1(reshape(dy',[],1),mfs);
			
			% median filtered dx
			dx_mf = zeros(sdx(1), sdx(2));
			dy_mf = zeros(sdx(1), sdx(2));
			
			for f=1:nframes
				si = (f-1)*sdx(2) + 1;
				ei = si + sdx(2) -1;
				dx_mf(f,:) = dxt(si:ei);
				dy_mf(f,:) = dyt(si:ei);
			end

			dx = dx_mf;
			dy = dy_mf;
		end

		% edge-correction? if so, do it (i.e., if previous frame says bottom/top will be out of 
		%  frame, skip)

		% adaptive?  if so, do adaptive correction of dx/dy based on single-line correlation 
    if (opt.post_proc(3))
			if (opt.debug  >= 1) ;  waitbar(0.5, wb,'imreg-wrapup::adaptive correction'); end

		  dx_pre_c = dx;
		  dy_pre_c = dy;
      % GOOD FOR ME: [dx_c dy_c] = correct_displacement_vectors(dx, dy, opt.lE, 0.8, 1, 10, opt.debug);
      % GOOD FOR NLX data: [dx_c dy_c] = correct_displacement_vectors(dx, dy, opt.lE, 0.9, 1, 10, opt.debug);
      % GOOD FOR NLX (best?) : [dx_c dy_c] = correct_displacement_vectors(dx, dy, opt.lE, .9, 1, 15, opt.debug);
      %[dx_c dy_c] = correct_displacement_vectors(dx, dy, opt.lE, 0.8, 1, 10, opt.debug);
      [dx_c dy_c] = correct_displacement_vectors(dx, dy, opt.lE, .9, 1, 20, opt.debug);

			dx = dx_c;
			dy = dy_c;
%	disp('saving to ~/Desktop/imreg_wrapup.mat'); save('~/Desktop/imreg_wrapup.mat', 'dx_o', 'dy_o', 'dx_mf', 'dy_mf', 'dx_c', 'dy_c', 'opt');
		end

		% best of? if so, get im_c from *multiple* approaches
		if (opt.post_proc(4))
		  % first, we do two more alpha values for adaptive...
			dx_c_p8 = dx_c;
			dy_c_p8 = dy_c;
			if (opt.debug  >= 1) ;  waitbar(0.7,wb, 'imreg-wrapup::adaptive correction 2'); end
      [dx_c_p9 dy_c_p9] = correct_displacement_vectors(dx_pre_c, dy_pre_c, opt.lE, 0.9, 1, 10, opt.debug);
			if (opt.debug  >= 1) ;  waitbar(0.9,wb, 'imreg-wrapup::adaptive correction 3'); end
      [dx_c_1 dy_c_1] = correct_displacement_vectors(dx_pre_c, dy_pre_c, opt.lE, 1, 1, 10, opt.debug);

			% now build
			im_c_p8 = build_piecewise_im_c(im_s, dx_c_p8, dy_c_p8, sim, sdx, nframes, dpc, dpsl);
			im_c_p9 = build_piecewise_im_c(im_s, dx_c_p9, dy_c_p9, sim, sdx, nframes, dpc, dpsl);
			im_c_1 = build_piecewise_im_c(im_s, dx_c_1, dy_c_1, sim, sdx, nframes, dpc, dpsl);
			im_c_o = build_piecewise_im_c(im_s, dx_o, dy_o, sim, sdx, nframes, dpc, dpsl);
			im_c_mf = build_piecewise_im_c(im_s, dx_mf, dy_mf, sim, sdx, nframes, dpc, dpsl);

			% for each one, get errors per-frame, and take best frame for each one -- use corr error ALWAYS
			if (opt.debug  >= 1) ;  waitbar(1, wb, 'imreg-wrapup::selecting best per frame'); end
			E_base = calculate_error(3, im_s, im_t, 1:nframes);
			E_p8 = calculate_error(3, im_c_p8, im_t, 1:nframes);
			E_p9 = calculate_error(3, im_c_p9, im_t, 1:nframes);
			E_1 = calculate_error(3, im_c_1, im_t, 1:nframes);
			E_o = calculate_error(3, im_c_o, im_t, 1:nframes);
			E_mf = calculate_error(3, im_c_mf, im_t, 1:nframes);
			im_c = zeros(size(im_s));
      for f=1:nframes
			  [irr mi] = max([E_base(f) E_p8(f) E_p9(f) E_1(f) E_o(f) E_mf(f)]);
				if (mi == 1) % source image is best
				  im_c(:,:,f) = im_s(:,:,f);
				elseif(mi == 2) % 0.8 adaptive
				  im_c(:,:,f) = im_c_p8(:,:,f);
				elseif(mi == 3) % 0.9 adaptive
				  im_c(:,:,f) = im_c_p9(:,:,f);
				elseif(mi == 4) % 1 adaptive
				  im_c(:,:,f) = im_c_1(:,:,f);
				elseif(mi == 5) % just correction
				  im_c(:,:,f) = im_c_o(:,:,f);
				elseif(mi == 6) % just median filter
				  im_c(:,:,f) = im_c_mf(:,:,f);
				end
			end
	  else 
			im_c = build_piecewise_im_c(im_s, dx, dy, sim, sdx, nframes, dpc, dpsl);
		end


    % --- post-correction median filter
    if (opt.post_proc(1))
			if (opt.debug  >= 1) ;  waitbar(0.25,wb, 'imreg-wrapup::median filtering'); end
			mfs = round(sdx(2)/10); % 1/10th of height
			dxt = medfilt1(reshape(dx',[],1),mfs);
			dyt = medfilt1(reshape(dy',[],1),mfs);
			
			% median filtered dx
			dx_mf = zeros(sdx(1), sdx(2));
			dy_mf = zeros(sdx(1), sdx(2));
			
			for f=1:nframes
				si = (f-1)*sdx(2) + 1;
				ei = si + sdx(2) -1;
				dx_mf(f,:) = dxt(si:ei);
				dy_mf(f,:) = dyt(si:ei);
			end

			dx = dx_mf;
			dy = dy_mf;
		end
	end

	% --- apply rotations
  if (length(dtheta) > 0)
	  if (length(dtheta) == 1)
		  dtheta = dtheta*ones(1,nframes);
		end
		for f=1:nframes
			if (opt.debug  >= 1) ;  waitbar(f/nframes,wb, 'imreg-wrapup::rotation'); end
			im_c(:,:,f) = imrotate(im_c(:,:,f),dtheta(f), 'bicubic','crop');
		end
	end

	% --- error measure
	if (opt.debug  >= 1) ;  waitbar(1, wb, 'imreg-wrapup::measuring error'); end
	E = calculate_error(opt.err_meth, im_c, im_t, 1:nframes);

	% --- plot movie? for debugging
	if (opt.debug  >= 1) ;  delete(wb); end
	if (opt.debug >= 2)
    plot_err_movie( im_s, im_t, im_c, dx, dy, sim, nframes, opt.debug)
	end
%
% This will use dx and dy on im_s to build a corrected image, assuming *multi*
%  lines per frame (piecewise rigid).  sim is size(im), sdx is size(dx). dpc
%  and dpsl are the number of image columns spanned by each dx entry and dpsl
%  is the number of divisions per scaline (stems from variable "p", which is the
%  vector containing dx in an optimization algorithm - p is the optimized variable
%  in classic optimization theory).
%
function im_c = build_piecewise_im_c(im_s, dx, dy, sim, sdx, nframes, dpc, dpsl)
	for f=1:nframes
		x_i = 1;
		for i=1:sdx(2)

			% The rows -- or y coordinates
			if (dpsl == 1)
				y_s = (i-1)*dpc+1:i*dpc;
			else
				y_s = ceil(i/dpsl);
			end
			y_c = round(y_s + dy(f,i));

			y_bad = find (y_c < 1 | y_c > sim(1));
			y_good = setdiff(1:length(y_s),y_bad);
			y_s = y_s(y_good);
			y_c = y_c(y_good);

			% The columns -- or x coordinates -- ASSUME integral number of lines
			x_s = (x_i-1)*(sim(2)/dpsl)+1:(x_i)*(sim(2)/dpsl);
			x_c = round(x_s + dx(f,i));

			x_i = x_i + 1; 
			if (x_i > dpsl) ; x_i = 1; end

			x_bad = find (x_c < 1 | x_c > sim(2));
			x_good = setdiff(1:length(x_s),x_bad);
			x_s = x_s(x_good);
			x_c = x_c(x_good);

			
			% build image subsection if there is something there -- at edges, may have nothing to do
			if (length(y_s) > 0 & length(x_s) > 0)
				im_c(y_c,x_c,f) = im_s(y_s,x_s,f);
			end
		end
	end
 
%
% This function computes the error for all frames specified
%  err_meth: 1: raw difference (L1) ; 2: diff only where > 0 
%            3: median-normalized correlation
%
function	E = calculate_error(err_meth, im_c, im_t, frames);
  E = zeros(length(frames),1);
	for f=frames
    if (err_meth == 1) % raw diff
			E(f) = sum(sum(im_t - im_c(:,:,f)));
    elseif (err_meth == 2) % diff where > 0
			imsl = reshape(im_t,1,[]);
			imcl = reshape(im_c(:,:,f),1,[]);
			val = find(imcl ~= -1);
			E(f) = sum(abs(imcl(val)-imtl(val)));
    elseif (err_meth == 3) % correlation normalized to MEDIAN
			imsl = reshape(im_t,1,[]);
			imcl = double(reshape(im_c(:,:,f),1,[]));
			val = find(imcl ~= -1);
			imsl = imsl(val)/median(imsl(val));
			imcl = imcl(val)/median(imcl(val));
      R = corrcoef(imsl,imcl);
			E(f) = R(1,2);
	  end
	end

%
% This will play difference images as a movie
%
% debug: 2 - show
%
function plot_err_movie( im_s, im_t, im_c, dx, dy, sim, nframes, debug)
  % invoke figure, defaults
	im1_rgb = zeros(sim(1),sim(2), 3);
	im2_rgb = zeros(sim(1),sim(2), 3);
	im_s = im_s/max(max(max(im_s)));
	im_t = im_t/max(max(max(im_t)));
	dim_c = im_c/max(max(max(im_c)));

	im1_rgb(:,:,1) = im_t;

	figure;
	sp1 = subplot(2,2,1);
	sp2 = subplot(2,2,2);

	for f=1:nframes
		% plot
		subplot(sp1);
		im2_rgb(:,:,2) = dim_c(:,:,f);
		imshow(im1_rgb + im2_rgb, 'Border','tight');
%		axis square;

		imsl = reshape(im_t,1,[]);
		imcl = reshape(im_c(:,:,f),1,[]);
		val = find(imcl ~= -1);
		imsl = imsl(val)/median(imsl(val));
		imcl = imcl(val)/median(imcl(val));
		R = corrcoef(imsl,imcl);
		Err = R(1,2);
		title(['f: ' num2str(f) ' corr: ' num2str(Err)]);

		subplot(sp2);
		im2_rgb(:,:,2) = im_s(:,:,f);
		imshow(im1_rgb + im2_rgb, 'Border', 'tight');
%		axis square;
		imsl = reshape(im_t,1,[]);
		imcl = reshape(im_s(:,:,f),1,[]);
		val = find(imcl ~= -1);
		imsl = imsl(val)/median(imsl(val));
		imcl = imcl(val)/median(imcl(val));
		R = corrcoef(imsl,imcl);
		Err = R(1,2);
		title(['f: ' num2str(f) ' corr: ' num2str(Err)]);

		% plot dx and dy vectors if needbe
		if (size(dx,1) >= f)
			subplot(2,2,3) ;
			plot(dx(f,:), 'b');
			title('dx');
			subplot(2,2,4) ;
			plot(dy(f,:));
			title('dy');
		end
		pause;
	end

%
% This performs adaptive correction on a dx/dy series given an error
%  Returns corrected dx and dy displacement vectors. lE is line-by-line
%  error vector.  The following parameters are important:
%    alpha: lE (errors -- correlations actually) below median(lE)*alpha are "bad"
%    beta: medial(lE)*beta is good
%    gbs: good block size (see below)
%  If you are "bad", you are replaced by the mean value of the first gbs good points
%    and the preceding good gbs points.  
%  debug: if 1, waitbar
%
function [dx_r dy_r] = correct_displacement_vectors(dx, dy, lE, alpha, beta, gbs, debug)
  nframes = size(dx,1);
	nlines = size(dx,2);

	lE = reshape(lE',[],1); % error throughout, line-by-line
	dx_c = reshape(dx',[],1); % dx and dy -- Corrected
	dy_c = reshape(dy',[],1);
	mE = median(lE);
%mE
%gbs
	% go line-by-line
	lE_good = find(lE >= mE*beta);
	Lg = length(lE_good);
	lE_bad = find(lE < mE*alpha);
	Lb = length(lE_bad);
	last_nxt = [];
	if (debug >= 1) ; wb = waitbar(0,'Applying adaptive correction ...'); end
	for b=1:length(lE_bad)
		if (debug >= 1) ; wb = waitbar(b/length(lE_bad), wb, 'Applying adaptive correction ...'); end
	  
		% for all bad points, find good points before and after
		prev = max(find(lE_good < lE_bad(b)));
		nxt = min(find(lE_good > lE_bad(b)));

		% blockiness
		if (gbs > 1)
			found = 0;	  
			% next block:
			while(found == 0 & length(nxt) > 0)
				if (nxt + gbs > Lg) 
					found = -1;
				else
					desired = lE_good(nxt)+(0:gbs-1);
					if (length(intersect(desired,lE_good)) == length(desired))
						found = 1;

%fidx = 1+floor(lE_good(nxt)/(nlines));
%lidx = lE_good(nxt)-((fidx-1)*(nlines));
%lidxb = lE_bad(b)-((fidx-1)*(nlines));
%disp(['f: ' num2str(fidx) ' from l: ' num2str(lidxb) ' to : ' num2str(lidx)]);

						nxt = nxt:nxt+gbs-1;
						% the current next will be the prev on the following run
						if (length(last_nxt) > 0) 
							prev = last_nxt;
						end
						last_nxt = nxt;
					else
						nxt = nxt+1;
					end
				end
			end
		end

		% make sure they exist -- at end and beginning you only have 1 
		Sdx = []; % sum of 2 good pts
		Sdy = []; % sum of 2 good pts
		if (length(prev) > 0)
			Sdx = [Sdx  dx_c(lE_good(prev))'];
			Sdy = [Sdy  dy_c(lE_good(prev))'];
		end
		if (length(nxt) > 0)
			Sdx = [Sdx  dx_c(lE_good(nxt))'];
			Sdy = [Sdy  dy_c(lE_good(nxt))'];
		end

		% and correct . . . 
		dx_c(lE_bad(b)) = mean(Sdx);
		dy_c(lE_bad(b)) = mean(Sdy);
	end
			
  if (debug >= 1) ; delete(wb) ; end

	% and construct reshaped vector
	dx_r = zeros(nframes,nlines);
	dy_r = zeros(nframes,nlines);
	
	for f=1:nframes
		si = (f-1)*nlines + 1;
		ei = si + nlines -1;
		dx_r(f,:) = dx_c(si:ei);
		dy_r(f,:) = dy_c(si:ei);
	end
