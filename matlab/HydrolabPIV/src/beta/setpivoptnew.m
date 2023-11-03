function opt = setpivopt(varargin) 
% SETPIVOPT create a structure of option for use with <a 
% href="matlab:help normalpass">normalpass</a> 
% and <a href="matlab:help distortedpass">distportedpass</a>.
% 
% opt = SETPIVOPT('name1',value1,'name2',value2,...) creates an 
% options structure opt in which the named properties have the
% specified values. Any unspecified properties have default values. 
% Case is ignored for property names. 
% 
% opt = SETPIVOPT(oldopt,'name1',value1,'name2',value2,...) alters 
% an existing options structure oldopt.
%        
%     
% List of SETPIVOPT properties
% -------------------------------------------------------------------------
%
% opt = SETPIVOPT(..., 'subwindow', height, width, overlap, ...) sets
% the subwindow size and overlap. 
% 
% % Example: 24x24 pixel subwindow with 50% overlap
% opt = SETPIVOPT('subwindow',24,24,.50)
% 
% -------------------------------------------------------------------------
%
% opt = SETPIVOPT(..., 'range', [umin umax vmin vmax], ...) sets the
% search range for finding the peak within the difference measure.
% 
% Example: Set search range to 1/3 of 24x24 subwindow
% SETPIVOPT('subwindow',24,24,.50,'range',[-8 8 -8 8])
% -------------------------------------------------------------------------
%
% 
% opt = SETPIVOPT(..., 'measure', @mfun, mparam, @prefun, ...) sets
% the difference measure function
% 
% [dm,mm] = mfun(f1,m1,f2,m2,idh,idw,mparam)
% 
% Given two subwindows with masks f1/m1 and f2/m2 with sizes (i,j) and (k,l) 
% respectivly, @mfun should return a difference measure, dm, with zero 
% displacement peak located in the center and  the masked fraction, mm,
% both of size (i+k-1)x(j+l-1). The indices (idh,idw) show which part of 
% the difference measure is needed.
% 
% mparam provideds an additional argument to the @mfun function, see
% spesific mfun function for it usage.
%
% The function @prefun is applied to the difference measure before
% subpixel interpolation.
%
% For a set of fft based difference measures, including cross-correlation
% and minimum quadratic difference (mqd), see <a 
% href="matlab:help maskedcc">maskedcc</a>, <a 
% href="matlab:help maskedncc">maskedncc</a>, 
% <a href="matlab:help maskedmqd">maskedmqd</a>, <a
% href="matlab:help maskednmqd">maskednmqd</a>, <a
% href="matlab:help phasercorr">phasecorr</a>.
%
% In addition most of the measures is also implemented as straight 
% summation in Java, which is often faster when the search range is small, 
% since only dm(idh,idw) is computed. Mqd is provides as a special
% case of the more general Minkowski differences,
% see <a href="matlab:help maskedccj">maskedccj</a>, <a
% href="matlab:help maskednccj">maskednccj</a>, <a
% href="matlab:help maskedmdj">maskedmdj</a>, <a
% href="matlab:help maskednmdj">maskednmdj</a>.
% 
% % Example 1: Masked cross-correlation with padding
% % assuming peak is gaussian by using log as prefun.
% opt = SETPIVOPT('measure', @maskedcc, true, @log)
%
% % Example 2: Minimum absolute difference
% opt = SETPIVOPT('measure', @maskedmdj, 1, @(x) -x):
% -------------------------------------------------------------------------
%
% opt = SETPIVOPT(...,'window', @blackmanharris, ...) set a window function
% with positive weights, which is multiplied with the subwindows in both 
% direction before calculating the difference measures, to reduce spectral 
% leakage. See <a href="matlab:help window">window</a> for more information.
% -------------------------------------------------------------------------
%
% opt = SETPIVOPT(..., 'subpixel', @subpixelfun, ...) sets the
% subpixel estimator function
%
% [x0,delta,out] = subpixelfun(F)
% 
% Subpixelfun should return the integer posiston from the top corner, x0, 
% the subpixel correction, delta, and an indication whether the method 
% belives the result is valid (out=1) or an outlier (out=0).
% See also <a
% href="matlab:help subpixel3x3">subpixel3x3</a>, <a
% href="matlab:help subpixel3x3ls">subpixel3x3ls</a>, <a
% href="matlab:help subpixel3x3lm">subpixel3x3lm</a>, <a
% href="matlab:help subpixel5x5">subpixel5x5</a>, <a
% href="matlab:help subpixel3x2">subpixel3x2</a>, <a
% href="matlab:help subpixelnone">subpixelnone</a>.
% -------------------------------------------------------------------------
%
% opt = SETPIVOPT(..., 'ensemble',@ensfun, ...) sets the averaging
% function used when compting the ensemble difference measure.
% Defaults to @nanmean , for better robustness use @nanmeadian.
% Only used with ensemble PIV with single call to <a 
% href="matlab:help normalpass">normalpass</a> 
% and <a href="matlab:help distortedpass">distportedpass</a>.
% -------------------------------------------------------------------------
%
% opt = SETPIVOPT(..., 'snr', area, ...) sets up signal-to-noise ratio 
% estimation. The area can be set to 'sub', 'full' or 'none'.
% The signal is defined as the difference between the peak value and the noise,
% which is estimated as the median absolute deviation of the difference
% measure, over the search range ('sub') or complete difference measure
% ('full'). If area is set to 'none' the signal-to-noise ratio is not
% estimated. Default value is 'sub'.
% -------------------------------------------------------------------------
%
% opt = SETPIVOPT(..., 'savepeaks', true/false, ...) saves the peaks 
% found by the difference measure, defaults to false as this uses quite 
% a bit of memory. Needs to be set to true when using ensemble PIV with
% multiple calls to <a 
% href="matlab:help normalpass">normalpass</a> or <a
% href="matlab:help distortedpass">distportedpass</a>.
% -------------------------------------------------------------------------
%
% opt = SETPIVOPT(..., 'tableau', tableau, alpha, ...) sets up 
% the intergration of displacement field to undistort the image pair
% by providing a implicit Runge-Kutta Butcher tableau, see <a
% href="matlab:help rungekutta">rungekutta</a>.
% To only partially undistort the image pair, alpha can be set to 
% some positive value less than one.
% -------------------------------------------------------------------------
%
% opt = SETPIVOPT(..., 'iminterp', @iminterp2fun, iparam, ...) specifies
% interpolation method to be used in <a href="matlab:help distortedpass">distportedpass</a>. 
% Given the image/mask pair and the integrated displacements, @iminterp2fun
% should return the undistorted image/mask pair and be of the form
%
% [im1i,mask1i,im2i,mask2i] = iminterp2fun(im1,mask1,x1,y1,im2,mask2,x2,y2).
%
% Extrapolation should result in the mask to be zero.
% Lanczos resampling is slower but gives best result for small/normal size
% particles, while bspline works better for large particles, 
% see <a
% href="matlab:help iminterp2lanczosj">iminterp2lanczosj</a>, <a
% href="matlab:help iminterp2lanczosj">iminterp2bsplinej</a> and <a
% href="matlab:help iminterp2lanczosj">iminterp2matlabj</a>.
% -------------------------------------------------------------------------
%
% opt = SETPIVOPT(..., 'savedistortions',true/false, ...) saves 
% the integrated displacement used to undistort the image pair, 
% in <a href="matlab:help distortedpass">distportedpass</a>, defaults to false. 
% -------------------------------------------------------------------------
%
% See also normalpass, distortedpass


% Other option remove in future?
% opt = SETPIVOPT(..., 'rkfun',@rkfun, ...)
% 
% Add options
% Use extended subwindows setpivopt('subwindow',[H1 H2],[W1 W2],overlap)?
    
  if (nargin>0 && isstruct(varargin{1}))
    % Add options to previous structure
    opt = varargin{1};
    k=2;
  else
    % Set default options
    
    % Max displacements
    opt.Umin    =  -6;
    opt.Vmin    =  -6;
    opt.Umax    =   6;
    opt.Vmax    =   6;
  
    % Subwindow size
    opt.Width   =  32;
    opt.Height  =  32;
    opt.Overlap = .75;
    
    % Difference measure
    %opt.measure = @maskednccj;  % inconsistant speed?
    opt.measure = @maskedncc; % safest option with respect to speed?
    opt.window = @rectwin;
    opt.mparam = 1;    
    opt.prefun = @(x) x;      
    opt.ensfun = @nanmean;
    
    % Subpixel interpolation
    opt.subpixel = @subpixel3x3;    
    opt.savepeaks = false;
    opt.snr = 'sub';
      
    % Distorted pass integration (Runge-Kutta standard, 4. order)
    opt.rkfun = @rungekutta;
    opt.tableau = [0 0 0 0 0; 1/2 1/2 0 0 0; 1/2 0 1/2 0 0; 1 0 0 1 0; 0 1/6 1/3 1/3 1/6];   
    opt.alpha   =  1;     
    opt.iminterp = @iminterp2lanczosj; % @iminterp2bsplinej;  
    opt.iparam = 4;
    opt.savedistortions = false;
        
    k=1;
  end
  
  while(k < nargin)
    switch (lower(varargin{k}))  
     case 'range'
      tmp = varargin{k+1};
      opt.Umin = tmp(1);
      opt.Umax = tmp(2);
      opt.Vmin = tmp(3);
      opt.Vmax = tmp(4);
      k=k+2;  
     case 'subwindow'       
      opt.Width = varargin{k+1};
      opt.Height = varargin{k+2};            
      opt.Overlap = varargin{k+3};
      k=k+4;     
     case 'measure'
      opt.measure = varargin{k+1};
      opt.mparam = varargin{k+2};
      opt.prefun = varargin{k+3};
      k=k+4;      
     case 'window'
      opt.window = varargin{k+1};
      k=k+2;
     case 'ensemble'
      opt.ensfun = varargin{k+1};
      k=k+2;
     case 'subpixel'
      opt.subpixel = varargin{k+1};      
      k=k+2;            
     case 'snr'
      % 'none','sub','full'
      opt.snr = varargin{k+1};
      k=k+2;
     case 'savepeaks'
      opt.savepeaks = varargin{k+1};
      k=k+2;
     case 'iminterp'
      opt.iminterp = varargin{k+1};
      opt.iparam = varargin{k+2};
      k=k+3;
     case 'rkfun' % To be removed/changed in future?
      opt.rkfun = varargin{k+1};
      k=k+2;
     case 'tableau'
      switch(lower(varargin{k+1}))
       case 'euler'
        % Euler (1. order)
        opt.tableau = [0 0; 0 1];
       case 'modeuler'
        % Modified Euler (2. order)
        opt.tableau = [0 0 0; 1/2 1/2 0; 0 0 1];
       case 'heun'
        % Heun's method (2. order)
        opt.tableau = [0 0 0; 1 1 0; 0 1/2 1/2]; 
       case 'rk4std'
        % Runge-Kutta standard (4. order)
        opt.tableau = [0 0 0 0 0; 1/2 1/2 0 0 0; 1/2 0 1/2 0 0; 1 0 0 1 0; 0 1/6 1/3 1/3 1/6];  
       case 'rk43/8'
        % Runge-Kutta 3/8 rule (4. order)
        opt.tableau = [0 0 0 0 0; 1/3 1/3 0 0 0; 2/3 -1/3 1 0 0; 1 1 -1 1 0; 0 1/8 3/8 3/8 1/8];    
       otherwise
        opt.tableau = varargin{k+1};
      end      
      opt.alpha = varargin{k+2};
      k=k+3;        
     case 'savedistortions'
      opt.savedistortions = varargin{k+1};
      k=k+2;
     otherwise
      warning('Unknown option "%s"',varargin{k});  
      k=k+1;
      while(k < nargin && ~ischar(varargin{k}))
        k=k+1;
      end
    end
  end
  if(false)
  % Error checking of options    
  if(opt.Umax <= opt.Umin)
      error('Invalid search range, Umax must be larger than Umin.')
  end
  if(opt.Vmax <= opt.Vmin)
      error('Invalid search range, Vmax must be larger than Vmin.')      
  end
  if (opt.Umax >  (opt.Width-1) || opt.Umin < -(opt.Width-1) || ...
      opt.Vmax > (opt.Height-1) || opt.Vmin < -(opt.Height-1))
      error('Invalid search range, range too large compared to subwindow size.')
  end
  if (opt.Umax >  (opt.Width/2) || opt.Umin < -(opt.Width/2) || ...
      opt.Vmax > (opt.Height/2) || opt.Vmin < -(opt.Height/2))
      warning('It recommanded to set maximum search range to no more than half of subwindow size.')
  end
  dx = floor(opt.Width*(1-opt.Overlap));
  dy = floor(opt.Height*(1-opt.Overlap)); 
  if(dx==0 || dy==0)
      error('Too high overlap, or too small subwindows.');
  end
  end