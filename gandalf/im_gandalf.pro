;######################################################################
;
; Copyright (C) 2006, M. Sarzi, J. Falcon-Barroso & R.F. Peletier
; E-mail: sarzi@star.herts.ac.uk
; E-mail: jfalcon@rssd.esa.int
; E-mail: peletier@astro.rug.nl
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;######################################################################
;+
; NAME:
;   GANDALF 
;   (Gas AND Absorption Line Fitting, formerly known as PPXF_GAS)
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; PRO GANDALF, templates, galaxy, noise, velScale, sol, emission_setup, l0_gal, lstep_gal,       $
;              GOODPIXELS=goodPixels, DEGREE=degree, MDEGREE=mdegree, INT_DISP=int_disp,            $
;              BESTFIT=bestFit, EMISSION_TEMPLATES=emission_templates, WEIGHTS=weights, ERROR=esol, $
;              PLOT=plot, QUIET=quiet, LOG10=log10, REDDENING=reddening, L0_TEMPL=l0_templ,         $
;              FOR_ERRORS=for_errors
;
; INPUT PARAMETERS:
;
;   TEMPLATES: vector containing the spectrum of a single template
;       star or array of dimensions nPixels_Templates x nTemplates
;       containing different templates to be optimized during the fit
;       of the kinematics.  nPixels_Templates has to be >= the number
;       of pixels sampling the galaxy spectrum nPixels_Galaxy.
;   GALAXY: vector containing the spectrum of the galaxy to be measured.
;
;       Both star and the galaxy spectra HAVE to be logarithmically
;       rebinned on a natural ln-base (or in log10, see LOG10 keyword).
;
;   NOISE: vector containing the 1*sigma error of the emission
;       spectrum. If this is not available an array of constant unity
;       values should be passed.
;   VELSCALE: velocity scale of the spectra in km/s per pixel. It has
;       to be the same for both the galaxy and the template spectra.
;   SOL: on INPUT it must be a vector containing the stellar
;       kinematics needed to convolve the input stellar spectrum or
;       template library. On OUTPUT it will contain the results of the
;       gas fit with weights assigned to the multiplicative polynomial
;       appended.
;   EMISSION_SETUP: A structure containing the index, the name, and
;       the wavelength of the fitted emission lines. It must also
;       contain the starting/input values for the line amplitudes,
;       velocities, and widths of the lines, as well as keywords
;       specifying whether each line is a part of a doublet, and
;       whether its position and width fit freely, to be hold at its
;       input values, or to be tied that of another line
;   L0_GAL: the ln-lambda value corresponding to the starting pixel in
;       the galaxy spectrum.
;   L0_STEP: the ln-lambda step of corresponding to the pixels
;       sampling the galaxy spectrum. 
;   
;       Wavelengths are assumed to be in Angstrom
;
; KEYWORDS:
;
;   DEGREE: degree of the Legendre polynomial used for correcting
;       the template continuum shape during the fit (default: -1).
;       This correction is ADDITIVE, and should NOT be
;       used. Eventually this keyword will be remouved.
;   MDEGREE: degree of the Legendre polynomial used for correcting
;       the template continuum shape during the fit (default: 0).
;       This correction is MULTIPLICATIVE. 
;   GOODPIXELS: integer vector containing the indices of the pixels in
;       the galaxy spectrum (in increasing order) that will be
;       included in the fit. IMPORTANT: in all likely situations this
;       keyword *has* to be specified.
;   INT_DISP: allows to input the instrumental dispersion (sigma) in km/s
;       of the instrument used. This is used to set a lower limit on the
;       velocity dispersion of the emission lines measured.
;       It is *NOT* used for anything else. The default value is 0.
;   LOG10: allows to deal with data that have been logarithmically
;       rebinned in lambda, using a base 10
;   REDDENING: allows to include in the fit the effect of reddening by
;       dust, by specifying a single E(B-V) guess for extinction, or a
;       two-element array of E(B-V) guesses. A single guess will
;       trigger the use of a single-screen dust model, affecting the
;       entire spectrum, whereas passing a two-elements array of
;       guesses will add a second dust component affecting only the
;       nebular fluxes. This second option is recommended only when
;       recombination lines are clearly detected, so that a
;       temperature-dependent prior (e.g. the decrement of the Balmer
;       lines) on their relative intensity can be use to constrain
;       such a second, internal component. 
;   L0_TEMPL: the ln-lambda value corresponding to the starting pixel
;       in the template spectra. Needed if using the REDDENING keyword
;       and the DUST_CALZETTI function. 
;   FOR_ERRORS: A keyword specifying whether we wish errors to be
;       properly estimated

;   /QUIET: set this keyword to mute the output.
;   /PLOT: set this keyword to plot the best fitting solution at the end of the fit.
;
; OUTPUT PARAMETERS:
;
;   SOL: A 4xNlines vector containing the values of the best fitting
;       flux, amplitude, Velocity, and velocity dispersion for each
;       emission line. If the keyword INT_DISP has been set the velocity
;       dispersion is already the intrinsic one, otherwise it will be the
;       observed width of the line.
;   BESTFIT: a named variable to receive a vector containing sum of
;       best fitting stellar (convolved by the input LOSVD) and
;       emission-line Gaussian templates. This model include bending
;       of the stellar templates by the best fitting multiplicative
;       polynomials, and of additive polynomials, if specified.
;   EMISSION_TEMPLATES: a named variable to receive an
;       [nPixels_Galaxy,Nlines] array containing the best fitting emission-lines
;       templates. The emission spectrum can be obtained from this
;       array by simply doing total(emission_templates,2).
;   WEIGHTS: a named variable to receive the value of the weights by which
;       each template was multiplied to best fit the galaxy spectrum.
;   ERROR: a named variable that will contain a vector of formal
;       errors (1 sigma) for the parameters in the output vector
;       SOL. 

;       If the FOR_ERRORS keyword IS NOT specified no errors for the
;       line amplitudes and fluxes will be returned and the
;       uncertainties on the position and width of the lines should be
;       regarded only as order of magnitude estimates. 
;       If the FOR_ERRORS keyword IS specified, a second fully
;       non-linear fit of all the emission-line parameter will provide
;       correct estimates for the errors on all emission-line
;       parameters.
;       Still, Keep in mind that these errors are meaningless unless
;       Chi^2/DOF~1. If a constant noise spectrum is provided, the
;       formal errors will be automatically rescaled under the
;       assumption the model is a good representation of the data.
;
; REQUIRED ROUTINES:
;
;       RANGE: by M. Cappellari from http://www.strw.leidenuniv.nl/~mcappell/idl/
;       BVLS:  by M. Cappellari from http://www.strw.leidenuniv.nl/~mcappell/idl/
;       MPFIT: by C.B. Markwardt from http://astrog.physics.wisc.edu/~craigm/idl/
;       ROBUST_SIGMA: by H. Freudenreich from http://idlastro.gfsc.nasa.gov/
;
; MODIFICATION HISTORY: 
;
; V1.0 Initially conceived by J. Falcon-Barroso & R.F. Peletier to fit
;      the Hb, [OIII], and [NI] lines in the SAURON spectra by
;      including Gaussian templates in the IDL routine PPXF of
;      Cappellari & Emsellem (2004, PASP, 116, 138) to derive the
;      stellar kinematics, it was then further generalised by M. Sarzi
;      to include any number emission lines and to allow different
;      sort of fits. The present version allow to derive emission-line
;      fluxes and kinematics fully consistent with those presented in
;      Sarzi et al. (2005, MNRAS, 366, 1151) and Falcon-Barroso et al.
;      (2006 MNRAS, 369, 529).
;
; V1.1 Fixed problem occurring in the absence of doublets. (Thanks to
;      JFB, la Palma, 15/12/06)
;
; V1.2 Added the possibility to tie only the position (v#) or only the
;      width (s#) of several lines in the emission-line setup file.
;
; V1.3 Added keyword LOG10 to deal with object and template spectra
;      that have been logarithmically rebinned in wavelength in
;      base-10.
;
; V1.4 Added keyword REDDENING to allow the use of reddening by dyst
;      to adjust of the stellar continuum and the emission line
;      fluxes. The keyword REDDENING allows to use either a single
;      dust component, affecting both the stellar continuum and the
;      emission-line fluxes, or additionally a second dust component
;      affecting only the emission-line templates. The importance of
;      reddening by dust is constraint by the observed decrement of
;      the Balmer-lines emission (or any other recombination series)
;      and, to some extent, by the shape of the stellar continuum. See
;      specific description of the REDDENING keyword for
;      details. Reddening by dust CANNOT be used in conjunction with a
;      multiplicative polynomial adjustment.
;
;      The function DUST_CALZETTI has been included to use a Calzetti
;      et al. (2000, ApJ, 533, 682) prescription for reddening by dust.
;
; V1.5 Now use    l_obs = l_rest * exp(v_sys/c) 
;      instead of l_obs = l_rest * (1+v_sys/c)
;
; V1.6 Major reajustment and resctructuring to include error
;      calculation using the keyword FOR_ERRORS
;
;
; J. Moustakas, 2009-Nov-19, UCSD - lots of small tweaks; major
;   changes include: (1) do not apply the reddening to the
;   emission-line templates


;------------------------------------------------------------------------------------
FUNCTION BVLSN_Solve_pxf, A, b, degree, $
                          FOR_ERRORS=for_errors, NLINES=nlines
compile_opt idl2, hidden

; No need to enforce positivity constraints if fitting one
; single template: use faster SVD solution instead of BVLS.
;
s = size(a)
IF s[2] EQ degree+2 THEN BEGIN ; Fitting single template
    SVDC, A, w, u, v, /COLUMN, /DOUBLE
    small = WHERE(w LT MAX(w)*1d-10,m)
    IF m NE 0 THEN w[small] = 0d
    soluz = SVSOL(u, w, v, b, /COLUMN, /DOUBLE)
ENDIF ELSE BEGIN               ; Fitting multiple templates
    BND = DBLARR(2,s[2],/NOZERO)
    mx = (MACHAR()).XMAX
    IF (degree GE 0) THEN BND[0,0:degree] = -mx   ; No constraints on the Legendre polynomials
    BND[0,degree+1:*] = 0d  	    	    	  ; Positivity constraints on the templates
    BND[1,*] = mx
    ; If we are fitting with MPFIT also the amplitudes of the emission lines, 
    ; then fix the weights for these components, by setting them to 1.0
    if keyword_set(for_errors) and keyword_set(nlines) then begin
        for i=0,nlines-1 do begin
            BND[0,s[2]-1-i] = 1.0d
            BND[1,s[2]-1-i] = 1.0d
        endfor
    endif
    BVLS, A, B, BND, soluz
ENDELSE

return, soluz
END

;------------------------------------------------------------------------------------
FUNCTION CREATE_GAUSSN, X=x, PARS=pars,INT_DISP_PIX=int_disp_pix
; The instrumental resolution is supposed to be in sigma and in pixels
; at this stage

npars=n_elements(pars)
npix = n_elements(x)
y=fltarr(npix)
FOR i=0, npars-3, 3 DO BEGIN
    pars[i+2] = sqrt(pars[i+2]*pars[i+2] + (int_disp_pix)^2)
    w=(findgen(n_elements(x))-pars[i+1])/pars[i+2]
    tmp=pars[i]*EXP(-w^2/2.)
    y=y+tmp
ENDFOR

return, y
END

;------------------------------------------------------------------------------------
FUNCTION CONVOLVE_TEMPLATES, templates,kinstars,velscale

vel   = kinstars[0]/velscale   ; in pixels
sigma = kinstars[1]/velscale   ; in pixels
dx = ceil(abs(vel)+4*sigma)    ; Sample the Gaussian and GH at least to vel+4*sigma
x  = range(dx,-dx)             ; Evaluate the Gaussian using steps of 1 pixel.
w  = (x - vel)/sigma
w2 = w^2
losvd = exp(-0.5d*w2)/(sqrt(2d*!dpi)*sigma) ; Normalized total(Gaussian)=1


; Hermite polynomials as in van der Marel & Franx (1993).
; Coefficients are given e.g. in Appendix C of Cappellari et al. (2002)
nkins = n_elements(kinstars)
IF nkins GT 2 THEN BEGIN
    poly = 1d + kinstars[2]/Sqrt(3d)*(w*(2d*w2-3d)) $       ; H3
              + kinstars[3]/Sqrt(24d)*(w2*(4d*w2-12d)+3d)   ; H4
    IF nkins EQ 6 THEN $
        poly = poly + kinstars[4]/Sqrt(60d)*(w*(w2*(4d*w2-20d)+15d)) $      ; H5
                    + kinstars[5]/Sqrt(720d)*(w2*(w2*(8d*w2-60d)+90d)-15d)  ; H6
    losvd = losvd*poly
ENDIF

s = size(templates)
ctemplates = dblarr(s[1],s[2])
IF s[0] EQ 2 THEN ntemp = s[2] ELSE ntemp = 1   ; Number of template spectra
FOR j=0,ntemp-1  DO BEGIN
        ctemplates[*,j] = convol(templates[*,j],losvd,/EDGE_TRUNCATE) ;
ENDFOR

return,ctemplates
END

FUNCTION DUST_CALZETTI,l0_gal,lstep_gal,npix,ebv,vstar,LOG10=log10
; This procedure uses the dust model of Calzetti et al. (2000, ApJ,
; 533, 682), and for a given E(B-V) value returns the flux attenuation
; array, which can be used to get reddened templates. Here the spectra
; are assumed to be binned on a ln-rebinned wavelentgh grid as defined
; by input l0_gal,lstep_gal,npix parameters. The input receiding
; velocity vstar, is used to derive the dust reddening in the galaxy
; rest-frame.
; 
; Can be used also to de-reddened the object spectra by the Milky-Way
; dust extinction, using as E(B-V) the opposite of the Schlegel et
; al. values found in NED and vstar = 0.
;
; Initial version kindly provided by S. Kaviray, Oxford, 2006.

; reconstruct the wavelength array in Anstroms, and compute rest-frame
; values
lambda = exp(dindgen(npix)*lstep_gal + l0_gal)
if keyword_set(log10) then lambda = 10^(dindgen(npix)*lstep_gal + l0_gal)
lambda = lambda/exp(vstar/299792.458d)

; array to hold k(lambda) values
k = fltarr(n_elements(lambda))           

for i=0,n_elements(lambda)-1 do begin
     ; convert wavelength units from angstroms to micrometres
     l = lambda(i)/1e4                   
     ; assign k values
     if (l ge 0.63 and l le 2.2) then k(i) = 2.659*(-1.857+1.040/l)+4.05
     if (l lt 0.63)              then k(i) = 2.659*(-2.156+1.509/l-0.198/l^2+0.011/l^3)+4.05
     if (l gt 2.2)               then k(i) = 0.0
endfor

return,(10^(-0.4*ebv*k))
; this should be then multiplied by the spectrum flux array
END



;------------------------------------------------------------------------------------
PRO SET_CONSTRAINTS, GALAXY=galaxy, NOISE=noise, CSTAR=cstar,                                     $
                     KINSTARS=kinstars, VELSCALE=velscale, DEGREE=degree, MDEGREE=mdegree,        $
                     GOODPIXELS=goodpixels, EMISSION_SETUP=emission_setup, START_PARS=start_pars, $
                     L0_GAL=l0_gal, LSTEP_GAL=lstep_gal, PARINFO=parinfo, FUNCTARGS=functargs,    $
                     INT_DISP=int_disp, LOG10=log10, REDDENING=reddening, L0_TEMPL=l0_templ,      $
                     FOR_ERRORS=for_errors, vmaxshift=vmaxshift, sigmamax=sigmamax

    npix = n_elements(galaxy)
    
; This subroutine sets up the constraints and boundaries for the
; variables to be fitted, preparing and returning the PARINFO and
; FUNCTARG structure for MPFIT

; Total number of parameters that MPFIT will deal with, i.e. the
; number of emission lines (counting multiplets only once) times two
; (or three if we are fitting also the amplitudes to get the errors on
; them), plus the order of the multiplicative polynomials and if
; needed the reddening parameters
; 
i_lines   = where(emission_setup.kind eq 'l')
nlines   = n_elements(i_lines)
if (keyword_set(for_errors) eq 0) then n_pars = nlines*2 + mdegree else n_pars = nlines*3 + mdegree
if (n_elements(reddening) ne 0) then n_pars = n_pars + n_elements(reddening)

; Setup the PARINFO structure that will allow us to control the limits
; of our parameters with MPFIT, decide whether we want to hold them at
; the input values, or tie some of them to others.

; jm09dec10ucsd - don't set STEP for the amplitude if /FOR_ERRORS
if keyword_set(for_errors) then begin
   parinfo = REPLICATE({step:velscale/100d,limits:[0d,0d],limited:[1,1],fixed:0, tied:' '}, n_pars)
   parinfo[lindgen(nlines)*3+0].step = 0.0
endif else begin
   parinfo = REPLICATE({step:velscale/100d,limits:[0d,0d],limited:[1,1],fixed:0, tied:' '}, n_pars)
endelse
    

; A) First of all, fix V_gas and S_gas to their input values for
; the lines with a 'h' (for hold) as their fit-kind tag
for i = 0,nlines-1 do begin
    j = i_lines[i]
    if (emission_setup[j].fit eq 'h') then begin
        if (keyword_set(for_errors) eq 0) then parinfo[2*i:2*i+1].fixed = 1 else parinfo[3*i+1:3*i+2].fixed = 1
    endif
endfor

if (n_elements(vmaxshift) eq 0) then vmaxshift = 500d
if (n_elements(sigmamax) eq 0) then sigmamax = 1000d

; B) Second, set the limits
; i) for V_gas and S_gas
if (keyword_set(for_errors) eq 0) then begin
    for i=0,nlines*2 -1, 2 do begin
; jm09nov18ucsd - make VMAXSHIFT and SIGMAMAX free parameters; also
; restrict the center to be >1 and <npix-2 to avoid edge effects
        parinfo[i].limits   = start_pars[i] + [-vmaxshift,vmaxshift]/velscale ; Limits for Vgas (from km/s to pixels)
;       parinfo[i].limits[0] = parinfo[i].limits[0]>1.0
;       parinfo[i].limits[1] = parinfo[i].limits[1]<(npix-2.0)
        parinfo[i+1].limits = [-sigmamax,sigmamax]/velscale ; Limits for Sgas (from km/s to pixels)
;       parinfo[i].limits   = start_pars[i] + [-1d3,1d3]/velscale ; Limits for Vgas (from km/s to pixels)
;       parinfo[i+1].limits = [-1d4,1d4]/velscale                ; Limits for Sgas (from km/s to pixels)
        ; to avoid problems near the previous boundaries
        if (start_pars[i]   le parinfo[i].limits[0]  ) then start_pars[i]   = parinfo[i].limits[0]  +0.0001
        if (start_pars[i]   ge parinfo[i].limits[1]  ) then start_pars[i]   = parinfo[i].limits[1]  -0.0001
        if (start_pars[i+1] le parinfo[i+1].limits[0]) then start_pars[i+1] = parinfo[i+1].limits[0]+0.0001
        if (start_pars[i+1] ge parinfo[i+1].limits[1]) then start_pars[i+1] = parinfo[i+1].limits[1]-0.0001
     endfor
endif else begin
    for i=0,nlines*3 -1, 3 do begin
        parinfo[i+1].limits = start_pars[i+1] + [-vmaxshift,vmaxshift]/velscale ; Limits for Vgas (from km/s to pixels)
;       parinfo[i+1].limits[0] = parinfo[i].limits[0]>1.0
;       parinfo[i+1].limits[1] = parinfo[i].limits[1]<(npix-2.0)
        parinfo[i+2].limits = [-sigmamax,sigmamax]/velscale                ; Limits for Sgas (from km/s to pixels)
;       parinfo[i+1].limits = start_pars[i+1] + [-1d3,1d3]/velscale ; Limits for Vgas (from km/s to pixels)
;       parinfo[i+2].limits = [-1d4,1d4]/velscale                  ; Limits for Sgas (from km/s to pixels)
        ; to avoid problems near the previous boundaries
        if (start_pars[i+1] le parinfo[i+1].limits[0]) then start_pars[i+1] = parinfo[i+1].limits[0]+0.0001
        if (start_pars[i+1] ge parinfo[i+1].limits[1]) then start_pars[i+1] = parinfo[i+1].limits[1]-0.0001
        if (start_pars[i+2] le parinfo[i+2].limits[0]) then start_pars[i+2] = parinfo[i+2].limits[0]+0.0001
        if (start_pars[i+2] ge parinfo[i+2].limits[1]) then start_pars[i+2] = parinfo[i+2].limits[1]-0.0001
    endfor
endelse

; ii) for the mult. polynomial (if needed) 
if mdegree ge 1 then begin
    parinfo[n_pars-mdegree:*].limits = [-1d,1d]
    parinfo[n_pars-mdegree:*].step = 1d-3
endif

; iii) and for the reddening parameters (if needed). These will follow the
; emission-line parameters in the parameter array. 
if (n_elements(reddening) ne 0) then begin
    if (keyword_set(for_errors) eq 0) then begin
        parinfo[nlines*2:nlines*2+n_elements(reddening)-1].limits = [0d,5d]
        parinfo[nlines*2:nlines*2+n_elements(reddening)-1].step = 1d-3
    endif else begin
        parinfo[nlines*3:nlines*3+n_elements(reddening)-1].limits = [0d,5d]
        parinfo[nlines*3:nlines*3+n_elements(reddening)-1].step = 1d-3
    endelse
endif

; iv) Also, if we are also fitting the amplitudes non-linearly
; with MPFIT to derive errors, constrain the latter to be
; non-negative, unless we have specified negative amplitudes in the
; emission-line setup. In this case we want non-positive values.
if keyword_set(for_errors) then begin
    for i = 0,nlines-1 do begin
        j = i_lines[i]
        if (emission_setup[j].a gt 0) then begin
            parinfo[i*3].limited = [1,0]
            parinfo[i*3].limits[0] = 0
        endif else begin
            parinfo[i*3].limited = [0,1]
            parinfo[i*3].limits[1] = 0
        endelse
    endfor
endif

; C) Finally, find the lines for which the kinematics needs to be tied
; to that of other lines. These have either 't#', 'v#' or 's#'
; fit-kind tags, where # indicates the index of the line to which they
; are tied. Notice that the # indexes are the ones listed in the .i
; tag of the emission-lines setup structure
;
; If 't#' we tie both the position and width of the lines, if 'v#'
; only the position, and if 's#' only the width.
for i = 0,nlines-1 do begin
    j = i_lines[i]    
    ; check if we have a 't' tag
    if (strmid(emission_setup[j].fit,0,1) eq 't') then begin
        ; find the reference line reference index, as given in the
        ; 'fit' tag of the emission setup structure
        k_refline = fix(strmid(emission_setup[j].fit,1)) 
        ; which correspond to the following position in emission setup
        ; structure that was passed to this function
        j_refline = where(emission_setup.i eq k_refline)
        if (j_refline eq -1) then $
          message,'Hey, you tied '+emission_setup[j].name+' to a line you are not even fitting...'
        ; and to the following position in the list of the emission
        ; lines to fit
        i_refline = where(emission_setup[i_lines].i eq k_refline)
        ; 
        l_line    = emission_setup[j].lambda         & str_l_line    = string(l_line)
        l_refline = emission_setup[j_refline].lambda & str_l_refline = string(l_refline)
        ;
        if (keyword_set(for_errors) eq 0) then begin
            parinfo[2*i].tied = $ 
              strcompress('P['+string(2*i_refline)+']-alog('+str_l_refline+'/'+str_l_line+')/')+ $
              strtrim(string(lstep_gal),2)
            ; to deal with log10-lambda rebinned data, instead of ln-lambda
            if keyword_set(log10) then parinfo[2*i].tied = $
              strcompress('P['+string(2*i_refline)+']-alog10('+str_l_refline+'/'+str_l_line+')/')+ $
              strtrim(string(lstep_gal),2)
            parinfo[2*i+1].tied =  strcompress('P['+string(2*i_refline+1)+']')        
        endif else begin
            parinfo[3*i+1].tied = $ 
              strcompress('P['+string(3*i_refline+1)+']-alog('+str_l_refline+'/'+str_l_line+')/')+ $
              strtrim(string(lstep_gal),2)
            ; to deal with log10-lambda rebinned data, instead of ln-lambda
            if keyword_set(log10) then parinfo[3*i+1].tied = $
              strcompress('P['+string(3*i_refline+1)+']-alog10('+str_l_refline+'/'+str_l_line+')/')+ $
              strtrim(string(lstep_gal),2)
            parinfo[3*i+2].tied =  strcompress('P['+string(3*i_refline+2)+']')        
        endelse        
    endif
    ; check if we have a 'v' tag
    if (strmid(emission_setup[j].fit,0,1) eq 'v') then begin
        ; find the reference line reference index, as given in the
        ; 'fit' tag of the emission setup structure
        k_refline = fix(strmid(emission_setup[j].fit,1)) 
        ; which correspond to the following position in emission setup
        ; structure that was passed to this function
        j_refline = where(emission_setup.i eq k_refline)
        if (j_refline eq -1) then $
          message,'Hey, you tied '+emission_setup[j].name+' to a line you are not even fitting...'
        ; and to the following position in the list of the emission
        ; lines to fit
        i_refline = where(emission_setup[i_lines].i eq k_refline)
        ; 
        l_line    = emission_setup[j].lambda         & str_l_line    = string(l_line)
        l_refline = emission_setup[j_refline].lambda & str_l_refline = string(l_refline)
        ;
        if not keyword_set(for_errors) then begin
            parinfo[2*i].tied = $
              strcompress('P['+string(2*i_refline)+']-alog('+str_l_refline+'/'+str_l_line+')/')+ $
              strtrim(string(lstep_gal),2)
            ; to deal with log10-lambda rebinned data, instead of ln-lambda
            if keyword_set(log10) then parinfo[2*i].tied = $
              strcompress('P['+string(2*i_refline)+']-alog10('+str_l_refline+'/'+str_l_line+')/')+ $
              strtrim(string(lstep_gal),2)
        endif else begin
            parinfo[3*i+1].tied = $
              strcompress('P['+string(3*i_refline+1)+']-alog('+str_l_refline+'/'+str_l_line+')/')+ $
              strtrim(string(lstep_gal),2)
            ; to deal with log10-lambda rebinned data, instead of ln-lambda
            if keyword_set(log10) then parinfo[3*i+1].tied = $
              strcompress('P['+string(3*i_refline+1)+']-alog10('+str_l_refline+'/'+str_l_line+')/')+ $
              strtrim(string(lstep_gal),2)
        endelse
    endif
    ; check if we have a 's' tag
    if (strmid(emission_setup[j].fit,0,1) eq 's') then begin
        ; find the reference line reference index, as given in the
        ; 'fit' tag of the emission setup structure
        k_refline = fix(strmid(emission_setup[j].fit,1)) 
        ; which correspond to the following position in emission setup
        ; structure that was passed to this function
        j_refline = where(emission_setup.i eq k_refline)
        if (j_refline eq -1) then $
          message,'Hey, you tied '+emission_setup[j].name+' to a line you are not even fitting...'
        ; and to the following position in the list of the emission
        ; lines to fit
        i_refline = where(emission_setup[i_lines].i eq k_refline)
        ; 
        l_line    = emission_setup[j].lambda         & str_l_line    = string(l_line)
        l_refline = emission_setup[j_refline].lambda & str_l_refline = string(l_refline)
        if not keyword_set(for_errors) then begin
            parinfo[2*i+1].tied =  strcompress('P['+string(2*i_refline+1)+']')        
        endif else begin
            parinfo[3*i+2].tied =  strcompress('P['+string(3*i_refline+2)+']')        
        endelse
    endif
endfor

; Gathering all the information and setting the FUNCTARGS structure for MPFIT
if not keyword_set(for_errors) then for_errors =0 else for_errors =1
if not keyword_set(log10) then begin
    if (n_elements(reddening) eq 0) then begin
        functargs={CSTAR:cstar, GALAXY:galaxy, NOISE:noise, EMISSION_SETUP:emission_setup, $
                   KINSTARS:kinstars, VELSCALE:velscale, DEGREE:degree, MDEGREE:mdegree, $
                   GOODPIXELS:goodpixels, L0_GAL:l0_gal, LSTEP_GAL:lstep_gal, $
                   INT_DISP:int_disp, $
                   FOR_ERRORS:for_errors}
    endif else begin
        functargs={CSTAR:cstar, GALAXY:galaxy, NOISE:noise, EMISSION_SETUP:emission_setup,$
                   KINSTARS:kinstars, VELSCALE:velscale, DEGREE:degree, MDEGREE:mdegree, $
                   GOODPIXELS:goodpixels, L0_GAL:l0_gal, LSTEP_GAL:lstep_gal, $
                   INT_DISP:int_disp, REDDENING:reddening, L0_TEMPL:l0_templ, $
                   FOR_ERRORS:for_errors}
    endelse
endif else begin
    if (n_elements(reddening) eq 0) then begin
        functargs={CSTAR:cstar, GALAXY:galaxy, NOISE:noise, EMISSION_SETUP:emission_setup, $
                   KINSTARS:kinstars, VELSCALE:velscale, DEGREE:degree, MDEGREE:mdegree, $
                   GOODPIXELS:goodpixels, L0_GAL:l0_gal, LSTEP_GAL:lstep_gal, $
                   INT_DISP:int_disp, LOG10:log10, $
                   FOR_ERRORS:for_errors}
    endif else begin
        functargs={CSTAR:cstar, GALAXY:galaxy, NOISE:noise, EMISSION_SETUP:emission_setup,$
                   KINSTARS:kinstars, VELSCALE:velscale, DEGREE:degree, MDEGREE:mdegree, $
                   GOODPIXELS:goodpixels, L0_GAL:l0_gal, LSTEP_GAL:lstep_gal, $
                   INT_DISP:int_disp, LOG10:log10, REDDENING:reddening, L0_TEMPL:l0_templ, $
                   FOR_ERRORS:for_errors}
    endelse
endelse

END

;------------------------------------------------------------------------------------
FUNCTION FITFUNC_GAS, pars, CSTAR=cstar, GALAXY=galaxy, NOISE=noise, EMISSION_SETUP=emission_setup, $
                      KINSTARS=kinstars, VELSCALE=velscale, DEGREE=degree, MDEGREE=mdegree,         $
                      GOODPIXELS=goodpixels, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal,                    $
                      BESTFIT=bestfit, WEIGHTS=weights, EMISSION_TEMPLATES=emission_templates,      $
                      INT_DISP=int_disp, LOG10=log10, REDDENING=reddening, L0_TEMPL=l0_templ,       $
                      FOR_ERRORS=for_errors
compile_opt idl2, hidden


npix = n_elements(galaxy)
x = range(-1d,1d,npix)             ; X needs to be within [-1,1] for Legendre Polynomials

nlines = n_elements(where(emission_setup.kind eq 'l'))
IF NOT KEYWORD_SET(for_errors) THEN npars = nlines*2 ELSE npars = nlines*3
;print,npars
; append the reddening parameters if needed
IF (n_elements(reddening) ne 0) THEN npars = npars + n_elements(reddening)

; The zero order multiplicative term is already included in the
; linear fit of the templates. The polinomial below has mean of 1.
mpoly= 1d       ; The loop below can be null if mdegree < 1
FOR j=1,mdegree DO mpoly=mpoly + legendre(x,j)*pars[npars+j-1]

; Emission Lines as given by the values in pars
s = size(cstar)
; passing only the emission-line parameters 
IF (KEYWORD_SET(for_errors) eq 0) THEN eml_pars = pars[0:nlines*2-1] ELSE eml_pars = pars[0:nlines*3-1]
int_disp_pix = int_disp/velscale
gaus = gandalf_create_templates(EMISSION_SETUP=emission_setup,LSTEP_GAL=lstep_gal, NPIX=npix, PARS=eml_pars, $
                        INT_DISP_PIX=int_disp_pix, LOG10=log10, $
                        FOR_ERRORS=for_errors)
; Stacking all the inputs together:
;   1.- Legendre polinomials of order 'degree'
;   2.- Convolved SSP models (pre-convolved by the best LOSVD in set_constraints)
;       adjusted by a multiplicative polinomials - or by reddening
;   3.- Emission Lines, also reddened - NO! jm09nov16ucsd
IF s[0] EQ 2 THEN ntemp = s[2] ELSE ntemp = 1   ; Number of template spectra
c = dblarr(npix,degree+nlines+ntemp+1,/NOZERO)  ; This array is used for estimating predictions
a = c               	    	    	    	; This array is used for the actual solution of the system
FOR j=0,degree   DO c[*,j] = legendre(x,j)      ; Legendre polinomials (additive)
IF (n_elements(reddening) eq 0) THEN BEGIN
    FOR j=0,ntemp-1  DO c[*,degree+1+j] = mpoly*cstar[0:npix-1,j] ; Convolved templates x mult. polinomials
    FOR j=0,nlines-1 DO c[*,degree+ntemp+1+j]=gaus[*,j]           ; Emission lines
ENDIF ELSE BEGIN
    ; redden both stellar and emission-line templates - NO! jm09nov16ucsd
    IF (KEYWORD_SET(for_errors) eq 0) THEN ebv = pars[nlines*2] ELSE ebv = pars[nlines*3]
    Vstar =  kinstars[0] + (l0_gal-l0_templ)*299792.458d
    if KEYWORD_SET(log10) then  Vstar =  kinstars[0] + (l0_gal-l0_templ)*299792.458d*alog(10.0d)
    reddening_attenuation = DUST_CALZETTI(l0_gal,lstep_gal,npix,ebv,Vstar,LOG10=log10)
    ; but also include extra internal reddening if requested 
    IF n_elements(reddening) EQ 2 THEN BEGIN
       message, 'Not supported!' ; jm09nov16ucsd
        IF (KEYWORD_SET(for_errors) eq 0) THEN int_ebv = pars[nlines*2+1] ELSE int_ebv = pars[nlines*3+1]
        int_reddening_attenuation = DUST_CALZETTI(l0_gal,lstep_gal,npix,int_ebv,Vstar,LOG10=log10)
    ENDIF ELSE int_reddening_attenuation = 1.0d
    FOR j=0,ntemp-1  DO c[*,degree+1+j] = cstar[0:npix-1,j]*reddening_attenuation
    FOR j=0,nlines-1 DO c[*,degree+ntemp+1+j]=gaus[*,j] ; jm09nov16ucsd
;   FOR j=0,nlines-1 DO c[*,degree+ntemp+1+j]=gaus[*,j]*reddening_attenuation*int_reddening_attenuation 
ENDELSE

FOR j=0,degree+ntemp+nlines DO a[*,j] = c[*,j]/noise          ; Weight all columns with errors

; Select the spectral region to fit and solve the overconditioned system
; using SVD/BVLS. Use unweighted array for estimating bestfit predictions.
sol = BVLSN_Solve_pxf(a[goodpixels,*],galaxy[goodpixels]/noise[goodpixels],degree, $
                      FOR_ERRORS=for_errors,NLINES=nlines)
bestfit = c # sol
err = (galaxy[goodpixels]-bestfit[goodpixels])/noise[goodpixels]

if (total(sol gt 1D20) ge 1) then splog, 'Possible problem!' ; message, 'Problem here' ; test

; output weights for the templates
weights = sol[degree+1:n_elements(sol)-1]
; Make the array containing each of the best matching emission-line templates
;
; Array with the Gaussian templates weigths. 
; In case we have used the keyword FOR_ERRORS, these should all be 1
; and gaus should already have the right amplitudes
sol_gas = sol[degree+ntemp+1:degree+ntemp+nlines]
emission_templates = gaus
IF (n_elements(reddening) eq 0) THEN BEGIN
    for i = 0,n_elements(sol_gas)-1 do emission_templates[*,i] = gaus[*,i]*sol_gas[i]
ENDIF ELSE BEGIN
    ; Provide the emission-line templates as observed, i.e. reddened
    for i = 0,n_elements(sol_gas)-1 do emission_templates[*,i] = gaus[*,i]*sol_gas[i] ; jm09nov16ucsd
;   for i = 0,n_elements(sol_gas)-1 do emission_templates[*,i] = gaus[*,i]*sol_gas[i]* $
;     reddening_attenuation*int_reddening_attenuation
ENDELSE

return, err
END

;------------------------------------------------------------------------------------
PRO REARRANGE_RESULTS, res, weights, chi2, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal,    $
                       VELSCALE=velscale, EMISSION_SETUP=emission_setup, SOL=sol, $
                       INT_DISP=int_disp, LOG10=log10, REDDENING=reddening,       $
                       ERR=in_err, ESOL=esol, FOR_ERRORS=for_errors
err = in_err
; Given the input res array from the MPFIT fit (with V_gas and S_gas
; best-fitting values - in pixels) and the weight from the BVLS fit
; (which with the emission-line basic amplitudes to get A_gas),
; construct the sol solution array, containing F_gas, A_gas, V_gas and
; S_gas (the latter in km/s). Also rearrange the MPFIT formal
; uncertainties in V_gas and S_gas, which should be consider only as
; lower estimates.

; If this routine is called after a second MPFIT fit for A_gas, V_gas
; and S_gas then we also rearrange the corresponding MPFIT formal
; uncertainties, which in this case will be the correct ones.

; NOTE: at the moment this routine does not rearrange mult. polynomial
; coefficient and corresponding uncertainties. This could be easily
; implemented by adding a MDEGREE=mdegree keyword and then append to
; the sol_final and esol_final array the corresponding numbers from
; the res and err input arrays

c = 299792.458d ; Speed of light in km/s
i_lines = where(emission_setup.kind eq 'l')
nlines = n_elements(where(emission_setup.kind eq 'l'))
lambda0 = emission_setup[i_lines].lambda
offset  = (alog(lambda0)-l0_gal)/lstep_gal
; to deal with log10-lambda rebinned data, instead of ln-lambda
if keyword_set(log10) then offset  = (alog10(lambda0)-l0_gal)/lstep_gal

; make final output solution array
sol_final = dblarr(nlines*4)
; make room for the E(B-V) reddening value(s) in the final solution vector
if (n_elements(reddening) ne 0) then sol_final = dblarr(nlines*4+n_elements(reddening))

k=0 & h=0
for i=0,nlines-1 do begin
    if (keyword_set(for_errors) eq 0) then begin
        ; processing outputs from quickest fits, where only V_gas and S_gas
        ; were solved with MPFIT whereas A_gas was left to BVLS
        ampl_i = emission_setup[i_lines[i]].a
        sol_final[k+1] = ampl_i*weights[n_elements(weights)-nlines+i] ; Final amplitude of the Emission line
        sol_final[k+2] = -(offset[i]-res[h])*velscale                 ; Radial velocity of the Emission line [km/s]
        sol_final[k+3] = abs(res[h+1]*velscale)                       ; Sigma (intrinsic!) of the Emission line [km/s]
    endif else begin
        ; processing outputs from fits meant to properly derive errors,
        ; where also A_gas was solved with MPFIT
        sol_final[k+1] = res[h]                                ; Amplitude of the Emission lines
        sol_final[k+2] = -(offset[i]-res[h+1])*velscale        ; Radial velocity of the Emission lines [km/s]
        sol_final[k+3] = abs(res[h+2]*velscale)                ; Sigma (intrinsic!) of the Emission lines [km/s]
    endelse
    ; Sigma (as observed!)
    sigma          = sqrt(sol_final[k+3]^2. + int_disp^2.)
    ; Flux of the Emission lines
    sol_final[k]   = sol_final[k+1]* sqrt(2*!pi) * sigma * lambda0[i] * exp(sol_final[k+2]/c)/c 
    k=k+4
    if (keyword_set(for_errors) eq 0) then h=h+2 else h=h+3
endfor
; Append reddening values to the final solution vector, after the
; emission-line parameters
if (keyword_set(for_errors) eq 0) then begin
   if (n_elements(reddening) ne 0) then sol_final[nlines*4:nlines*4+n_elements(reddening)-1] = $
      res[nlines*2:nlines*2+n_elements(reddening)-1]
endif else begin
    if (n_elements(reddening) ne 0) then sol_final[nlines*4:nlines*4+n_elements(reddening)-1] = $
      res[nlines*3:nlines*3+n_elements(reddening)-1]
endelse

;; --------------------------------------------------
;; jm09dec08ucsd - lines/parameters that are tied have zero error
;; returned by MPFIT; fix that here
;if (n_elements(err) ne 0) then begin
;   for i = 0,nlines-1 do begin
;      if (strmid(emission_setup[i_lines[i]].fit,0,1) eq 't') then begin
;         i_refline = where(emission_setup[i_lines].i eq $
;           fix(strmid(emission_setup[i_lines[[i]]].fit,1)))
;         if (keyword_set(for_errors) eq 0) then begin
;            err[2*i+0] = err1[2*i_refline+0] ; redshift
;            err[2*i+1] = err1[2*i_refline+1] ; sigma-width
;         endif else begin
;;           factor = res[3*i_refline+2]/(res[3*i+2]+(res[3*i+2] eq 0))*(res[3*i+2] ne 0)
;            factor = res[3*i+0]/(res[3*i_refline+0]+(res[3*i_refline+0] eq 0))*(res[3*i_refline+0] ne 0)
;            err[3*i+0] = err1[3*i_refline+0]*factor ; scaled amplitude
;            err[3*i+1] = err1[3*i_refline+1] ; redshift
;            err[3*i+2] = err1[3*i_refline+2] ; sigma-width
;         endelse
;      endif
;   endfor
;endif
;; --------------------------------------------------
   
if (n_elements(err) ne 0) then begin
    ; make final output error array
    esol_final = dblarr(nlines*4)
    ; make room for errors in the reddening paramenter(s)
    if (n_elements(reddening) ne 0) then esol_final = dblarr(nlines*4+n_elements(reddening))
    k=0 & h=0

    if (keyword_set(for_errors) eq 0) then begin
        ; MPFIT errors only from V_gas and S_gas fit 
        for i=0,nlines-1 do begin
            esol_final[k]   = 0d
            esol_final[k+1] = 0d
            esol_final[k+2] = err[h]  *velscale ; these are almost certain lower limits for
            esol_final[k+3] = err[h+1]*velscale ; the real uncertainties on these parameters
            k=k+4 & h=h+2
        endfor
    endif else begin
        ; MPFIT errors only from for non-linear fit of A_gas, V_gas and S_gas fit 
       nmc = 100 ; number of Monte Carlo iterations
       err1 = err
       for i=0,nlines-1 do begin
;; jm09nov17ucsd - mpfit does not return errors on lines that have been
;; tied to other lines; therefore the velocity and velocity dispersion
;; error on most lines is zero; put the right values here!
            if (strmid(emission_setup[i_lines[i]].fit,0,1) eq 't') then begin ; is this line tied?
               refline = where(emission_setup[i_lines].i eq fix(strmid(emission_setup[i_lines[i]].fit,1)))
; scale the amplitude error by the ratio of the amplitudes
               factor = res[h+0]/(res[3*refline+0]+(res[3*refline+0] eq 0))*(res[3*refline+0] ne 0)
               if (factor eq 0.0) then err1[h+0] = err[h+0] else $ ; amplitude
                 err1[h+0] = factor*err[refline*3+0]
               err1[h+1] = err[refline*3+1] ; redshift
               err1[h+2] = err[refline*3+2] ; sigma-width
            endif
            esol_final[k+1] = err1[h+0]          ; amplitude
            esol_final[k+2] = err1[h+1]*velscale ; velocity
            esol_final[k+3] = err1[h+2]*velscale ; dispersion
            ; Simple MC error propagation to derive errors in F_gas
            ; This implicitly assumes errors in A_gas ,V_gas and S_gas
            ; are uncorrelated...

; jm09dec10ucsd - the original code requires a non-zero velocity and
; sigma-width *error*; however, if either of these parameters hit
; their boundaries (e.g., sigma-width=300 km/s) then MPFIT returns
; null errors, which means that we also don't get an estimate of the
; *flux error*, even if the line is indeed significant; therefore,
; loosen the restrictions here, recognizing that the flux error will
; be underestimated for these lines 

; jm09dec10ucsd - another occasional problem is that the parameters
; can be reasonable, but the errors are ginormous (like >10^5);
; so...if the S/N on the amplitude is <10^-3 then set the amplitude,
; amplitude error, flux, and flux error equal to zero
;           if esol_final[k+2] gt 0 and esol_final[k+3] gt 0 then begin
            if (sol_final[k+1]/(esol_final[k+1]+(esol_final[k+1] eq 0))*(esol_final[k+1] ne 0) gt 1E-3) and $
              (esol_final[k+1] gt 0) then begin ; AMP_ERR>0 *and* AMP/AMP_ERR>1E-3
               ampls_i   = sol_final[k+1] + esol_final[k+1]*randomn(seed,nmc)
               vels_i    = sol_final[k+2] + esol_final[k+2]*randomn(seed,nmc)
               sigmas_i  = sol_final[k+3] + esol_final[k+3]*randomn(seed,nmc)
               sigmas_i  = sqrt(sigmas_i^2. + int_disp^2.)
               
               fluxes_i = ampls_i* sqrt(2*!pi) * sigmas_i * lambda0[i] * exp(vels_i/c)/c 
               esol_final[k+0]  = im_robust_sigma(fluxes_i)
               if (esol_final[k+0] eq 0) then message, 'Error is zero!'
            endif else begin
               sol_final[k+0] = 0.0
               sol_final[k+1] = 0.0
               esol_final[k+0] = 0.0
               esol_final[k+1] = 0.0
            endelse
            if (sol_final[k+1] gt 0) and (esol_final[k+1] eq 0) then stop
            k=k+4 & h=h+3
         endfor
    endelse
    
    ; Add reddening errors
    if (keyword_set(for_errors) eq 0) then begin
        if (n_elements(reddening) ne 0) then esol_final[nlines*4:nlines*4+n_elements(reddening)-1] = $
          err[nlines*2:nlines*2+n_elements(reddening)-1]
    endif else begin
        if (n_elements(reddening) ne 0) then esol_final[nlines*4:nlines*4+n_elements(reddening)-1] = $
          err[nlines*3:nlines*3+n_elements(reddening)-1]
    endelse

endif

sol=sol_final
if (n_elements(err) ne 0) then esol=esol_final

END

;------------------------------------------------------------------------------------
PRO SHOW_FIT,galaxy, bestfit, emission, best_pars, sol, GOODPIXELS=goodpixels, MDEGREE=mdegree, $
             REDDENING=reddening, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal, L0_TEMPL=l0_templ, LOG10=log10, $
             KINSTARS=kinstars, NLINES=nlines
; in progress...

; plot final data-model comparison if required.
; the colors below were chosen for the sauron colormap.
    s2 = size(galaxy)
    sauron_colormap & device,decompose=0
    mn = min(galaxy, max=mx) & mn = mn -0.1*(mx-mn)
    resid = mn + galaxy - bestfit
    y1 = min(resid[goodpixels]) & y2 = mx & diff = y2-y1
    y1 = y1 - 0.15*diff & y2 = y2 + 0.1*diff
    plot, galaxy, xtitle='pixels', ytitle='counts', xstyle=1, /ynozero, psym=10, $
      yrange=[y1,y2], ystyle=1,xrange=[-0.02,1.02]*s2[1]
    oplot, bestfit, color=210, thick=2
    oplot, bestfit-emission, color=210, linestyle=2
    oplot, emission + mn, color=80
    oplot, goodpixels, resid[goodpixels], psym=3, color=128, symsize=0.3
    n = n_elements(goodpixels)
    oplot, goodpixels, replicate(mn,n), psym=3, color=128, symsize=2
    w = where(goodpixels - shift(goodpixels,1) gt 1, m)
    if m gt 0 then w = [0,w-1,w,n-1] else w = [0,n-1] ; add first and last point
    for j=0,n_elements(w)-1 do oplot, goodpixels[w[[j,j]]], [mn,galaxy[goodpixels[w[j]]]], color=128
    resid_noise = im_robust_sigma(resid[goodpixels]-mn, /zero)
    oplot, emission*0 + resid_noise +mn,linestyle=2

    ; plots the polinomial correction and the unadjusted templates
    x = range(-1d,1d,n_elements(galaxy)) ; x needs to be within [-1,1] for legendre polynomials
    mpoly= 1d                            ; the loop below can be null if mdegree < 1
    npars=n_elements(best_pars) - mdegree
    if (n_elements(reddening) eq 0) then begin
        for j=1,mdegree do mpoly = mpoly + legendre(x,j)*best_pars[npars+j-1]
        oplot,(bestfit-emission)/mpoly,col=140,linestyle=1
    endif else begin
        c = 299792.458d         ; speed of light in km/s
        vstar =  kinstars[0] + (l0_gal-l0_templ)*c
        if keyword_set(log10) then  vstar =  kinstars[0] + (l0_gal-l0_templ)*c*alog(10.0d)
        reddening_attenuation = $
          dust_calzetti(l0_gal,lstep_gal,n_elements(galaxy),sol[nlines*4],vstar,log10=log10)        
        if n_elements(reddening) eq 2 then begin
            int_reddening_attenuation = $
              dust_calzetti(l0_gal,lstep_gal,n_elements(galaxy),sol[nlines*4+1],vstar,log10=log10)
        endif else int_reddening_attenuation = 1.0 
        oplot,bestfit/reddening_attenuation-emission,col=140,linestyle=1 ; jm09nov16ucsd
;       oplot,((bestfit-emission)/reddening_attenuation),col=140,linestyle=1
    endelse
    legend,['data','fit','res.-noise'], $
      color=[255,210,255],linestyle=[0,0,2],box=0,/bottom,/left
    legend,['residuals','emission-lines','unadjusted continuum'], $
      color=[128,80,140],linestyle=[1,0,1],box=0,/bottom,/right

END


;------------------------------------------------------------------------------------
PRO im_GANDALF, templates, galaxy, noise, velScale, sol, emission_setup, l0_gal, lstep_gal,       $
  GOODPIXELS=goodPixels, DEGREE=degree, MDEGREE=mdegree, INT_DISP=int_disp,            $
  BESTFIT=bestFit, EMISSION_TEMPLATES=emission_templates, WEIGHTS=weights, ERROR=esol, $
  PLOT=plot, QUIET=quiet, LOG10=log10, REDDENING=reddening, L0_TEMPL=l0_templ,         $
  FOR_ERRORS=for_errors, vmaxshift_iter1=vmaxshift_iter1, vmaxshift_iter2=vmaxshift_iter2, $
  sigmamax_iter1=sigmamax_iter1, sigmamax_iter2=sigmamax_iter2, chi2=chi2
compile_opt idl2
;on_error, 2

; ------------------------------------
; Do some initial input error checking
s1 = size(templates)
s2 = size(galaxy)
s3 = size(noise)
IF (s1[0] GT 2 OR s2[0] NE 1 OR s3[0] NE 1) THEN message, 'Wrong input dimensions'
IF (s2[1] NE s3[1]) THEN message, 'GALAXY and NOISE must have the same size'
IF (s1[1] LT s2[1]) THEN message, 'STAR length cannot be smaller than GALAXY'
IF n_elements(degree)  LE 0 THEN degree = -1 ELSE degree = degree > (-1)
IF n_elements(mdegree) LE 0 THEN mdegree = 0 ELSE mdegree = mdegree > 0
IF n_elements(goodPixels) LE 0 THEN goodPixels = indgen(s2[1]) ELSE goodPixels = goodPixels
IF max(goodPixels) GT s2[1]-1 THEN message, 'GOODPIXELS are outside the data range'
IF n_elements(int_disp)  LE 0 THEN int_disp = 0d ELSE int_disp = int_disp
;  do not allow use simultaneous use of reddening and polynomials.
IF n_elements(reddening) GT 0 AND degree NE -1 THEN message, 'Reddening & polynomial adjust. cannot be used together'
IF n_elements(reddening) GT 0 AND mdegree NE 0 THEN message, 'Reddening & polynomial adjust. cannot be used together'
;F n_elements(reddening) ne 1 THEN message, 'Reddening must be a scalar' ; jm09nov16ucsd
;IF n_elements(reddening) GT 2 THEN message, 'Sorry, can only deal with two dust components...'

; ------------------------------------
; First of all find the emission-lines which we are effectively going
; to fit.  That is, exclude from the input structure the lines that
; are either being masked or not fitted.

; jm09nov19ucsd - also ignore lines that are outside the wavelength
; range because it can cause mpfit problems on the edges when lines
; are not tied to one another

;npad = 4
;vpixel = (alog(emission_setup.lambda)-l0_gal)/lstep_gal ; position in pixels
;for ii = 0, n_elements(emission_setup.name)-1 do begin
;   good = where((emission_setup.action[ii] eq 'f') and $
;     (vpixel[ii] gt npad) and (vpixel[ii] lt (n_elements(galaxy)-npad-1)),$
;     ngood,comp=ignore,ncomp=nignore)
;   if (nignore ne 0) then begin
;      emission_setup.action[ii] = 'i' ; modifying the input structure here!!!!!
;; check to see whether [NII] 6584 or [OIII] 5007 were dropped, in
;; which case we drop the weaker component, too
;      istied = where(strmid(emission_setup.kind,1) eq emission_setup.i[ii],nistied)
;      if (nistied ne 0) then emission_setup.action[istied] = 'i' ; note!
;   endif 
;endfor
;
;; now we have to loop back through and check to see if any of the blue
;; rest-wavelength Balmer or forbidden lines are tied to a line that
;; has just been dropped (e.g., H-gamma and H-alpha tied to H-beta)
;for ii = 0, n_elements(emission_setup.name)-1 do begin
;   istied = where((emission_setup.action[ii] eq 'i') and $
;     strmid(emission_setup.fit,1) eq emission_setup.i[ii],nistied)
;   if (nistied ne 0) then begin
;      this = where((emission_setup.action[istied] eq 'f'),nthis)
;      if (nthis ne 0) then begin
;         this = this[nthis-1]   ; pick the reddest not-ignored line
;         emission_setup.fit[istied] = 't'+strtrim(emission_setup.i[istied[this]],2)
;      endif
;   endif
;endfor

dummy = emission_setup ; This will help us restoring later the input emission_setup structure
i_f = where(emission_setup.action eq 'f')
emission_setup = dummy[i_f]
;emission_setup = create_struct('i',dummy.i[i_f],'name',dummy.name[i_f],$
;                               'lambda',dummy.lambda[i_f],'action',dummy.action[i_f],$
;                               'kind',dummy.kind[i_f],'a',dummy.a[i_f],$
;                               'v',dummy.v[i_f],'s',dummy.s[i_f],$
;                               'fit',dummy.fit[i_f])
; Count the number of single lines or the number of multiplets
i_lines  = where(emission_setup.kind eq 'l')
nlines   = n_elements(i_lines)

; ------------------------------------
; Make sure that the input amplitudes of each single line or of the
; main lines of each multiplsts are either 1 or -1. 
for i=0,nlines-1 do begin
    ampl_i = emission_setup[i_lines[i]].a
    if ampl_i gt 0 then emission_setup[i_lines[i]].a =  1.0
    if ampl_i lt 0 then emission_setup[i_lines[i]].a = -1.0
endfor

; ------------------------------------
; Declare and fill in the array with the starting guesses for the
; parameter to fit with MPFIT, namely V_gas and S_gas for the emission
; lines parameters, the coefficients of the mult. polynomials and, if
; necessary the reddening parameters. The latter are placed
; right after the emission-line parameters.
if (n_elements(reddening) eq 0) then $
  start_pars = dblarr(2*nlines+mdegree) else $
  start_pars = dblarr(2*nlines+n_elements(reddening)+mdegree) 

; Loop over the lines, and assign starting V_gas and S_gas as from the
; emission-line setup structure. These are set as starting position
; and width in pixel units.
h = 0
for i = 0,nlines-1 do begin 
    ; current emission-line index in the input setup structure	
    j = i_lines[i]
    ; to deal with log10-lambda rebinned data, instead of ln-lambda
    if not keyword_set(log10) then $
      offset = (alog(emission_setup[j].lambda)-l0_gal)/lstep_gal else $
      offset = (alog10(emission_setup[j].lambda)-l0_gal)/lstep_gal
    ;
    start_pars[h+0] = emission_setup[j].v/velscale + offset
    start_pars[h+1] = emission_setup[j].s/velscale
    h = h + 2
endfor
; add if necessary the starting E(B-V) values
if (n_elements(reddening) ne 0) then start_pars[2*nlines:2*nlines+n_elements(reddening)-1] = reddening

; ------------------------------------
; Convolve the input stellar templates with the input stellar kinematics
kinstars = sol[0:5]
cstar=convolve_templates(templates,kinstars,velscale)

; #########################
; jm09nov19ucsd - iterate the line-fitting twice here; on the second
; iteration use the previous fit as the starting parameter values, and
; tighten the constrains on VMAXSHIFT and SIGMAMAX

if (n_elements(vmaxshift_iter1) eq 0) then vmaxshift_iter1 = 500.0D
if (n_elements(vmaxshift_iter2) eq 0) then vmaxshift_iter2 = 300.0D
if (n_elements(sigmamax_iter1) eq 0) then sigmamax_iter1 = 500.0D
if (n_elements(sigmamax_iter2) eq 0) then sigmamax_iter2 = 500.0D

for iter = 0, 1 do begin

   if (keyword_set(quiet) eq 0) then splog, 'Iteration ', iter, vmaxshift
   
   if (iter eq 0) then begin
      vmaxshift = vmaxshift_iter1
      sigmamax = sigmamax_iter1
   endif else begin
      vmaxshift = vmaxshift_iter2
      sigmamax = sigmamax_iter2
;     splog, 'Starting parameters:'
;     print, reform(start_pars[0:2*nlines-1],2,nlines)*rebin([1.0,velscale],2,nlines)
;     splog, 'Best-fit parameters:'
;     print, reform(best_pars[0:2*nlines-1],2,nlines)*rebin([1.0,velscale],2,nlines)
      start_pars = best_pars
; on the second iteration, loosen the fitting restrictions
      emission_setup.fit = emission_setup.fit_iter2
   endelse

; ------------------------------------
; Set the limits and the appropriate inter-dependencies for fitting
; emission-line Gaussians and prepare the FUNCTARGS and PARINFO
; arguments to be passed to MPFIT
   set_constraints, GALAXY=galaxy, NOISE=noise, CSTAR=cstar, KINSTARS=kinstars, $
     VELSCALE=velscale, DEGREE=degree, MDEGREE=mdegree, GOODPIXELS=goodpixels, $
     EMISSION_SETUP=emission_setup, START_PARS=start_pars, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal, $
     PARINFO=parinfo, FUNCTARGS=functargs, INT_DISP=int_disp, $
     LOG10=log10, REDDENING=reddening, L0_TEMPL=l0_templ, vmaxshift=vmaxshift, sigmamax=sigmamax

; ------------------------------------
; This is where the GANDALF fit is actually performed. Call MPFIT to
; find the best-fitting position and width of the emission lines while
; using BVLS at each iteration to solve for the relative weights of
; the stellar and emission-line templates. The latter weights
; correspond to the emission-line amplitudes. Solve also for the
; mdegree mult. polynomial coeffients and, if needed, also for the
; best reddening parameters.
;
; Note that we evalutate also the errors on the best parameters, but
; as regards the position and width of the lines these should only be
; considered as lower estimates for the real uncertainties.
;  if iter eq 1 then stop

   best_pars = mpfit('FITFUNC_GAS',start_pars, FUNCTARGS=functargs, PARINFO=parinfo, $
     FTOL=1d-2, NFEV=ncalls, ERRMSG=errmsg, PERROR=errors, STATUS=status, QUIET=1)

   if errmsg ne '' then begin
      print, errmsg
      error = best_pars*0.0
   endif
endfor
; #########################

; ------------------------------------
; Call again FITFUNC_GAS with the best paramenters, not only to
; compute the final fit residuals and hence assess the quality of the
; fit, but also to retrieve:
; a) the best fitting template weights   (WEIGHTS)
; b) the best fitting overall model      (BESTFIT)
; c) the best fitting emission templates (EMISSION_TEMPLATES)

; jm10arp28ucsd - KINSTARS=sol[0:5] should be KINSTARS=KINSTARS!!!
resid = fitfunc_gas(best_pars,CSTAR=cstar, GALAXY=galaxy, NOISE=noise, $
                    KINSTARS=kinstars, VELSCALE=velscale, DEGREE=degree, MDEGREE=mdegree, $
                    GOODPIXELS=goodpixels, BESTFIT=bestfit, WEIGHTS=weights, $
                    EMISSION_SETUP=emission_setup, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal, $
                    EMISSION_TEMPLATES=emission_templates, INT_DISP=int_disp, $
                    LOG10=log10, REDDENING=reddening, L0_TEMPL=l0_templ)
    chi2   = im_robust_sigma(resid, /ZERO)^2 
;
if total(noise) eq n_elements(galaxy) then begin
    ; If you have input as errors on the fluxes an array of constant unity vales
    ; compute Chi^2/DOF and use this instead of bestnorm/dof to rescale the formal uncertainties
    chi2   = im_robust_sigma(resid, /ZERO)^2 
    errors  = errors*sqrt(chi2)
endif
; ------------------------------------
; Add up the best-fitting emission templates to get the emission spectrum
if ((size(emission_templates))[0] eq 1) then emission = emission_templates          
if ((size(emission_templates))[0] eq 2) then emission = total(emission_templates,2)

; ------------------------------------
; Rearrange the final results (both best_pars and weights) in the
; output array SOL, which includes also line fluxes. Fill in also the
; ESOL error array.  

rearrange_results, best_pars, weights, chi2, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal, $
              VELSCALE=velscale, EMISSION_SETUP=emission_setup, SOL=sol, $
              ERR=errors, ESOL=esol, INT_DISP=int_disp, $
              LOG10=log10,REDDENING=reddening
; Appends to the best-fitting gas results also the mdegree polynomial
; coefficients.
if mdegree ne 0 then sol = [sol,best_pars[n_elements(best_pars)-mdegree:n_elements(best_pars)-1]]

; Show the fit if requested
if keyword_set(plot) then show_fit, galaxy, bestfit, emission, best_pars, sol, $
  GOODPIXELS=goodpixels, MDEGREE=mdegree, $
  REDDENING=reddening, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal, L0_TEMPL=l0_templ, LOG10=log10, $
  KINSTARS=kinstars, NLINES=nlines

;ww = where(emission_setup.kind eq 'l',nww)
;niceprint, emission_setup.name[ww], sol[indgen(nww)*4+0], $
;  sol[indgen(nww)*4+1], sol[indgen(nww)*4+2], sol[indgen(nww)*4+3]

;plot, galaxy, xr=[1600,1900], ps=10
;djs_oplot, bestfit, ps=10, color='red'
; ------------------------------------
; Properly compute error estimates on all emission-line parameters, by
; solving non-linearly also for the line amplitudes with MPFIT, not
; only for the line positions and widths, as done previously. BVLS
; will now deal only with the weight of the stellar templates. We will
; start such new fit from the previous solution.
IF KEYWORD_SET(FOR_ERRORS) THEN BEGIN
;   print,'computing errors...'

    ; -----------------
    ; Set up the starting guesses for the new fit, including now the amplitudes.
    if (n_elements(reddening) eq 0) then $
      start_pars = dblarr(3*nlines+mdegree) else $
      start_pars = dblarr(3*nlines+n_elements(reddening)+mdegree) 
    ; Fill this in using the values in the solution array sol which
    ; lists, in the order:
    ; (F_gas, A_gas, V_gas, S_gas) for each line
    ; Reddening E(B-V) value(s)     - if any
    ; Mult. polynomial coefficients - if any 
    h = 0
    for i = 0,nlines-1 do begin 
        ; Best-fitting emission-line parameters for this line in the 
        sol_i = sol[i*4+1:i*4+3]
        ; determine offset in pixels
        j = i_lines[i]
        if not keyword_set(log10) then $
          offset = (alog(emission_setup[j].lambda)-l0_gal)/lstep_gal else $
          offset = (alog10(emission_setup[j].lambda)-l0_gal)/lstep_gal
        start_pars[h+0] = sol_i[0] 
        start_pars[h+1] = sol_i[1]/velscale + offset ; back to pixels positions
        start_pars[h+2] = sol_i[2]/velscale          ; back to pixels widths
        h = h + 3
    endfor
    ; If needed, add the starting reddening guesses
    if (n_elements(reddening) ne 0) then $
      start_pars[3*nlines:3*nlines+n_elements(reddening)-1] = sol[4*nlines:4*nlines+n_elements(reddening)-1]
    ; If needed, add the starting mult. polynomial coefficients
    ; which are at the end of the sol solution array
    if mdegree ne 0 then $
      start_pars[n_elements(start_pars)-mdegree:n_elements(start_pars)-1] = $
      sol[n_elements(sol)-mdegree:n_elements(sol)-1]

    ; -----------------
    ; Set up again the limits and the appropriate inter-dependencies
    ; for the parameters to be fitted. This time use the FOR_ERRORS 
    ; keyword in SET_CONSTRAINTS.
    set_constraints, GALAXY=galaxy, NOISE=noise, CSTAR=cstar, KINSTARS=kinstars, $
      VELSCALE=velscale, DEGREE=degree, MDEGREE=mdegree, GOODPIXELS=goodpixels, $
      EMISSION_SETUP=emission_setup, START_PARS=start_pars, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal, $
      PARINFO=parinfo, FUNCTARGS=functargs, INT_DISP=int_disp, $
      LOG10=log10, REDDENING=reddening, L0_TEMPL=l0_templ, $
      FOR_ERRORS=for_errors, vmaxshift=vmaxshift, sigmamax=sigmamax

    ; -----------------
    ; Re-run MPFIT starting from previous solution and using now the
    ; FOR_ERRORS keyword to specify that we solve non-linearly also for
    ; the amplitudes, and not only for the line position and width.
    best_pars_2 = mpfit('FITFUNC_GAS',start_pars, FUNCTARGS=functargs, PARINFO=parinfo, $
                        FTOL=1d-1, $
      NFEV=ncalls, ERRMSG=errmsg, PERROR=errors_2, STATUS=status, QUIET=1)
    ; -----------------
    ; Re-evaluate the fit residuals to re-assess the fit quality and
    ; rescale the errors. The last MPFIT fit should have always
    ; converged to the input best solution Also, the weights here are
    ; set to unity for the emission-line templates, as their amplitude
    ; is determined by MPFIT.
; jm10arp28ucsd - KINSTARS=sol[0:5] should be KINSTARS=KINSTARS!!!
    resid_2 = fitfunc_gas(best_pars_2,CSTAR=cstar, GALAXY=galaxy, NOISE=noise, $
                          KINSTARS=kinstars, VELSCALE=velscale, DEGREE=degree, MDEGREE=mdegree, $
                          GOODPIXELS=goodpixels, BESTFIT=bestfit_2, WEIGHTS=weights_2, $
                          EMISSION_SETUP=emission_setup, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal, $
                          EMISSION_TEMPLATES=emission_templates_2, INT_DISP=int_disp, $
                          LOG10=log10, REDDENING=reddening, L0_TEMPL=l0_templ, $
                          FOR_ERRORS=for_errors)
    ;

    if total(noise) eq n_elements(galaxy) then begin
        ; If as errors on the fluxes you have input an array of
        ; constant unity vales compute Chi^2/DOF and use this instead
        ; of bestnorm/dof to rescale the formal uncertainties
;        print,'rescaling...'
        chi2_2   = im_robust_sigma(resid_2, /ZERO)^2 
        errors_2 = errors_2*sqrt(chi2_2)
    endif

    ; -----------------
    ; Add up the best-fitting emission templates to get the emission spectrum
    if ((size(emission_templates_2))[0] eq 1) then emission_2 = emission_templates_2          
    if ((size(emission_templates_2))[0] eq 2) then emission_2 = total(emission_templates_2,2)

    ; -----------------
    ; Rearrange the final results in the output array SOL, which
    ; includes also line fluxes. This time evaluate also the errors on
    ; these last values, using for now a simple MC error propagation
    rearrange_results, best_pars_2, weights_2, chi2_2, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal, $
      VELSCALE=velscale, EMISSION_SETUP=emission_setup, SOL=sol_2, $
      ERR=errors_2, ESOL=esol_2, INT_DISP=int_disp, $
      LOG10=log10,REDDENING=reddening, $
      FOR_ERRORS=for_errors
    ; Appends to the best-fitting gas results also the mdegree polynomial
    ; coefficients.
    if mdegree ne 0 then sol_2 = [sol_2,best_pars_2[n_elements(best_pars_2)-mdegree:n_elements(best_pars_2)-1]]

    ; -----------------
    ; Rewrite on the final solution array
    sol = sol_2           & esol = esol_2
    best_pars = best_pars & bestfit = bestfit_2 
    emission = emission_2 & emission_templates = emission_templates_2
    weights = weights_2

    ; -----------------
    ; Show the fit if requested
    if keyword_set(plot) then show_fit, galaxy, bestfit, emission, best_pars, sol, $
      GOODPIXELS=goodpixels, MDEGREE=mdegree, $
      REDDENING=reddening, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal, L0_TEMPL=l0_templ, LOG10=log10, $
      KINSTARS=kinstars, NLINES=nlines
    
ENDIF

; jm09nov19ucsd - this is the wrong thing to do, so ignore

;;; ------------------------------------
;;; If we used reddening, recover the reddened amplitudes.  In other
;;; words, make sure we output the emission-line amplitudes as observed.
;;; This is the right thing to later compute the amplitude-over-noise
;;; ratio of the lines and decide whether they are detected
;;if (n_elements(reddening) ne 0) then begin
;;    ; make the spectrum wavelength array
;;    if not keyword_set(log10) then $
;;      ob_lambda = exp(dindgen(n_elements(galaxy))*lstep_gal + l0_gal) else $
;;      ob_lambda = 10^(dindgen(n_elements(galaxy))*lstep_gal + l0_gal)
;;    ; receding velocity
;;    c = 299792.458d
;;    if not keyword_set(log10) then $ 
;;      Vstar = kinstars[0] + (l0_gal-l0_templ)*c else $
;;      Vstar = kinstars[0] + (l0_gal-l0_templ)*c*alog(10.0d)
;;    ; total reddening attenuation that was applied to the emission lines
;;    ; in FITFUNC_GAS
;;    reddening_attenuation = $
;;      dust_calzetti(l0_gal,lstep_gal,n_elements(galaxy),sol[nlines*4],Vstar,log10=log10)        
;;    if n_elements(reddening) eq 2 then begin
;;        int_reddening_attenuation = $
;;          dust_calzetti(l0_gal,lstep_gal,n_elements(galaxy),sol[nlines*4+1],Vstar,log10=log10)
;;    endif else int_reddening_attenuation = 1.0 
;;    ; get the reddening attenuation at the line wavelength 
;;    rf_lambda = ob_lambda/exp(Vstar/c)
;;    reddening_attenuation_emission = $
;;      interpol(reddening_attenuation*int_reddening_attenuation,rf_lambda,emission_setup.lambda[i_lines])
;;    l = indgen(nlines)
;;    ; Finally, attenuate the output amplitude
;;    sol[l*4+1] = sol[l*4+1]*reddening_attenuation_emission
;;    ; and corresponding errors on the amplitudes
;;    if keyword_set(for_errors) then esol[l*4+1] = esol[l*4+1]*reddening_attenuation_emission
;;endif

; ------------------------------------
if (keyword_set(quiet) eq 0) then begin
    resid = galaxy - bestfit
    resid_noise = im_robust_sigma(resid[goodpixels], /ZERO)
    if not keyword_set(reddening) then print,'Line','Flux','Ampl.','V   ','sig','A/N',f='(6A10)'
    if keyword_set(reddening)     then print,'Line','Deredd. Flux','Ampl.','V   ','sig','A/N',f='(A10,A14,A7,A9,A9,A12)'
    l = indgen(nlines)
    forprint,emission_setup[i_lines].name,sol[l*4],sol[l*4+1],sol[l*4+2],sol[l*4+3],$
      sol[l*4+1]/resid_noise,textout=2,f='(A10,2F10.2,2F10.1,F10.2)'
    if keyword_set(reddening) then begin
        print,'E(B-V)=',sol[nlines*4],f='(A10,F10.2)'
        if n_elements(reddening) eq 2 then print,'E(B-V)_int=',sol[nlines*4+1],f='(A10,F10.2)'
    endif
    if keyword_set(for_errors) then begin
        print,['']
        print,['=========================== Errors =========================']
        print,['']
        l = indgen(nlines)
        forprint,emission_setup[i_lines].name,esol[l*4],esol[l*4+1],esol[l*4+2],esol[l*4+3],$
          esol[l*4+1]/resid_noise,textout=2,f='(A10,2F10.2,2F10.1,F10.2)'
        if keyword_set(reddening) then begin
            print,'E(B-V)',esol[nlines*4],f='(A10,F10.3)'
            if n_elements(reddening) eq 2 then print,'E(B-V)_int',esol[nlines*4+1],f='(A10,F10.3)'
        endif
        print, 'RoN =', resid_noise,f='(A10,F10.1)'
    endif

    print, 'feval', ncalls, FORMAT='(A10,I10)'
endif
; ------------------------------------

; Restore the input emission-line setup structure
emission_setup = dummy

END
;------------------------------------------------------------------------------------
