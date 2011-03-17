;+
; NAME:
;   long_reduce_work
;
; PURPOSE:
;
;   Main program for the Low-redux pipeline.  This set of algorithms
;   runs mainly as a black box.
;
; CALLING SEQUENCE:
;  long_reduce, planfile, /clobber, /NOZAP, /NOFLEX, /NOHELIO 
;
; INPUTS:
;  planfile  -- File created by long_plan which guides the reduction
;               process
;
; OPTIONAL INPUTS:
; /NOFLEX  -- Do not apply flexure correction [necessary if your setup
;             has not been calibrated.  Contact JH or JXP for help if
;             this is the case.]
;  HAND_FWHM -- Set the FWHM of the object profile to this value (in
;  PROF_NSIGMA= -- Extend the region to fit a profile by hand 
; /NOHELIO -- Do not correct to heliocentric velocities
; /NOZAP   -- Do not flag CRs
;
; OUTPUTS:
;  (1) Various calibration files
;  (2) One multi-extension FITS file in Science per exposure containing
;  the extracted data and processed images
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;   
; REVISION HISTORY:
;   11-Mar-2005  Written by JH + SB
;-  
;-----------------------------------------------------------------------------
; BUGS:
;   CR zapping should be improved???
;
;-----------------------------------------------------------------------------
; The code that follows is the science frame reduction code

PRO atlas_long_ccdproc, filename, scifile, ISLIT = ISLIT $
                      , slitfile = slitfile, wavefile = wavefile $
                      , biasfile = biasfile, pixflatfile = pixflatfile $
                      , illumflatfile = illumflatfile $
                      , verbose = verbose, box_rad = box_rad1 $
                      , maxobj = maxobj $
                      , reduxthresh = reduxthresh, SIG_THRESH=sig_thresh $
                      , NOFLEX = noflex1, NOSHIFT = NOSHIFT1 $
                      , NOHELIO = NOHELIO , MAXFLEX=maxflex1 , MAXGOOD=maxgood $
                      , NOZAP = NOZAP, HAND_X = HAND_X, HAND_Y = HAND_Y $
                      , HAND_FWHM = HAND_FWHM, STD = STD, FILESTD = FILESTD1 $
                      , CHK = CHK, TRCCHK = TRCCHK, SKYTRACE = SKYTRACE1 $
                      , PROF_NSIGMA=prof_nsigma1, NOLOCAL=nolocal, $
  niter=niter, novac=novac, badpixfile=badpixfile ; jm09dec18ucsd

t0 = systime(1)

;if (NOT keyword_set(profile_filename)) then profile_filename = 0
if (NOT keyword_set(maxobj)) then maxobj = 0
  
;;----------
;; Read the raw science image
long_proc, filename, sciimg, sciivar, hdr = scihdr $
           , biasfile = biasfile, pixflatfile = pixflatfile, bin = bin
;;  What telescope are we using?
telescope = strcompress(sxpar(scihdr[*, 0], 'TELESCOP'), /rem)
;;  Determine detector specific parameters for reduction
long_reduce_params, scihdr, bin, skyfile = skyfile, anamorph = anamorph $
                    , bsp = bsp, SN_GAUSS = SN_GAUSS, SKYTRACE = SKYTRACE $
                    , NOSHIFT = NOSHIFT, NCCD = NCCD, PEAK_SMTH = PEAK_SMTH $
                    , FILESTD = FILESTD, FWHM = FWHM, BOX_RAD = BOX_RAD

; repair bad pixels
if (n_elements(badpixfile) ne 0) then begin
   if (file_test(badpixfile) eq 0) then message, $
     'Bad pixel file '+badpixfile+' not found'

   splog, 'Reading '+badpixfile
   readcol, badpixfile, x1, x2, y1, y2, format='L,L,L,L', $
     /silent, comment='#'
;  splog, 'Identified '+strn(n_elements(x1))+$
;    ' rows in the bad pixel file '+badpixfile
; IDL is zero-indexed but the bad pixel file is one-indexed (IRAF
; convention)
   dims = size(sciimg,/dim)
   badpixmask = sciimg*0.0
   for ii = 0, n_elements(x1)-1 do $
     badpixmask[(y1[ii]-1)>0:(y2[ii]-1)<(dims[0]-1),$
     (x1[ii]-1)>0:(x2[ii]-1)<(dims[1]-1)] = 1
   sciimg = djs_maskinterp(sciimg,badpixmask,iaxis=1,/const)
   sciivar = djs_maskinterp(sciivar,badpixmask,iaxis=1,/const) ; wrong!
endif


;; override param default if standard is set
IF KEYWORD_SET(FILESTD1) THEN FILESTD = FILESTD1
IF KEYWORD_SET(FILESTD) THEN BEGIN
    splog, 'Using standard star trace as crutch from ' + filestd
    stdstruct = xmrdfits(filestd, 5, /silent)
    stdmax = max(stdstruct.PEAKFLUX, stdind)
    stdtrace = stdstruct[stdind].XPOS
ENDIF ELSE stdtrace = 0

IF KEYWORD_SET(STD) THEN BEGIN
   NOFLEX = 1
   SKYTRACE = 0
   box_rad = 50
   NOSHIFT = 1
   NOLOCAL=1
   NITER = 1
   maxobj = 1
   sigrej = 25.0
   PROF_NSIGMA = 20.0
;;  Over-ridden by JXP
;   reduxthresh=0.001 ;; makes sure we reduce other objects on the slit. 
ENDIF ELSE BEGIN
   ;; If any of these were passed in, overwrite long_reduce_params values
   IF n_elements(skytrace1) GT 0 THEN SKYTRACE = SKYTRACE1
   IF n_elements(noshift1)  GT 0 THEN NOSHIFT = NOSHIFT1
   IF n_elements(noflex1) GT 0 THEN NOFLEX = NOFLEX1
   IF n_elements(maxflex1) GT 0 THEN MAXFLEX = MAXFLEX1
   IF n_elements(box_rad1)  GT 0 THEN BOX_RAD = BOX_RAD1
   IF KEYWORd_SET(PROF_NSIGMA1) THEN PROF_NSIGMA=PROF_NSIGMA1
ENDELSE

dims = size(sciimg, /dimens)
nx = dims[0]
ny = dims[1]

    ;; Read in slitmask structure
tset_slits = xmrdfits(slitfile, 1, silent = (keyword_set(verbose) EQ 0))
;; FLEXURE SHIFT SLITS FOR LRIS
IF NOT KEYWORD_SET(NOSHIFT) THEN $
  xshift = long_xcorr_slits(sciimg, tset_slits, /shift)
;;   Generate slit position, wavelengths, and slit illumination 
ximg = long_slits2x(tset_slits, edgmask = edgmask, nslit = nslit)
slitmask = long_slits2mask(tset_slits) * (ximg GT 0. AND ximg LT 1.)
IF strcmp(telescope, 'Gemini-North') THEN BEGIN
    IF KEYWORD_SET(illumflatfile) THEN $
      slit_illum = xmrdfits(illumflatfile, silent = (keyword_set(verbose) EQ 0))
    waveimg = xmrdfits(wavefile, silent = (keyword_set(verbose) EQ 0), 0)
    piximg = xmrdfits(wavefile, silent = (keyword_set(verbose) EQ 0), 0)
ENDIF ELSE BEGIN
   ;; TESTING
   ;IF KEYWORD_SET(illumflatfile) THEN $
   ;   slit_illum = long_slitillum(illumflatfile, slitmask, ximg, edgmask)
    ;;   Reconstruct PIXIMG WAVEIMG and using coefficient sets
   slitillum = 0.0*ximg + 1.0D
    pixset  = xmrdfits(wavefile, silent = (keyword_set(verbose) EQ 0), 1)
    wavesvfile =  repstr(wavefile, '.fits', '.sav')
    restore, wavesvfile
    piximg = long_wpix2image(pixset, tset_slits, XFIT = xfit $
                             , waveimg = waveimg)
ENDELSE
;   Read in wavelength solution structures
fwhmset = xmrdfits(wavefile, silent = (keyword_set(verbose) EQ 0), 2)

;--------- 
;  If there is no slit illumination function assume that it is unity
IF NOT KEYWORD_SET(slit_illum) THEN slit_illum = 0.0*sciimg +1.0
splog, 'Applying slit illumination'
; Allocate images which we will need later
objimage = fltarr(nx, ny)
;; Don't apply illumination corrections larger than 30%
gdpix = WHERE(slit_illum GT 0.6D AND slit_illum LT 1.4)
sciimg[gdpix]  = sciimg[gdpix]/slit_illum[gdpix]
sciivar[gdpix] = sciivar[gdpix]*slit_illum[gdpix]^2 

;; Trace sky lines for sky subtraction?
IF KEYWORD_SET(SKYTRACE) THEN BEGIN
    wavesvfile =  repstr(wavefile, '.fits', '.sav')
    restore, wavesvfile
    wstruct = long_wstruct(scihdr)
    pixset = long_wavepix(sciimg, tset_slits, fwhm = fwhmset.MEDIAN $
                          , box_radius = wstruct.radius $
                          , sig_thresh = wstruct.sig_wpix $
                          , pkwdth = wstruct.pkwdth $
                          , TOLER = wstruct.TOLER, CHK = TRCCHK $
                          , piximg_in = piximg)
    piximg = long_wpix2image(pixset, tset_slits, XFIT = xfit $
                             , waveimg = waveimg)
ENDIF
splog, 'Finding objects on the slits: First pass'
;----------
;  Find objects on the slits
objstruct1 = long_objfind(sciimg, tset_slits = tset_slits, invvar = sciivar $
                          , skymask = skymask1, objmask = objmask1 $
                          , nperslit = maxobj, peakthresh = reduxthresh $
                          , fwhm = FWHM, PEAK_SMTH = PEAK_SMTH $
                          , SIG_THRESH = SIG_THRESH $
                          , HAND_X = HAND_X, HAND_Y = HAND_Y $
                          , HAND_FWHM = HAND_FWHM, STDTRACE = STDTRACE $
                          , ISLIT = ISLIT)
splog, 'Aperture masked sky subtraction'
skyimage = long_skysub(sciimg, sciivar, piximg, slitmask, skymask1, edgmask $
                       , bsp = bsp, ISLIT = ISLIT)
IF NOT KEYWORD_SET(NOZAP) THEN BEGIN
    splog, 'Cosmic ray rejection'
    IF KEYWORD_SET(FWHMSET) THEN sigma_psf = $
      djs_median(fwhmset.median)/2.35482D $
    ELSE sigma_psf = 3.0D/2.35482D
    ;;   Take the PSF width to be that of the spectral direction.  
    ;;   This prevents the routine from rejecting sky lines
    ;;   This is a description of the 3x3 core of the 2D PSF for reject_cr.pro
    ;;    
    ;;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1] 
    ;;    PSFVALS[0]          1.   PSFVALS[0]
    ;;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1]
    psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
    crmask  = psf_reject_cr(sciimg-skyimage, sciivar, psfvals $
                            , satmask = (sciimg-skyimage) GT 8d4)
    sciivar = sciivar*(crmask EQ 0)
    ;; Do we need to sky-subtract before CR rejection???
 ENDIF
splog, 'Finding objects in sky-subtracted image: Second pass'
; Redo object finding on sky subtracted image
;IF KEYWORD_SET(OBJSTRUCT1) THEN FWHM = djs_median(objstruct1.FWHM)
objstruct = long_objfind(sciimg-skyimage, tset_slits = tset_slits $
                         , invvar = sciivar, skymask = skymask $
                         , objmask = objmask, nperslit = maxobj $
                         , peakthresh = reduxthresh $
                         , SIG_THRESH = SIG_THRESH $
                         , fwhm = FWHM, PEAK_SMTH = PEAK_SMTH $
                         , HAND_X = HAND_X, HAND_Y = HAND_Y $
                         , HAND_FWHM = HAND_FWHM, STDTRACE = STDTRACE $
                         , ISLIT = ISLIT)
splog, 'Redoing global sky subtraction'
skyimage = long_skysub(sciimg, sciivar, piximg, slitmask, skymask, edgmask $
                       , bsp = bsp, ISLIT = ISLIT)

;; Keep spatial flexure shift in objstruct
IF NOT KEYWORD_SET(NOSHIFT) AND KEYWORD_SET(xshift) $
  THEN objstruct.FLX_SHFT_SPA = xshift
final_struct = 0

;----------
; Write output file

modelivar = sciivar ; jm09dec28ucsd; note!
outmask = objmask
final_struct = objstruct

splog, 'Writing FITS file ', scifile
mwrfits, float(sciimg), scifile, scihdr[*, 0], /create
mwrfits, float(modelivar)*float(slitmask GT 0), scifile
mwrfits, float(skyimage)*float(slitmask GT 0), scifile
mwrfits, float(objimage)*float(slitmask GT 0), scifile
mwrfits, float(outmask)*float(slitmask GT 0), scifile
mwrfits, final_struct, scifile

splog, 'Compressing ', scifile
spawn, 'gzip -f '+scifile
;long_plotsci, scifile, hard_ps = repstr(scifile, '.fits', '.ps')

splog, 'Elapsed time = ', systime(1)-t0, ' sec'

RETURN
END
