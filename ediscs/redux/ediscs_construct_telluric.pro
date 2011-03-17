;+
; NAME:
;       EDISCS_CONSTRUCT_TELLURIC
;
; PURPOSE:
;       Construct a telluric absorption-line spectrum from
;       one-dimensional telluric standard star spectra.
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;       tlist - 
;
; OPTIONAL INPUTS:
;       datapath  - I/O path
;       pixspace  - 
;       ncoeff    - 
;       sigclip   - sigma-clipping threshold (default 3.0)
;       plottitle -
;       psname    - 
;       tellfile  - 
;
; KEYWORD PARAMETERS:
;       dospline - 
;       debug    - 
;       doplot   - 
;       write    - 
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       tellwave - 
;       tellcorr - 
;
; PROCEDURES USED:
;       RD1DSPEC(), ICLEANUP, SPLOG, BSPLINE_ITERFIT(),
;       BSPLINE_VALU(), FUNC_FIT(), FLEGENDRE(), DJS_REJECT(),
;       DJS_PLOT, DJS_OPLOT, MWRFITS, DFPSPLOT, DFPSCLOSE, ANGSTROM(),
;       REPSTR() 
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 July 12, U of A, written, based entirely on
;          David Schlegel's TELLURIC_CORR
;       jm04sep10uofa - added SIGCLIP optional input
;-

pro ediscs_construct_telluric, tlist, datapath=datapath, bkspacing=bkspacing, $
  plottitle=plottitle, psname=psname, tellwave=tellwave, tellcorr=tellcorr, $
  tellfile=tellfile, _extra=extra, debug=debug, doplot=doplot, write=write

    ntell = n_elements(tlist)
    if (ntell eq 0L) then begin
       print, 'Syntax - ediscs_construct_telluric '
       return
    endif

    splog, 'Using '+string(ntell,format='(I0)')+' telluric standards.'
    
    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(bkspacing) eq 0) then bkspacing = 1.0
    if (n_elements(plottitle) eq 0) then plottitle = ''
    if (n_elements(psname) eq 0) then psname = 'qaplot_telluric_spectrum.ps'
    if (n_elements(tellfile) eq 0) then tellfile = 'telluric_spectrum.fits'

    if keyword_set(write) then begin
       doplot = 0L
       debug = 0L
       postthick = 5.0
    endif else postthick = 2.0

; initialize the telluric spectrum then read all the telluric standard
; star spectra

    tfirst = rd1dspec(tlist[0],datapath=datapath,/silent)
    npix = tfirst.npix

    flux = fltarr(npix,ntell)+1.0
    ivar = fltarr(npix,ntell)
    wave = fltarr(npix,ntell)

    continuum = im_fitcontinuum(tfirst.wave,tfirst.spec,$
      doplot=debug,/silent,_extra=extra)

    telluric_mask = telluric_mask(tfirst.wave,bad=telluric_pixels)
    
    flux[telluric_pixels,0] = tfirst.spec[telluric_pixels]/continuum[telluric_pixels]
    ivar[telluric_pixels,0] = 1.0/(tfirst.sigspec[telluric_pixels]/continuum[telluric_pixels])^2
    wave[*,0] = tfirst.wave
    
    if keyword_set(debug) then cc = get_kbrd(1)
    
    icleanup, tfirst
    
    for i = 1L, ntell-1L do begin

       tnext = rd1dspec(tlist[i],datapath=datapath,/silent)
       if (tnext.npix ne npix) then begin
          splog, 'All the telluric standards must have the same number of pixels!'
          splog, 'Incompatible standard '+datapath+tlist[i]+'.'
          return
       endif

       continuum = im_fitcontinuum(tnext.wave,tnext.spec,$
         doplot=debug,/silent,_extra=extra)

       telluric_mask = telluric_mask(tnext.wave,bad=telluric_pixels)
       
       flux[telluric_pixels,i] = tnext.spec[telluric_pixels]/continuum[telluric_pixels]
       ivar[telluric_pixels,i] = 1.0/(tnext.sigspec[telluric_pixels]/continuum[telluric_pixels])^2.0
       wave[*,i] = tnext.wave

       icleanup, tnext

       if keyword_set(debug) then cc = get_kbrd(1)

    endfor

; select the data points that describe the telluric factor, and sort
; by wavelength.

    indx = where(ivar NE 0)
    isort = sort(wave[indx])
    fitwave = wave[indx[isort]]
    fitflux = flux[indx[isort]]
    fitivar = ivar[indx[isort]]

; now bspline the features in contflux

;   everyn = 1
;   everyn = ntell
;   everyn = ntell*2
;   everyn = ntell/2.0

    tellwave = reform(wave[*,0])
    bkspace = bkspacing*(tellwave[1]-tellwave[0])
    
    tellset = bspline_iterfit(fitwave,fitflux,nord=nord,$
      maxiter=30,lower=lower,upper=upper,invvar=fitivar, $
      everyn=everyn,bkspace=bkspace,rejper=0.1,/silent)
    tellcorr = bspline_valu(tellwave,tellset)

    if keyword_set(doplot) or keyword_set(write) then begin

       if keyword_set(write) then begin
          splog, 'Generating postscript output '+datapath+psname+'.'
          dfpsplot, datapath+psname, /color, /landscape
          colors = ['navy','red']
       endif else colors = ['cyan','red']
       
       indx = where(ivar ne 0L)
       indx = indx[sort(wave[indx])]
       xrange = minmax(wave[indx])
       
       plot, wave[indx], flux[indx], xrange=xrange, xsty=3, ysty=3
         
stop       
       
       if keyword_set(write) then dfpsclose

    endif 

    if keyword_set(write) then begin

       mkhdr, tellhead, tellwave
       sxdelpar, tellhead, 'COMMENT' & sxdelpar, tellhead, 'DATE'
       sxaddpar, tellhead, 'CRVAL1', tellwave[0], ' central wavelength of first pixel'
       sxaddpar, tellhead, 'CD1_1', tellwave[1]-tellwave[0], ' dispersion [Angstrom/pixel]'
       sxaddpar, tellhead, 'CRPIX1', 1, ' starting pixel (1-indexed)'
       sxaddpar, tellhead, 'CTYPE1', 'LINEAR'
       sxaddpar, tellhead, 'DC-FLAG', 0, ' log-linear flag'
       sxaddhist, "'Telluric correction spectrum generated "+im_today()+"'", tellhead

       splog, 'Writing '+datapath+tellfile+'.'
       mwrfits, tellcorr, datapath+tellfile, tellhead, /create

    endif

stop
    
return
end
    
