;+
; NAME:
;       KPNO_ATLAS_SKYSPEC
;
; PURPOSE:
;       Generate an average sky spectrum.
;
; CALLING SEQUENCE:
;
; INPUTS:
;       speclist - list of one-dimensional input spectra
;       
;
; OPTIONAL INPUTS:
;       datapath - full path name to SPECLIST (default PWD)
;       skyfile  - file with a list of sky lines with which to shift
;                  the sky spectra to zero velocity (default
;                  SKYLINES.DAT)
;       skypath  - corresponding path name for SKYFILE (default to
;                  ${ISPEC_DIR}/etc/)
;       wave1    - output starting wavelength for SKYSPEC (Angstrom) 
;       wave2    - output ending wavelength for SKYSPEC (Angstrom)
;       dwave    - output dispersion for SKYSPEC (Angstrom/pixel)
;       normwave - wavelength around which to normalize the sky
;                  spectra (default to the mean wavelength of all the
;                  spectra)
;       normwidth - normalize about NORMWAVE +/- NORMWIDTH
;       title    - title for the plot and the header history note
;       psname   - output root name for the postscript output and for
;                  the FITS file output of SKYSPEC if WRITE has been
;                  set 
;
; KEYWORD PARAMETERS:
;       write    - generate postscript output and write a FITS file
;                  (see PSNAME)
;
; OUTPUTS:
;       skyspec - mean sky spectrum
;       header  - corresponding FITS header for SKYSPEC
;
; OPTIONAL OUTPUTS:
;       See keyword WRITE (above).
;
; PROCEDURES USED:
;       GET_ELEMENT, DJS_ITERSTAT, SPLOG, READCOL, MPFITPEAK(),
;       DJS_MEDIAN(), IFORAGE(), RD1DSPEC(), COMBINE1FIBER, DFPSPLOT,
;       DFPSCLOSE, MWRFITS
;
; INTERNAL ROUTINES:
;       SPEC_NORMALIZE(), FIT_SKYLINE()
;
; COMMENTS:
;
;
; EXAMPLE:
;       IDL> speclist = ['fdwcra.0001.fits','fdwcra.0002.fits']
;       IDL> imakeskyspec, speclist, psname='kpno_sky', /write, $
;       IDL>  title='KPNO Sky Spectrum'
;
; DATA FILES:
;       ${ISPEC_DIR}/etc/skylines.dat
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 December 4, U of A
;-

function spec_normalize, wave, spec, normwave, normwidth, norm=norm

    get_element, wave, normwave+[-normwidth/2.0,+normwidth/2.0], xx
    
    djs_iterstat, spec[xx[0]:xx[1]], sigrej=2.0, median=norm
    newspec = spec / norm
    
return, newspec
end
    
function fit_skyline, wave, spec, linewave, medshift=medshift
; fit a Gaussian to the sky lines to zero the wavelength scale (should
; be very small shifts)

    nline = n_elements(linewave)
    shift = fltarr(nline)

    for k = 0L, nline-1L do begin
       
       get_element, wave, linewave[k], y12
    
       local = where((wave gt wave[y12-10]) and (wave lt wave[y12+10]),nlocal)
       gsfit = mpfitpeak(double(wave[local]),spec[local],a,nterm=4,/gaussian,/positive)
       shift[k] = a[1]-linewave[k]
       
    endfor

    medshift = djs_median(shift)
    newwave = wave + medshift

;   plot, newwave[local], skyvector[local,0], ps=10, xsty=3, ysty=3
;   djs_oplot, newwave[local], gsfit, color='red', ps=10

return, newwave
end    

pro imakeskyspec, speclist, skyspec, header, datapath=datapath, skyfile=skyfile, $
  skypath=skypath, wave1=wave1, wave2=wave2, dwave=dwave, normwave=normwave, $
  normwidth=normwidth, title=title, psname=psname, write=write

    nspec = n_elements(speclist)
    if nspec eq 0L then begin
       print, 'Syntax - imakeskyspec, speclist, skyspec, header, '
       return
    endif

    if n_elements(datapath) eq 0L then datapath = cwd()
    pushd, datapath

    if n_elements(skypath) eq 0L then skypath = $
      filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')
    if n_elements(skyfile) eq 0L then skyfile = 'skylines.dat'

;;    forage = iforage(speclist)
;;    npix = min(forage.naxis1)
;;    if n_elements(dwave) eq 0L then dwave = forage[0].cd1_1
;;    if n_elements(wave1) eq 0L then wave1 = max(forage.crval1)
;;    if n_elements(wave2) eq 0L then wave2 = dwave*npix+wave1
;;    skywave = wave1+dwave*findgen(npix)

    dwave = 2.75
    wave1 = 3800.0
    wave2 = 6700.0
    skywave = findgen((wave2-wave1)/dwave)*dwave+wave1
    npix = n_elements(skywave)
    
    if file_test(skypath+skyfile,/regular) eq 0L then begin
       splog, 'Sky line file '+skypath+skyfile+' not found.'
       return
    endif else splog, 'Reading '+skypath+skyfile+'.'

    readcol, skypath+skyfile, linewaves, /silent, format='D'
    good = where((linewaves gt min(skywave)) and (linewaves lt max(skywave)),nline)
    if nline ne 0L then linewaves = linewaves[good] else begin
       splog, 'No sky lines in '+skypath+skyfile+' are in the wavelength range.'
       return
    endelse
    
    if n_elements(psname) eq 0L then psname = 'sky' ; root name

; could put some checks in here for WAVE1 and WAVE2
    
    splog, 'Reading '+strn(nspec)+' spectra.'

    skyvector = dblarr(npix,nspec)

; choose the normalization wavelengths

    if n_elements(normwave) eq 0L then normwave = djs_mean(skywave) ; Angstrom
    if n_elements(normwidth) eq 0L then normwidth = 200.0           ; Angstrom

    for i = 0L, nspec-1L do begin

       spec = rd1dspec(speclist[i],/silent)
       sky = spec_normalize(spec.wave,spec.sky,normwave,normwidth,norm=norm1)
       skyivar = 1.0/(spec.sigspec/norm1)^2.0

       wave = spec.wave
;;       wave = fit_skyline(spec.wave,sky,linewaves,medshift=medshift1)
;      splog, 'Median shift (Angstrom): '+strn(medshift1)

       combine1fiber, alog10(wave), sky, skyivar, newloglam=alog10(skywave), newflux=newflux
       skyvector[*,i] = newflux

;      skyvector[*,i] = interpol(sky,wave,skywave) ; linearly interpolate

;      if i eq 0L then $
;        djs_plot, skywave, skyvector[*,i], color='red', xsty=3, ysty=3, ps=10 else $
;        djs_oplot, skywave, skyvector[*,i], color='green', ps=10
;      cc = get_kbrd(1)

       print, format='("Reading spectrum ",I3,"/",I3,".",A1,$)', i+1, nspec, string(13b)
       
    endfor
    print
    
; compute the mean and the one-sigma sky spectra
    
    skyspec = total(skyvector,2)/nspec
;   skyspec = djs_median(skyvector,2)

;   skysig = sqrt(total((skyvector-skyspec#replicate(1,nspec))^2.0,2)/(nspec-1))
    skysig = skyspec*0.0
    for j = 0L, npix-1L do begin
       djs_iterstat, skyvector[j,*], sigma=sigma1, sigrej=3.0
       skysig[j] = sigma1
;      skysig[j] = stddev(skyvector[j,*])
    endfor
       
; plot and write out

    if keyword_set(write) then itermax = 1L else itermax = 0L
    etcpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')
    
    for k = 0L, itermax do begin

       if k eq 1L then $
         dfpsplot, etcpath+psname+'.ps', /color, /isolatin1, /square else $
         window, 0, xs=550, ys=550
       
       djs_plot, skywave, skyspec, ps=10, xsty=3, ysty=3, yrange=[0,5], $
         ytitle='f_{\lambda} (arbitrary units)', xtitle='Wavelength (\AA)', $
         charsize=2.0, charthick=2.0, xthick=2.0, ythick=2.0, thick=3.0, title=title
       djs_oplot, skywave, skyspec-3*skysig, line=1, color='green'
       djs_oplot, skywave, skyspec+3*skysig, line=1, color='green'
;      legend, ['KPNO night sky spectrum'], /left, /top, box=0, charsize=2.0, charthick=2.0         
       
       if k eq 1L then dfpsclose

    endfor

; write out

    skyfits = etcpath+psname+'.fits'
    mkhdr, header, skyspec, /extend
    sxaddpar, header, 'CTYPE1', 'LINEAR'
    sxaddpar, header, 'CRPIX1', 1.0
    sxaddpar, header, 'CRVAL1', wave1
    sxaddpar, header, 'CD1_1', dwave
    sxaddpar, header, 'ZUNITS', 'erg/s/cm^2/A', 'normalized flux units'
    sxaddpar, header, 'NCOMBINE', nspec, 'Number of combined spectra'
    if n_elements(title) ne 0L then sxaddhist, title, header

    if keyword_set(write) then begin
    
       splog, 'Writing '+skyfits+'.'

       mwrfits, float(skyspec), skyfits, header, /create
       mwrfits, float(skysig), skyfits

    endif

    popd
    
return
end
