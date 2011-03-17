pro stelib_templates, starflux, write=write, doplot=doplot
; jm03may7uofa - adopted from an earlier version

; i would like to generate bins of metallicity and assign each bin to
; a different FITS extension, sequenced by effective temperature
    
; also, make sure to reject Wolf-Rayet stars and also stars that do
; not have information at every pixel in the wavelength range of
; interest
    
    starinfo = read_stelib(/silent)
    starcount = n_elements(starinfo)

    keep = where(starinfo.teff gt -900.0,starcount)
    starinfo = starinfo[keep]
    srt = reverse(sort(starinfo.teff))
    starinfo = starinfo[srt]

    print_struct, starinfo, ['STAR','TYPE','TEFF','LOGG','FE_H']

; bin the STELIB stars by metallicity:  solar, metal_poor, and 
; super_metal_poor
    
    bin1 = where(starinfo.fe_h ge 0.0,nbin1)
    bin2 = where((starinfo.fe_h ge -1.0) and (starinfo.fe_h lt 0.0),nbin2)
    bin3 = where(starinfo.fe_h lt -1.0,nbin3)
    
;; select a subset of stars - generate 20 bins of effective temperature
;; see Allen p. 388
;
;; choose main sequence BAF stars and GKM giants and no supergiants
;
;;   plothist, starinfo.teff, bin=1000, xrange=[2000,11000], xsty=3, ysty=3
;  
;    tminmax = [$
;      [10000.0,30000.0], $
;      [8000.0,9999.0], $  
;      [7000.0,7999.0], $  
;      [6000.0,6999.0], $  
;      [5000.0,5999.0], $
;      [4000.0,4999.0], $
;      [3000.0,3999.0] $
;      ]
;    nbins = n_elements(tminmax[0,*])
;    nspecbin = 3L ; number of spectra per bin
;
;    indx = lonarr(nbins*nspecbin)
;    for j = 0L, nbins-1L do begin
;
;       w = where((starinfo.teff gt tminmax[0,j]) and (starinfo.teff lt tminmax[1,j]),nw)
;       rindx = floor(randomu(seed,nspecbin)*nw)
;       indx[nspecbin*j:nspecbin*(j+1)-1L] = w[rindx]
;       
;    endfor
;
;;   print, indx
;    starinfo = starinfo[indx]
;    starcount = n_elements(starinfo)
    
; compute the color excess, E(B-V)

;   glactc, 15*im_hms2dec(starinfo.ra), im_hms2dec(starinfo.dec), 2000.0, gl, gb, 1, /degree
;   starinfo.ebv = dust_getval(gl,gb,/interp)

    minwave = 3200.0
    maxwave = 9000.0
    disp = 1.0       ; Angstrom per pixel

    refwave = starinfo[0].wave
    npix = starinfo[0].npix
    starflux = fltarr(npix,starcount)
    
    for k = 0L, starcount-1L do begin

       flux = starinfo[k].spec
       newflux = flux
;      ccm_unred, refwave, flux, starinfo[k].ebv, newflux
       
       starflux[*,k] = newflux/max(newflux)

       if keyword_set(doplot) then begin
          djs_plot, refwave, starflux[*,k], xsty=3, ysty=3, charsize=2.0, charthick=2.0, $
            xtitle='Wavelength ('+angstrom()+')', ytitle='Relative f_{\lambda}', $
            xthick=2.0, ythick=2.0, yrange=[-0.1,1.1]
          legend, [starinfo[k].star,starinfo[k].type,'T = '+string(starinfo[k].teff,format='(F7.1)'),$
            '[Fe/H] = '+string(starinfo[k].fe_h,format='(F5.2)'),'log g = '+$
            string(starinfo[k].logg>(-1),format='(F4.2)')], $
            /right, /top, box=0, charsize=2.0, charthick=2.0
          cc = get_kbrd(1)
       endif
       
    endfor

; generate the FITS header
    
    crpix1 = 0L
    crval1 = minwave
    cd1_1 = disp
       
    mkhdr, header, starflux
    sxaddpar, header, 'CRPIX1', crpix1
    sxaddpar, header, 'CRVAL1', crval1
    sxaddpar, header, 'CD1_1', cd1_1
;   sxaddpar, header, 'NSTAR', starcount

    for j = 0L, starcount-1L do sxaddpar, header, 'STAR'+string(j,format='(I2.2)'), $
      strn(starinfo[j].star)
    for j = 0L, starcount-1L do sxaddpar, header, 'TYPE'+string(j,format='(I2.2)'), $
      strn(starinfo[j].type)
    
; write out

    if keyword_set(write) then begin
    
       tpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='templates')
       tfile = 'stelib_templates.fits'
       print, 'Writing '+tpath+tfile+'.'
       mwrfits, starflux, tpath+tfile, header, /create

    endif
    
return
end
