pro bcgimf_a2261
; jm11oct31ucsd -

; a2261_bcg_apmag_total.txt = aperture photometry of 2D model fits with aperture radius = 1 arcmin.
; a2261_bcg_apmag_3arcsec.txt = aperture photometry within 3 arcsec diameter to match SDSS.
; a2261_bcg_apmag_9arcsec.txt = aperture photometry within 8.85 arcsec (295 pixel) radius.

; X-ray analysis (from Megan):
; M2500 = 2.9 +/- 0.5 e14 Msun
;  
; Also attached is a profile from an isothermal model fit to the SDSS
; BCG velocity dispersion of v = (388 +/- 19) km/s (measured by
; Claudio). 
; 
; Rines10 dynamical mass estimates:
; M100 = 7.13 +/- 0.83 e14 Msun (virial estimate)
; M100 = 5.10 +/- 2.07 e14 Msun (caustic estimate)
; 
; z = 0.224
; 1" ~ 3.59 kpc
; Mvir ~ M115
    
    imfpath = clash_path(/bcgimf)    
    isedpath = clash_path(/ised)    
    scale = 3.59 ; [kpc/arcsec]

    ised = mrdfits(isedpath+'bcgs_fsps_chab_charlot_sfhgrid10.fits.gz',1)
    
    psfile = imfpath+'bcgimf_a2261.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.3
    djs_plot, [0], [0], position=pos, xsty=1, ysty=1, $
      xrange=[1,1000], yrange=[9.5,15], /xlog, xtitle='Radius !8r!6 (kpc)', $
      ytitle='Mass within Sphere M(<!8r!6) (M'+sunsymbol()+')'
    im_legend, 'Abell 2261', /left, /top, box=0, charsize=1.8

    rr = [500,800,800,500] & mm = [10.3,10.3,10.5,10.5]
    polyfill, rr, mm, /data, color=im_color('midnight blue'), noclip=0, /fill
    polyfill, rr, mm-0.4, /data, color=im_color('orange red'), noclip=0, /fill
    xyouts, rr[0]-100, mm[0], 'Strong+Weak Lensing', /data, align=1, charsize=1.6
    xyouts, rr[0]-100, mm[0]-0.4, 'X-ray + HSE', /data, align=1, charsize=1.6
    
; X-ray analysis
    xray = rsex(imfpath+'MNFW_Xray.dat')
;   keep = where(xray.rkpc lt 600)
;   xray = xray[keep]
    polyfill, [xray.rkpc,reverse(xray.rkpc)], alog10([xray.mlo,reverse(xray.mhi)]), $
      /data, color=im_color('orange red'), noclip=0, /fill

; weak + strong lensing mass    
    nfw = rsex(imfpath+'MNFW-5.dat')
    polyfill, [nfw.rkpc,reverse(nfw.rkpc)], alog10([nfw.mlo,reverse(nfw.mhi)]), $
      /data, color=im_color('midnight blue'), noclip=0, /fill

;; BCG analysis
;    xray = rsex(imfpath+'MNFW_Xray.dat')
;    polyfill, [xray.rkpc,reverse(xray.rkpc)], alog10([xray.mlo,xray.mhi]), $
;      /data, /color=im_color('dodger blue'), noclip=0, /fill
    
    djs_oplot, [1.5*scale], [ised[0].mass_50], psym=symcat(15,thick=5), $
      symsize=4, color=im_color('dodger blue')
    djs_oplot, [1.5*scale], [ised[0].mass_50], psym=symcat(6,thick=10), $
      symsize=4, color=im_color('black')
    djs_oplot, [1.5*scale], [ised[0].mass_50]+0.25, psym=symcat(6,thick=10), $
      symsize=4, color=im_color('black')

    djs_oplot, [8.85*scale], [ised[1].mass_50], psym=symcat(16,thick=5), $
      symsize=4, color=im_color('dodger blue')
    djs_oplot, [8.85*scale], [ised[1].mass_50], psym=symcat(9,thick=10), $
      symsize=4, color=im_color('black')
    djs_oplot, [8.85*scale], [ised[1].mass_50]+0.25, psym=symcat(9,thick=10), $
      symsize=4, color=im_color('black')

    djs_oplot, [60.0*scale], [ised[2].mass_50], psym=symcat(14,thick=5), $
      symsize=4.5, color=im_color('dodger blue')
    djs_oplot, [60.0*scale], [ised[2].mass_50], psym=symcat(4,thick=10), $
      symsize=4.5, color=im_color('black')
    djs_oplot, [60.0*scale], [ised[2].mass_50]+0.25, psym=symcat(4,thick=10), $
      symsize=4.5, color=im_color('black')

;   djs_oploterr, 3.0*scale/2.0, ised[0].mass_50, yerr=ised[0].mass_err, $
;     psym=symcat(6,thick=5), symsize=3.5
;   djs_oploterr, 8.85*scale/2.0, ised[1].mass_50, yerr=ised[1].mass_err, $
;     psym=symcat(16,thick=5), symsize=3.5
;   djs_oploterr, 60.0*scale/2.0, ised[2].mass_50, yerr=ised[2].mass_err, $
;     psym=symcat(14,thick=5), symsize=3.5

    im_plotconfig, psfile=psfile, /psclose

return
end
    
