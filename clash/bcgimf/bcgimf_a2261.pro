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

    vdisp_sdss = '387.5+/-19 km/s'
    vdisp_gemini = '387+/-16 km/s'
    
    imfpath = getenv('BCGIMF_DATA')+'/'
    isedpath = imfpath
    scale = 3.59 ; [kpc/arcsec]

    phot = mrdfits(isedpath+'bcgimf_phot.fits.gz',1)
    ised1 = mrdfits(isedpath+'bcgimf_fsps_salp_sfhgrid01.fits.gz',1)
    ised2 = mrdfits(isedpath+'bcgimf_fsps_salp_sfhgrid02.fits.gz',1)
    ised3 = mrdfits(isedpath+'bcgimf_bc03_salp_sfhgrid02.fits.gz',1)
    ised4 = mrdfits(isedpath+'bcgimf_fsps_salp_sfhgrid03.fits.gz',1)
;   phot1 = mrdfits(isedpath+'bcgimf_phot_bcgmodel.fits.gz',1)
;   ised1 = mrdfits(isedpath+'bcgmodel_fsps_salp_sfhgrid01.fits.gz',1)

    xrange = [0.3,100]
    yrange = [9.0,14]
    
    psfile = imfpath+'bcgimf_a2261.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.3
    djs_plot, [0], [0], position=pos, xsty=9, ysty=1, $
      xrange=xrange, yrange=yrange, /xlog, xtitle='Radius !8r!6 (kpc)', $
      ytitle='Log Mass within Cylinder M(<!8r!6) (M'+sunsymbol()+')'
;     ytitle='Log Mass within Sphere M(<!8r!6) (M'+sunsymbol()+')'
    axis, /xaxis, xrange=xrange/scale, xsty=1, xtitle='Radius !8r!6 (arcsec)'
;   im_legend, 'Abell 2261', /left, /top, box=0, charsize=1.8
    djs_oplot, 6.0*scale*[1,1], !y.crange

    plots, 30.0, 13.0, psym=8, symsize=3
    plots, 20.0, alog10(5D12), psym=8, symsize=3
    plots, 5.0, alog10(4D11), psym=8, symsize=3

    restore, imfpath+'bcg3d.sav'
    djs_oplot, r, alog10(m>1D9), line=0, color=im_color('blue')
    
    rr = [1.2,2.0,2.0,1.2] & mm = [13.5,13.5,13.7,13.7]
    polyfill, rr, mm, /data, color=im_color('midnight blue'), noclip=0, /fill
    polyfill, rr, mm-0.4, /data, color=im_color('orange red'), noclip=0, /fill
    xyouts, 2.3, mm[0]+0.03, 'Strong+Weak Lensing', /data, align=0, charsize=1.6
    xyouts, 2.3, mm[0]+0.03-0.4, 'X-ray + HSE', /data, align=0, charsize=1.6

; X-ray analysis
    xray = rsex(imfpath+'MNFW_Xray.dat')
;   keep = where(xray.rkpc lt 600)
;   xray = xray[keep]

;   polyfill, [xray.rkpc,reverse(xray.rkpc)], alog10([xray.mlo,reverse(xray.mhi)]), $
;     /data, color=im_color('orange red'), noclip=0, /fill

; weak + strong lensing mass    
    nfw = rsex(imfpath+'MNFW-5.dat')
    polyfill, [nfw.rkpc,reverse(nfw.rkpc)], alog10([nfw.mlo,reverse(nfw.mhi)]), $
      /data, color=im_color('midnight blue'), noclip=0, /fill

;   ww = where(nfw.rkpc lt 10)
;   print, linfit(nfw[ww].rkpc,alog10(nfw[ww].mlo))
    
;; BCG analysis
;    xray = rsex(imfpath+'MNFW_Xray.dat')
;    polyfill, [xray.rkpc,reverse(xray.rkpc)], alog10([xray.mlo,xray.mhi]), $
;      /data, /color=im_color('dodger blue'), noclip=0, /fill

;   djs_oplot, phot.radius*scale, ised.mass_50, line=0, thick=6
;   djs_oplot, phot.radius*scale, ised.mass_50+0.3, line=5, thick=6

    djs_oplot, phot.radius*scale, ised1.mass_50, line=0, thick=6, color='orange'
    djs_oplot, phot.radius*scale, ised1.mass_50+0.3, line=5, thick=6, color='orange'
    djs_oplot, phot.radius*scale, ised2.mass_50, line=0, thick=6, color='green'
    djs_oplot, phot.radius*scale, ised2.mass_50+0.3, line=5, thick=6, color='green'
    djs_oplot, phot.radius*scale, ised3.mass_50, line=0, thick=6, color='red'
    djs_oplot, phot.radius*scale, ised3.mass_50+0.3, line=5, thick=6, color='red'
    djs_oplot, phot.radius*scale, ised4.mass_50, line=0, thick=6, color='purple'
    djs_oplot, phot.radius*scale, ised4.mass_50+0.3, line=5, thick=6, color='purple'
    
;    djs_oplot, [1.5*scale], [ised[0].mass_50], psym=symcat(15,thick=5), $
;      symsize=4, color=im_color('dodger blue')
;    djs_oplot, [1.5*scale], [ised[0].mass_50], psym=symcat(6,thick=10), $
;      symsize=4, color=im_color('black')
;    djs_oplot, [1.5*scale], [ised[0].mass_50]+0.25, psym=symcat(6,thick=10), $
;      symsize=4, color=im_color('black')
;
;    djs_oplot, [8.85*scale], [ised[1].mass_50], psym=symcat(16,thick=5), $
;      symsize=4, color=im_color('dodger blue')
;    djs_oplot, [8.85*scale], [ised[1].mass_50], psym=symcat(9,thick=10), $
;      symsize=4, color=im_color('black')
;    djs_oplot, [8.85*scale], [ised[1].mass_50]+0.25, psym=symcat(9,thick=10), $
;      symsize=4, color=im_color('black')
;
;    djs_oplot, [60.0*scale], [ised[2].mass_50], psym=symcat(14,thick=5), $
;      symsize=4.5, color=im_color('dodger blue')
;    djs_oplot, [60.0*scale], [ised[2].mass_50], psym=symcat(4,thick=10), $
;      symsize=4.5, color=im_color('black')
;    djs_oplot, [60.0*scale], [ised[2].mass_50]+0.25, psym=symcat(4,thick=10), $
;      symsize=4.5, color=im_color('black')

;   djs_oploterr, 3.0*scale/2.0, ised[0].mass_50, yerr=ised[0].mass_err, $
;     psym=symcat(6,thick=5), symsize=3.5
;   djs_oploterr, 8.85*scale/2.0, ised[1].mass_50, yerr=ised[1].mass_err, $
;     psym=symcat(16,thick=5), symsize=3.5
;   djs_oploterr, 60.0*scale/2.0, ised[2].mass_50, yerr=ised[2].mass_err, $
;     psym=symcat(14,thick=5), symsize=3.5

    im_plotconfig, psfile=psfile, /psclose, /pdf

return
end
    
