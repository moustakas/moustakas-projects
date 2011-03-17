pro bootes_irac_apercor
; jm10oct08ucsd - test empirical aperture correction schemes for IRAC 

    common com_apercor, iband, ch1

    if (n_elements(iband) eq 0L) then begin
       bootespath = getenv('RESEARCHPATH')+'/data/bootes/2010b/'
       iband = mrdfits(bootespath+'bootes_I.fits.gz',1)
       ch1 = mrdfits(bootespath+'bootes_ch1.fits.gz',1)
       keep = where((iband.mag_auto gt 0) and (iband.mag_auto lt 90.0))
       iband = iband[keep]
       ch1 = ch1[keep]
    endif

; point-source aperture correction    
    readcol, ages_path(/psfs)+'I1_psf_flux.dat', $
      diam, frac, format='F,F', /silent
    psfdiff = 2.5*alog10(interpol(frac,diam,4.0)/interpol(frac,diam,10.0))

; find isolated sources    
    iso = where((iband.mag_auto lt 22.0) and (iband.class_star lt 0.5) and $
      (ch1.segflags_aper_10 eq 0) and (iband.segflags_aper_10 eq 0))
    iso_stars = where((ch1.segflags_aper_10 eq 0) and $
      (iband.segflags_aper_10 eq 0) and (iband.class_star gt 0.95) and (iband.mag_auto lt 18.5))
    diff = ch1[iso].mag_aper_10-ch1[iso].mag_aper_04
    auto = iband[iso].mag_auto
    size = sqrt(iband[iso].a_world*iband[iso].b_world)*3600.0 ; [arcsec]

    diff_stars = ch1[iso_stars].mag_aper_10-ch1[iso_stars].mag_aper_04
    auto_stars = iband[iso_stars].mag_auto
    size_stars = sqrt(iband[iso_stars].a_world*iband[iso_stars].b_world)*3600.0 ; [arcsec]

    psfile = ages_path(/qaplots)+'bootes_irac_apercor.ps'
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.6, $
      xmargin=[1.2,0.3]
; vs mag_auto
    hogg_scatterplot, auto, diff, xsty=1, ysty=1, position=pos[*,0], $
      xrange=[15.0,22.5], yrange=[-1.3,0], ytitle='[3.6] m(10")-m(4") (mag)', $
      levels=[0.1,0.5,0.75,0.9], $
      outliers=1, /internal, xnpix=20, ynpix=20, $
      xtitle=textoidl('I_{auto} (Vega mag)')
    djs_oplot, auto_stars, diff_stars, psym=6, sym=0.4;, color='orange'
;   djs_oplot, auto_contam, diff_contam, psym=6, sym=0.1, color='orange'
    djs_oplot, !x.crange, psfdiff*[1,1], line=0, thick=6, color='blue'
; vs size
    hogg_scatterplot, alog10(size>0.1), diff, xsty=1, ysty=1, /noerase, position=pos[*,1], $
      xrange=alog10([0.15,3.0]), yrange=[-1.3,0], ytitle='', ytickname=replicate(' ',10), $
      levels=[0.1,0.5,0.75,0.9], outliers=1, /internal, xnpix=20, ynpix=20, $
      xtitle=textoidl('log (size) (arcsec)')
    djs_oplot, alog10(size_stars>0.1), diff_stars, psym=6, sym=0.4;, color='orange'
;   djs_oplot, alog10(size_contam>0.1), diff_contam, psym=6, sym=0.1, color='orange'
    djs_oplot, !x.crange, psfdiff*[1,1], line=0, thick=6, color='blue'
    im_plotconfig, psfile=psfile, /psclose, /gzip

stop    

return
end
    
