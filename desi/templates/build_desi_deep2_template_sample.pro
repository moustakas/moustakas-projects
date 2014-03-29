function get_ugriz, cat, bri=bri, errbri=brierr, ugrizerr=ugrizerr
; compute dust-corrected ugriz magnitudes
    deep2_to_maggies, cat, mm, ii
    mag = maggies2mag(mm,ivarmaggies=ii,magerr=magerr)
    bri = mag[0:2,*]
    brierr = magerr[0:2,*]
    ugriz = mag[3:7,*]
    ugrizerr = magerr[3:7,*]
return, ugriz
end

pro build_desi_deep2_template_sample, out
; jm13dec18siena - build the sample of DEEP2 galaxies we will use to
;   construct DESI templates
; jm14mar11siena - major update

    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/'

    snrcut = 3.0
    oiifluxcut = 8D-17
    
    zcat = read_deep2_zcat(photo=photo)
    kised = mrdfits(templatepath+'desi_deep2_fsps_v2.4_miles_'+$
      'chab_charlot_sfhgrid01_kcorr.z0.0.fits.gz',1)
    ugriz = get_ugriz(photo,bri=bri,errbri=brierr,ugrizerr=ugrizerr)
    
; use the fixed-[OII] catalog; we'll add some jitter in the
; doublet-ratio in BUILD_DESI_DEEP2_TEMPLATES
    line = read_deep2(/ppxf,/fixoii)
    
;; use the fixed-[OII] catalog only when one of the components of the
;; doublet has been dropped
;    line = read_deep2(/ppxf)
;    linefixoii = read_deep2(/ppxf,/fixoii)
;    these = where((line.oii_3727_1[1] gt 0 and line.oii_3727_2[1] eq -1) or $
;      (line.oii_3727_1[1] eq -1 and line.oii_3727_2[1] gt 0),nthese)
;    line[these] = linefixoii[these]

; [OII] continuum flux, CFLUX_3727
    cflux_3727_rest = kised.cflux_3727            ; [erg/s/cm2/A]
    cflux_3727_obs = cflux_3727_rest/(1.0+zcat.z) ; [erg/s/cm2/A]
    
; a tiny number of these have CFLUX_3727=0 because of insufficient
; bands of photometry; apply the S/N cut on the *amplitude* of [OII]
; 3729, the stronger of the two lines
    oii = deep2_get_oiiflux(line,cflux_3727_obs=cflux_3727_obs,$
      cflux_3727_rest=cflux_3727_rest)
    
    ewoiisnr = oii.oii_3727_2_ew[0]/oii.oii_3727_2_ew[1]
    oiisnr = oii.oii_3727_2_amp[0]/(oii.oii_3727_2_amp[1]+$
      (oii.oii_3727_2_amp[1] eq 0))*(oii.oii_3727_2_amp[1] ne 0)
    oiicontinuum = line.oii_3727_2_continuum[0]
;   oiisnr = oii.oii_3727[0]/(oii.oii_3727[1]+(oii.oii_3727[1] eq 0))*(oii.oii_3727[1] ne 0)
    
    scale = 1D18
    djs_plot, scale*oii.oii_3727[0], oiisnr, psym=3, $
      xsty=1, ysty=1, /ylog, /xlog, xr=scale*[1D-18,1D-14], yr=[0.1,100]
    djs_oplot, 10^!x.crange, snrcut*[1,1], color='red'
    djs_oplot, scale*oiifluxcut*[1,1], 10^!y.crange, color='red'

;   ww = where(oiisnr lt snrcut and scale*oii.oii_3727[0] gt 500)
;   qaplot_deep2_gandalf_specfit_dr4, line[ww], psfile='~/junk.ps'
    
; apply the final sample cuts; note that the redshift cut is
; explicitly needed because of DESI_DEEP2_ISEDFIT
    these = where($
      zcat.zbest gt 0.7 and zcat.zbest lt 1.5 and $
      total(ugriz gt 0,1) gt 3 and $
      ugriz[2,*] gt 18 and ugriz[2,*] lt 24.0 and $
      oiisnr gt snrcut and ewoiisnr gt 1.0,ngal)
    splog, 'Sample ', ngal

; build an output data structure with all the information we need
    out = im_struct_trimtags(zcat[these],select=['objno','ra','dec','zbest'],$
      newtags=['objno','ra','dec','z'])
    out = struct_addtags(temporary(out),replicate({bri: fltarr(3), $
      brierr: fltarr(3), ugriz: fltarr(5), ugrizerr: fltarr(5), $
      cflux_3727_obs: 0.0},ngal))
    out.bri = bri[*,these]
    out.brierr = brierr[*,these]
    out.ugriz = ugriz[*,these]
    out.ugrizerr = ugrizerr[*,these]
    out.cflux_3727_obs = cflux_3727_obs[these]

    out = struct_addtags(temporary(out),struct_trimtags(kised[these],$
      select=['isedfit_id','maggies','ivarmaggies','absmag']))
    out = struct_addtags(temporary(out),struct_trimtags(oii[these],$
      except='*limit*'))
    im_mwrfits, out, templatepath+'desi_deep2_template_sample.fits', /clobber
    
return
end
