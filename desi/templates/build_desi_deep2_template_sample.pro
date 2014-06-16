function get_ugriz, cat, bri=bri, errbri=brierr, ugrizerr=ugrizerr, $
  wise=wise, errwise=errwise
; compute dust-corrected ugriz magnitudes
    deep2_to_maggies, cat, mm, ii, /unwise, ratag='ra_deep', dectag='dec_deep'
    mag = maggies2mag(mm,ivarmaggies=ii,magerr=magerr)
    bri = mag[0:2,*]
    brierr = magerr[0:2,*]
    ugriz = mag[3:7,*]
    ugrizerr = magerr[3:7,*]
    wise = mag[8:9,*]
    errwise = magerr[8:9,*]
return, ugriz
end

pro build_desi_deep2_template_sample, out
; jm13dec18siena - build the sample of DEEP2 galaxies we will use to
;   construct DESI templates; DESI_DEEP2_ISEDFIT needs to have been
;   run first!
; jm14mar11siena - major update

    version = 'v1.1'
    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/'+version+'/'

    oiisnrcut = 3.0
    oiifluxcut = 8D-17
    
    zcat = read_deep2_zcat(photo=photo)
    kised = mrdfits(templatepath+'desi_deep2_fsps_v2.4_miles_'+$
      'chab_charlot_sfhgrid02_kcorr.z0.0.fits.gz',1)
    photo = deep2_get_ugriz(photo,/unwise)
    
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

; a tiny number of these have CFLUX_3727=0 because of insufficient
; bands of photometry; apply the S/N cut on the *amplitude* of [OII]
; 3729, the stronger of the two lines
    oii = deep2_get_oiiflux(line,cflux_3727_rest=kised.cflux_3727) ; [erg/s/cm2/A]
    
    ewoiisnr = (oii.oii_3727_2_ew[0]/oii.oii_3727_2_ew[1])*(oii.oii_3727_2[1] ne -2)
    oiisnr = (oii.oii_3727_2_amp[0]/oii.oii_3727_2_amp[1])*(oii.oii_3727_2_amp[1] ne -2)
;   oiisnr = oii.oii_3727[0]/(oii.oii_3727[1]+(oii.oii_3727[1] eq 0))*(oii.oii_3727[1] ne 0)

    scale = 1D18
    djs_plot, scale*oii.oii_3727[0], oiisnr, psym=3, $
      xsty=1, ysty=1, /ylog, /xlog, xr=scale*[1D-18,1D-14], yr=[0.1,100]
    djs_oplot, 10^!x.crange, oiisnrcut*[1,1], color='red'
    djs_oplot, scale*oiifluxcut*[1,1], 10^!y.crange, color='red'

;   ww = where(oiisnr lt oiisnrcut and scale*oii.oii_3727[0] gt 500)
;   qaplot_deep2_gandalf_specfit_dr4, line[ww], psfile='~/junk.ps'
    
; apply the final sample cuts; note that the (soft) redshift cut is
; explicitly needed because of DESI_DEEP2_ISEDFIT
    these = where($
      zcat.zbest gt 0.75 and zcat.zbest lt 1.45 and $
;     zcat.zbest gt 0.1 and zcat.zbest lt 2.0 and $
      (total(photo.ugriz gt 0,1) gt 3) and $
      (total(photo.wise_err gt 0,1) ge 1) and $
      photo.ugriz[2] gt 18.5 and photo.ugriz[2] lt 24.0 and $
      oiisnr gt oiisnrcut and ewoiisnr gt 1.0,ngal)
    splog, 'Sample ', ngal

; build an output data structure with all the information we need
    out = im_struct_trimtags(zcat[these],select=['objno','mask','ra','dec','zbest'],$
      newtags=['objno','mask','ra','dec','z'])
    out = struct_addtags(struct_addtags(temporary(out),struct_trimtags(photo[these],$
      select=['ugriz*','bri*','wise*'])),replicate({cflux_3727_rest: 0.0},ngal))
    out.cflux_3727_rest = kised[these].cflux_3727

    out = struct_addtags(temporary(out),struct_trimtags(kised[these],$
      select=['isedfit_id','maggies','ivarmaggies','absmag']))
    out = struct_addtags(temporary(out),struct_trimtags(oii[these],$
      except='*limit*'))

    im_mwrfits, out, templatepath+'desi_deep2elg_template_sample_'+version+'.fits', /clobber
    
return
end
