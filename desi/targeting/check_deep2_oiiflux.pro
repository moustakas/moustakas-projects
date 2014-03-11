pro check_deep2_oiiflux
; jm14mar10siena - compare [OII] fluxes of objects in the DEEP2 Groth
; Strip (Field 1) to those measured from some special SDSS plates and
; from the VLT

    deep2path = getenv('IM_PROJECTS_DIR')+'/desi/deep2/'
    targpath = getenv('IM_PROJECTS_DIR')+'/desi/targeting/'

; match to BOSS
    deep2 = read_deep2_zcat()
    ppxf = read_deep2(/ppxf)
    kised = mrdfits(deep2path+'desi_deep2_fsps_v2.4_miles_'+$
      'chab_charlot_sfhgrid01_kcorr.z0.0.fits.gz',1)

    boss = mrdfits(targpath+'sdss-plug-grothstrip.fits.gz',1)
    line = mrdfits(targpath+'sdss-zline-grothstrip.fits.gz',1)
    line = reform(line,31,n_elements(line)/31)
    zans = mrdfits(targpath+'sdss-zans-grothstrip.fits.gz',1)

    spherematch, deep2.ra, deep2.dec, boss.ra, boss.dec, 3D/3600.0, m1, m2
    deep2oii = deep2_get_oiiflux(ppxf[m1],kised[m1],oiierr=deep2oiierr,$
      oiilimit=deep2oiilimit)
    these = where(deep2oii ne -2.0 and zans[m2].zwarning eq 0 and $
      abs(3E5*(deep2[m1].zbest-zans[m2].z)) lt 2000.0 and $
      (line[6,m2].linearea_err gt 0 or line[7,m2].linearea_err gt 0),ngal)
    splog, 'BOSS-DEEP2 sample ', ngal

    zobj = deep2[m1[these]].zbest
    deep2oii = deep2oii[these]
    deep2oiierr = deep2oiierr[these]
    deep2oiilimit = deep2oiilimit[these]

    bossoii = 1D-17*reform(line[6,m2[these]].linearea+line[7,m2[these]].linearea)
    bossoiierr = 1D-17*reform(sqrt(line[6,m2[these]].linearea_err^2+line[7,m2[these]].linearea_err^2))

; make the plot
    psfile = targpath+'qa_check_oiiflux.ps'
    im_plotconfig, 6, pos, psfile=psfile, height=[3.0,2.5], width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xtitle='Redshift', ytitle='log_{10} ([OII]_{DEEP2}/[OII]_{BOSS})', $
      xrange=[0.7,1.2], yrange=2.5*[-1,1];, yminor=3
    djs_oplot, !x.crange, [0,0], line=5

    good = where(deep2oiilimit eq 0 and bossoii gt 0)
    ratio = alog10(bossoii[good]/deep2oii[good])
    ratioerr = im_compute_error(bossoii[good],bossoiierr[good],$
      deep2oii[good],deep2oiierr[good],/log)
    oploterror, zobj[good], ratio, ratioerr, psym=symcat(16)
    
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
stop    

return
end



