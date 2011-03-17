pro qaplot_sfrm_sdss_galex
; jm10apr25ucsd - QAplots for the joint SDSS-GALEX sample

    sfrmpath = ages_path(/projects)+'sfrm/'
    paperpath = ages_path(/papers)+'sfrm/'

    parent = read_sfrm_sample(/sdss)

    insurvey = where(parent.galex_insurvey,ninsurvey)
    nuv = where(parent.galex_detect and parent.galex_insurvey,nnuv)
    fuv = where(parent.galex_detect and parent.galex_insurvey and $
      (parent.k_maggies[0] gt 0.0),nfuv)

    psfile = sfrmpath+'qaplots/sfrm_sdss_galex.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.2]

; --------------------------------------------------
; fraction of NUV/FUV detections    
    binsize = 0.2
    minmag = 12.0
    maxmag = 18.0

    rmag = -2.5*alog10(parent.k_maggies[4])
    im_plothist, rmag, bin=binsize, $
      min=minmag, max=maxmag, xall, yall, /noplot
    im_plothist, rmag[insurvey], bin=binsize, $
      min=minmag, max=maxmag, xinsurvey, yinsurvey, /noplot
    im_plothist, rmag[nuv], bin=binsize, $
      min=minmag, max=maxmag, xgalex, ynuv, /noplot
    im_plothist, rmag[fuv], bin=binsize, $
      min=minmag, max=maxmag, xfuv, yfuv, /noplot

    good = where(yall gt 0)
    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=1, ysty=1, yrange=[0,1], xrange=[minmag,maxmag], $
      ytitle='Fraction of Galaxies', xtitle='r (AB mag)'
    djs_oplot, xall[good], yinsurvey[good]/yall[good], $
      psym=10, color='red', line=0, thick=6
    djs_oplot, xall[good], ynuv[good]/yall[good], $
      psym=10, color='blue', line=5, thick=6
    djs_oplot, xall[good], yfuv[good]/yall[good], $
      psym=10, color='orange', line=3, thick=6
    im_legend, ['GALEX Footprint','NUV Detection','FUV Detection'], $
      /right, /top, box=0, line=[0,5,3], color=['red','blue','orange'], $
      thick=6, pspacing=2.0

; --------------------------------------------------
; redshift histograms
    binsize = 0.001
    zmin = 0.0
    zmax = 0.051

    zobj = parent.zdist
    im_plothist, zobj, bin=binsize, xall, yall, /noplot
    
    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=1, ysty=1, yrange=[0,max(yall)*1.05], $
      xrange=[zmin,zmax], ytitle='Number of Galaxies', $
      xtitle='Redshift'
    im_plothist, zobj, bin=binsize, /overplot, line=0, thick=6
    im_plothist, zobj[insurvey], bin=binsize, /overplot, $
      line=5, color='dark green', thick=6
    im_plothist, zobj[nuv], bin=binsize, /overplot, $
      line=3, color='red', thick=6
    im_plothist, zobj[fuv], bin=binsize, /overplot, $
      line=1, color='blue', thick=6
    im_legend, ['SDSS','GALEX Footprint','NUV Detection','FUV Detection'], $
      /left, /top, box=0, line=[0,5,3,1], thick=6, pspacing=2.0, $
      color=['default','dark green','red','blue','orange'], $
      charsize=1.6
    
    im_plotconfig, /psclose, psfile=psfile, /gzip
    
return
end
    
