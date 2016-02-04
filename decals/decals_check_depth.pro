pro decals_check_depth
; jm16jan29siena - check the depth in the annottated CCDs file

    common com_depth, jj, rr

    if n_elements(jj) eq 0L then jj = mrdfits(getenv('LEGACYPIPE_DIR')+$
      '/decals-ccds-annotated.fits',1)
    if n_elements(rr) eq 0L then rr = mrdfits(getenv('HOME')+$
      '/decals-obs-2015oct19.fits',1)
    ww = where(jj.tilepass eq 1 and jj.filter eq 'g' and jj.photometric eq 'T')

    psfile = 'arjun_vs_dustin.eps'
    im_plotconfig, 0, pos, height=5.0, psfile=psfile, xmargin=[1.3,0.2]
    djs_plot, jj[ww].seeing, rr.depth[ww]-jj[ww].psfdepth, xsty=3, ysty=3, $
      color=cgcolor('grey'), psym=cgsymcat(9), xtitle='FWHM Seeing (arcsec)', $
      ytitle='Depth (Arjun minus Dustin, AB mag)', yrange=[-1,0.7]
    im_plotconfig, psfile=psfile, /psclose, /png

    psfile = 'exptime_vs_seeing.eps'
    im_plotconfig, 0, pos, height=5.0, psfile=psfile, xmargin=[1.3,0.2]
    djs_plot, jj[ww].exptime, jj[ww].seeing, psym=cgsymcat(9), symsize=0.5, $
      ysty=3, xsty=3, xtitle='Exposure Time (s)', ytitle='FWHM Seeing (arcsec)', $
      position=pos, color=cgcolor('grey')
    im_legend, ['Pass 1','g-band'], /left, /top, box=0, margin=0
    djs_oplot, 100*[1,1], !y.crange, line=0
    im_plotconfig, psfile=psfile, /psclose, /png

    long = where(jj[ww].exptime gt 100)
    psfile = 'seeing_vs_depth.eps'
    im_plotconfig, 0, pos, height=5.0, psfile=psfile, xmargin=[1.3,0.2]
    djs_plot, jj[ww].seeing, jj[ww].psfdepth+0.37*0, psym=cgsymcat(9), symsize=0.5, $
      ysty=3, xsty=3, ytitle='Point-Source Depth (5\sigma, AB mag)', $
      xtitle='FWHM Seeing (arcsec)', position=pos, color=cgcolor('grey'), $
      yrange=[22.5,25.0]
    djs_oplot, jj[ww[long]].seeing, jj[ww[long]].psfdepth+0.37*0, psym=cgsymcat(9), $
      symsize=0.5, color=cgcolor('dodger blue')
    djs_oplot, !x.crange, 24.7*[1,1], line=0
    djs_oplot, 1.5*[1,1], !y.crange, line=0
    im_legend, ['Pass 1','g-band'], /right, /top, box=0, margin=0
    im_legend, ['Exptime>100 s'], /left, /bottom, box=0, margin=0, $
      psym=9, color='dodger blue'
    im_plotconfig, psfile=psfile, /psclose, /png


stop

return
end
