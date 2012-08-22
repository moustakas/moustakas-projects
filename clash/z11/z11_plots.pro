pro z11_plots
; jm12aug14siena

    path = clash_path(/z11)
    loz = mrdfits(path+'isedfit/z11_loz_fsps_chab_calzetti_sfhgrid01.fits.gz',1)
    hiz = mrdfits(path+'isedfit/z11_hiz_fsps_chab_calzetti_sfhgrid01.fits.gz',1)

    zobj = [loz.zobj,hiz.zobj]
    chi2 = [loz.chi2,hiz.chi2]
    chi2min = min(chi2)
    prob = exp(-0.5*((chi2-chi2min)))
    prob = prob/total(prob,/double)

    srt = sort(zobj)
    zobj = zobj[srt]
    prob = prob[srt]
    
    psfile = path+'isedfit_pofz.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, xmargin=[1.4,0.4], $
      width=6.6

    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=[0.5,12], $ ; /ylog, $
      yrange=[1E-6,0.02], position=pos, xtitle='Redshift', ytitle='Probability'
    djs_oplot, zobj, prob, thick=8, line=0, color='blue'
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
return
end
    
