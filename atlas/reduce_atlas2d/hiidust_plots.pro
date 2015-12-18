pro hiidust_plots
; jm06mar27uofa

    optir = rsex('hiidust_optir.txt')
    hii = mrdfits('53821_hiidust_specdata_speclinefit_nodust.fits.gz',1,/silent)
    hii = hii[[6,7,8,9,10,11,12,13,15,16]]
    niceprint, hii.specfile, optir.galaxy, optir.id

    im_window, 0, xratio=0.6, /square
    plotsym, 0, 1.5, /fill
    djs_plot, [0], [0], /nodata, xthick=2.0, ythick=2.0, xsty=3, ysty=3, charsize=2.0, $
      charthick=2.0, xtitle='log (24/H\alpha)', ytitle='E(H\beta-H\alpha)', xrange=[0.8,2.6], $
      yrange=[0,0.9]
    oploterror, alog10(optir.flux24/optir.fluxha), hii.ehbha, hii.ehbha_err, ps=8, thick=2.0
    niceprint, alog10(optir.flux24/optir.fluxha), hii.ehbha, hii.ehbha_err
    
    

return
end
    
