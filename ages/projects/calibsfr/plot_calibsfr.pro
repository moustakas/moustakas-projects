pro plot_calibsfr, ps=ps
; jm10feb21ucsd - build the basic CALIBSFR plots

    calibsfrpath = ages_path(/projects)+'calibsfr/'
    paperpath = ages_path(/papers)+'calibsfr/'

    sample = read_calibsfr_sample()

    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'
    lumaxis = im_array(35.0,45.0,0.1)
    
    xrange = [39.5,43.8]
    yrange = [41.0,45.0]

; --------------------------------------------------
; L(24) vs L(Ha)_obs

    xtitle = 'log_{10} [L(H\alpha)_{obs}] (erg s^{-1})'
    ytitle = 'log_{10} [L(24 \mu'+'m)] (erg s^{-1})'
    
    coeff = robust_linefit(sample.lha,sample.l24_ce01,/bisect)

    psfile = paperpath+'l24_vs_lha_obs'+suffix
    im_plotconfig, 0, pos, psfile=psfile
    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle
    djs_oplot, sample.lha, sample.l24_ce01, $
      psym=symcat(6,thick=8), symsize=0.5
    djs_oplot, lumaxis, poly(lumaxis,coeff), line=0, $
      thick=6, color='red'
;   djs_oplot, lumaxis, lumaxis, line=0, thick=6.0
    djs_oplot, lumaxis, poly(lumaxis,[1.35,1.0]), line=2, $
      color='blue', thick=6.0 ; [Zhu+08]
    im_plotconfig, /psclose

; --------------------------------------------------
    xtitle = 'log_{10} [L(H\alpha)_{cor}] (erg s^{-1})'
    coeff = robust_linefit(sample.lha_cor,sample.l24_ce01,/bisect)

    psfile = paperpath+'l24_vs_lha_cor'+suffix
    im_plotconfig, 0, pos, psfile=psfile
    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle
    djs_oplot, sample.lha_cor, sample.l24_ce01, $
      psym=symcat(6,thick=8), symsize=0.5
    djs_oplot, lumaxis, poly(lumaxis,coeff), line=0, $
      thick=6, color='red'
    djs_oplot, lumaxis, poly(lumaxis,[1.35,1.0]), line=2, $
      color='blue', thick=6.0 ; [Zhu+08]
    im_plotconfig, /psclose

return
end
    
