pro talk_12nov_faculty_lecture
; jm12nov07siena - make a blackbody curve for my public talk at Siena 

    path = getenv('IM_RESEARCH_DIR')+'/talks/2012/12nov_faculty_lecture/'

    wave = range(0.5,6,5000,/log) ; millimeters
    bb = planck(wave*1D7,2.74) ; wave must be in Angstrom
    bb = bb/max(bb)

    dwave = range(0.5,6,40,/log)    ; millimeters
    dwave = dwave+dwave*randomn(seed,40)*0.001
    dbb = planck(dwave*1D7,2.74)  ; wave must be in Angstrom
    dbb = dbb/max(dbb)

;   ww = im_color('white',10)
    
    psfile = path+'cmb_blackbody.ps'
    im_plotconfig, 0, pos, height=4.5, psfile=psfile, charsize=3.0, $
      thick=10, charthick=6, xthick=10, ythick=10
    plot, [0], [0], /nodata, xsty=9, ysty=9, $
      ytickname=replicate(' ',10), yminor=-1, $
      xminor=-1, yticklen=1D-10, position=pos, $
      xrange=minmax(wave), yrange=minmax(bb), color=im_color('white',10)
    djs_oplot, wave, bb, color=im_color('cyan'), thick=10
    djs_oplot, dwave, dbb, psym=symcat(15), symsize=1.2, color=10
    im_plotconfig, psfile=psfile, /psclose, /pdf

return
end
