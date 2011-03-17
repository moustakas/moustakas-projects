pro ppxf_oplot, specdata, specfit, pos1, pos2
; driver routine to plot a spectrum and the best fit; see below for
; the main routine
    
; define some convenient internal variables
    galaxy = strtrim(specdata.nice_galaxy,2)
    
    good = where(specfit.wave gt 0.0)
    wave = exp(specfit.wave[good])
    flux = specfit.flux[good]
    ferr = specfit.ferr[good]
    linefit = specfit.linefit[good]
    continuum = specfit.continuum[good]
    smooth_continuum = specfit.smooth_continuum[good]
    flux_noc = flux-continuum-smooth_continuum

    xtitle1 = 'Rest Wavelength (\AA)'
    power = ceil(abs(alog10(abs(median(flux)))))
    scale = 10.0^power
    ytitle1 = 'Flux (10^{-'+string(power,format='(I0)')+'} '+flam_units()+')'
    
; -------------------------
; top panel -- full spectral range
    xrange = minmax(wave)
    get_element, wave, xrange, xx
    stats = im_stats(flux[xx[0]:xx[1]],sigrej=2.8)
    yrange = [-0.05,1.0]*stats.maxrej

    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=xrange, $
      yrange=scale*yrange, position=pos1, ytitle=ytitle, $
      xtitle=xtitle
    legend, galaxy, /left, /top, box=0, margin=0
    djs_oplot, wave, scale*flux, ps=10, color=fsc_color('grey',101), thick=4
    djs_oplot, wave, scale*continuum, ps=10, thick=3, color='black'
    djs_oplot, wave, scale*(flux-continuum), ps=10, $
      color=fsc_color('grey',101), thick=4
    djs_oplot, wave, scale*linefit, $
      ps=10, color=fsc_color('firebrick',101), thick=4
;   djs_oplot, wave, scale*(smooth_continuum+linefit), $
;     ps=10, color=fsc_color('firebrick',101), thick=4
    
    xyouts, pos1[0]-0.05, pos1[1], textoidl(ytitle1), align=0.5, $
      orientation=90, /norm
    
; bottom panels -- zoom around interesting sets of emission lines
    xrange = fltarr(2,3)
    xrange[*,0] = 3727.0+50*[-1,1]
    xrange[*,1] = [4861,5007]+50*[-1,1]
    xrange[*,2] = [6548,6731]+[-30,30]
    label = ['[O II]','H\beta+[O III]','H\alpha+[N II]+[S II]']
    xxtick = [50,100,100]
    for ii = 0, 2 do begin
       get_element, wave, xrange[*,ii], xx
       stats = im_stats(flux[xx[0]:xx[1]],sigrej=3.0)
;      stats = im_stats(flux_noc[xx[0]:xx[1]],sigrej=3.0)
       yrange = [0.9*stats.min,1.05*stats.max]

       djs_plot, [0], [0], /nodata, /noerase, /xsty, /ysty, xrange=xrange[*,ii], $
         yrange=scale*yrange, position=pos2[*,ii], ytitle=ytitle, $
         xtitle=xtitle1, xtickinterval=xxtick[ii]
       im_legend, label[ii], left=(ii lt 2), right=(ii eq 2), /top, box=0, margin=0, charsize=1.4

; plot the data       
       djs_oplot, wave, scale*flux, ps=10, color=fsc_color('grey',101), thick=6
       djs_oplot, wave, scale*(continuum+linefit), ps=10, thick=8, color=fsc_color('firebrick',101)
       djs_oplot, wave, scale*continuum, ps=10, thick=6, color='black'
;      djs_oplot, wave, scale*(continuum+smooth_continuum+linefit), ps=10, thick=6, color='black'
;      djs_oplot, wave, scale*flux_noc, ps=10, color='grey', thick=6
;      djs_oplot, wave, scale*linefit, ps=10, thick=6, color='black'
    endfor

return
end

pro plotsings_ppxf_example, ps=ps
; jm10jul02ucsd - plot some example PPXF fitting results for the paper 

    pspath = sings_path(/papers)+'log12oh/FIG_LOG12OH/'
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'
    
    specdata = read_sings_gandalf(/drift20)
    specdata = specdata[0]
    specfit = read_sings_gandalf_specfit(specdata,/drift20)
    nobj = n_elements(specdata)
    
; set up the plot
    psfile = pspath+'ppxf_example'+suffix
    im_plotconfig, 0, pos1, psfile=psfile, height=3.5, $
      width=10D, ymargin=[0.6,5.0], xmargin=[1.1,0.4], $
      charsize=1.6
    im_plotconfig, 13, pos2, xspace=[0.5D,0.5D], width=3.0D*[1,1,1], $
      height=3.5D, ymargin=[5.0,1.1], xmargin=[1.1,0.4]

    for ii = 0, nobj-1 do ppxf_oplot, specdata[ii], specfit[ii], pos1, pos2

    im_plotconfig, psfile=psfile, /psclose

stop
    
return
end
    
