pro ppxf_oplot, specdata, specfit, ancillary, pos1, pos2
; driver routine to plot a spectrum and the best fit; see below for
; the main routine
    
; define some convenient internal variables
    galaxy = hogg_iau_name(ancillary.ra,ancillary.dec,'AGES')

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
    xrange = [3600.0,5100.0]
;   xrange = minmax(wave)
    get_element, wave, xrange, xx
    stats = im_stats(flux[xx[0]:xx[1]],sigrej=2.8)
    yrange = [-0.13,1.5]*stats.maxrej

    nsmooth = 5
    
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=xrange, $
      yrange=scale*yrange, position=pos1, ytitle=ytitle, $
      xtitle=xtitle
    djs_oplot, wave, scale*smooth(flux,nsmooth), ps=10, $
      color=im_color('grey'), thick=6
    djs_oplot, wave, scale*smooth(continuum,nsmooth), ps=10, $
      thick=3, color='black'

    djs_oplot, wave, scale*smooth((flux-continuum-smooth_continuum),nsmooth), $
      ps=10, color=im_color('dark grey'), thick=4
    djs_oplot, wave, scale*smooth(linefit,nsmooth), $
      ps=10, color=im_color('firebrick'), thick=4
;   djs_oplot, wave, scale*(smooth_continuum+linefit), $
;     ps=10, color=im_color('firebrick',101), thick=4
    
    xyouts, pos1[0]-0.05, pos1[1], textoidl(ytitle1), align=0.5, $
      orientation=90, /norm

;   legend, galaxy, /left, /top, box=0, margin=0, $
;     position=[xrange[0]+250,scale*yrange[1]-0.15], /data

; bottom panels -- zoom around interesting sets of emission lines
    xrange = fltarr(2,3)
    xrange[*,0] = 3727.0+30*[-1,1]
    xrange[*,1] = [4861,5007]+30*[-1,1]
    xrange[*,2] = [6548,6584]+[-30,30]
    label = ['[O II]','H\beta+[O III]','H\alpha+[N II]']
    xxtick = [25,50,100]
    for ii = 0, 1 do begin
;   for ii = 0, 2 do begin
       get_element, wave, xrange[*,ii], xx
       stats = im_stats(flux[xx[0]:xx[1]],sigrej=3.0)
;      stats = im_stats(flux_noc[xx[0]:xx[1]],sigrej=3.0)
;      yrange = [0.9*stats.min,1.05*stats.max]
       yrange = [0.7,5.5]/scale

       djs_plot, [0], [0], /nodata, /noerase, /xsty, /ysty, xrange=xrange[*,ii], $
         yrange=scale*yrange, position=pos2[*,ii], ytitle=ytitle, $
         xtitle=xtitle1, xtickinterval=xxtick[ii]
       im_legend, label[ii], left=(ii lt 2), right=(ii eq 2), /top, box=0, margin=0, charsize=1.4

; plot the data       
       djs_oplot, wave, scale*flux, ps=10, color=im_color('grey'), thick=8
       djs_oplot, wave, scale*(continuum+linefit), ps=10, thick=6, color=im_color('firebrick')
       djs_oplot, wave, scale*continuum, ps=10, thick=6, color='black'
;      djs_oplot, wave, scale*(continuum+smooth_continuum+linefit), ps=10, thick=6, color='black'
;      djs_oplot, wave, scale*flux_noc, ps=10, color='grey', thick=6
;      djs_oplot, wave, scale*linefit, ps=10, thick=6, color='black'
    endfor

return
end

pro mzplot_ppxf_example, pdf=pdf
; jm12aug22siena - example PPXF fit for the paper

    mzpath = mz_path()
    qapath = mzpath+'qaplots/'
    if keyword_set(pdf) then begin
       pspath = qapath
       suffix = '.pdf'
    endif else begin
       pspath = mz_path(/paper)
       suffix = '.eps'
    endelse

    ancillary = read_mz_sample(/mzhii_ancillary)
    specdata = read_mz_sample(/mzhii_ispec)
;   ww = where(ancillary.pass eq 304 and ancillary.aper eq 114,nobj)
    ww = where(ancillary.pass eq 605 and ancillary.aper eq 208,nobj)
    ancillary = ancillary[ww]
    specdata = specdata[ww]
    
    specfit = read_ages_gandalf_specfit(specdata)
    
    psfile = pspath+'ppxf_example'+suffix
    im_plotconfig, 8, pos1, psfile=psfile, height=3.5, $
      width=9.5D, ymargin=[0.6,5.0], xmargin=[1.2,0.3], $
      charsize=1.8
    im_plotconfig, 1, pos2, xspace=0.5D, width=[3.5D,5.0D], $
      height=3.5D, ymargin=[5.0,1.1], xmargin=[1.2,0.3]
    for ii = 0, nobj-1 do ppxf_oplot, specdata[ii], specfit[ii], ancillary[ii], pos1, pos2
    im_plotconfig, psfile=psfile, /psclose, pdf=pdf

return
end
    
