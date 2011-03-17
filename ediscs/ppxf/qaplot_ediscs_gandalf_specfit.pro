pro render_gandalf_lineplot, wave, flux_noc, linefit, specdata, $
  noerase=noerase, xrange=xrange, scale=scale, pos=pos, $
  linelabel=linelabel

    get_element, wave, xrange, xx
    max_norej = max(linefit[xx[0]:xx[1]])
;   max_norej = max(flux_noc[xx[0]:xx[1]])
    max_rej = im_max(flux_noc[xx[0]:xx[1]],sigrej=3.0)
    sigma_rej = djsig(flux_noc[xx[0]:xx[1]],sigrej=3.0)

    yrange = [(max_rej*(-0.12))<(-2.5*sigma_rej),$
      (max_norej*1.25)>(3.0*sigma_rej)]
;   yrange = (im_max(flux_noc[xx[0]:xx[1]],sigrej=3.0)>$
;     max(linefit[xx[0]:xx[1]]))*[-0.12,1.2]
    
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=xrange, $
      yrange=scale*yrange, position=pos, ytitle='', noerase=noerase, $
      xtickinterval=25
    djs_oplot, wave, scale*flux_noc, ps=10, color='grey', thick=3
    djs_oplot, wave, scale*linefit, ps=10, color='red', thick=2.5
    oplot_gandalf_linewave, specdata
    im_legend, linelabel, /left, /top, box=0, charsize=1.2

return
end

pro oplot_gandalf_linewave, specdata
; overplot the wavelength of the line, adjusting for the possible
; differences in the emission/absorption-line redshifts
    for jj = 0, n_elements(specdata.linename)-1 do begin
       line = strupcase(strtrim(specdata.linename[jj],2))
       linewave = specdata.(tag_indx(specdata,line+'_WAVE'))
       linez = (specdata.(tag_indx(specdata,line+'_LINEZ')))[0]
       if (linez eq 0.0) then factor = 1.0 else $
         factor = (1.0+linez)/(1.0+specdata.zabs)
       djs_oplot, factor*linewave*[1,1], !y.crange, line=5, color='blue'
    endfor
    djs_oplot, !x.crange, [0,0], line=0, color='red'
return
end

pro qaplot_gandalf, specdata, specfit, pos1
; make the actual QAplot - see below for the main routine
    
; define some convenient internal variables
    galaxy = strtrim(specdata.galaxy,2)
    
    good = where(specfit.wave gt 0.0)
    wave = exp(specfit.wave[good])
    flux = specfit.flux[good]
    linefit = specfit.linefit[good]
    continuum = specfit.continuum[good]
    smooth_continuum = specfit.smooth_continuum[good]
    flux_noc = flux - continuum - smooth_continuum

    xtitle1 = 'Rest Wavelength (\AA)'
    power = ceil(abs(alog10(abs(median(flux)))))
    scale = 10.0^power
    ytitle1 = 'Flux (10^{-'+string(power,format='(I0)')+'} '+flam_units()+')'
    
; ---------------------------------------------------------------------------    
; page1 - continuum plot

; first third    
    xrange = [min(wave)-10.0,(max(wave)-min(wave))/3.0+min(wave)]
;   xrange = [min(wave),mean(wave)]
    get_element, wave, xrange, xx
    stats = im_stats(flux[xx[0]:xx[1]],sigrej=5.0)
    yrange = [-0.12,1.2]*stats.maxrej

    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=xrange, $
      yrange=scale*yrange, position=pos1[*,0], ytitle='', $
      xtitle='', title=galaxy
    djs_oplot, wave, scale*flux, ps=10, color='grey'
    djs_oplot, wave, scale*continuum, ps=10, color='red', thick=2.5
    djs_oplot, wave, scale*(continuum+smooth_continuum), ps=10, $
      color='blue', thick=2.5
    djs_oplot, wave, scale*smooth_continuum, ps=10, color='dark green', thick=2.5

    im_legend, /right, /top, box=0, margin=0, charsize=1.2, $
      ['\chi^2_{\nu}='+strtrim(string(specdata.continuum_chi2,format='(F12.2)'),2),$
      'E(B-V)='+strtrim(string(specdata.continuum_ebv,format='(F12.3)'),2),$;+'\pm'+$
;     strtrim(string(specdata.continuum_ebv_err,format='(F12.3)'),2),$
      '\sigma='+strtrim(string(specdata.vdisp,format='(I0)'),2)+' km s^{-1}',$
;     '\sigma='+strtrim(string(specdata.vdisp,format='(F12.1)'),2)+'\pm'+$
;     strtrim(string(specdata.vdisp_err,format='(F12.1)'),2)+' km s^{-1}']
      'S/N='+strtrim(string(specdata.continuum_snr,format='(F12.1)'),2)]

    im_legend, /left, /top, box=0, margin=0, charsize=1.2, $
      ['Z='+strtrim(string(specdata.continuum_Z,format='(F5.3)'),2),$
      'z_{NED}='+strtrim(string(specdata.z,format='(F12.5)'),2),$
      'z_{abs}='+strtrim(string(specdata.zabs,format='(F12.5)'),2)]
;     'z_{Balmer}='+strtrim(string(specdata.zline_balmer,format='(F12.5)'),2),$
;     'z_{forbid}='+strtrim(string(specdata.zline_forbidden,format='(F12.5)'),2)]

; second third
    xrange = [(max(wave)-min(wave))/3.0+min(wave),(max(wave)-min(wave))*2.0/3.0+min(wave)]
    get_element, wave, xrange, xx
    stats = im_stats(flux[xx[0]:xx[1]],sigrej=5.0)
    yrange = [-0.12,1.2]*stats.maxrej

    djs_plot, [0], [0], /nodata, /noerase, /xsty, /ysty, $
      xrange=xrange, yrange=scale*yrange, position=pos1[*,1], $
      ytitle=ytitle1, xtitle=''
    djs_oplot, wave, scale*flux, ps=10, color='grey'
    djs_oplot, wave, scale*continuum, ps=10, color='red', thick=2.5
    djs_oplot, wave, scale*(continuum+smooth_continuum), ps=10, $
      color='blue', thick=2.5
    djs_oplot, wave, scale*smooth_continuum, ps=10, color='dark green', thick=2.5
    
; third third    
    xrange = [(max(wave)-min(wave))*2.0/3.0+min(wave),max(wave)]
;   xrange = [mean(wave),max(wave)]
    get_element, wave, xrange, xx
    stats = im_stats(flux[xx[0]:xx[1]],sigrej=5.0)
    yrange = [-0.12,1.2]*stats.maxrej

    djs_plot, [0], [0], /nodata, /noerase, /xsty, /ysty, $
      xrange=xrange, yrange=scale*yrange, position=pos1[*,2], $
      ytitle='', xtitle=xtitle1
    djs_oplot, wave, scale*flux, ps=10, color='grey'
    djs_oplot, wave, scale*continuum, ps=10, color='red', thick=2.5
    djs_oplot, wave, scale*(continuum+smooth_continuum), ps=10, $
      color='blue', thick=2.5
    djs_oplot, wave, scale*smooth_continuum, ps=10, color='dark green', thick=2.5
    
; ---------------------------------------------------------------------------    
; page2 - emission-line plot

    npanel = 10
    im_plotconfig, 23, pos2, xspace=0.5, yspace=0.3D*[1,1,1,1], $
       width=[3.25,3.25], height=1.55D*[1,1,1,1,1]
    linelabel = ['[Ne V]','[O II]','[Ne III]','H5','H\delta','H\gamma',$
      'H\beta','[O III]','H\alpha+[N II]','[S II]']
    xrange = fltarr(2,npanel)
    xrange[*,0] = 3426.0
    xrange[*,1] = 3727.0
    xrange[*,2] = 3869.69
    xrange[*,3] = 3889.05
    xrange[*,4] = 4101.0
    xrange[*,5] = 4340.0
    xrange[*,6] = 4861.0
    xrange[*,7] = [4959.0,5007.0]
    xrange[*,8] = [6548.0,6584.0]
    xrange[*,9] = [6716.0,6731.0]
    xrange = xrange+rebin(30.0*[-1,1],2,npanel)

;    npanel = 8 ; 10
;    im_plotconfig, 19, pos2, xspace=0.5, yspace=[0.5,0.5,0.5], $
;       width=[3.25,3.25], height=2.0*[1,1,1,1]
;    linelabel = ['[Ne V]','[O II]','H\delta','H\gamma',$
;      'H\beta','[O III]','H\alpha+[N II]','[S II]']
;    xrange = fltarr(2,npanel)
;    xrange[*,0] = 3426.0
;    xrange[*,1] = 3727.0
;    xrange[*,2] = 4101.0
;    xrange[*,3] = 4340.0
;    xrange[*,4] = 4861.0
;    xrange[*,5] = [4959.0,5007.0]
;    xrange[*,6] = [6548.0,6584.0]
;    xrange[*,7] = [6716.0,6731.0]
;    xrange = xrange+rebin(30.0*[-1,1],2,npanel)

    for ii = 0, npanel-1 do begin
       render_gandalf_lineplot, wave, flux_noc, linefit, $
         specdata, xrange=xrange[*,ii], scale=scale, $
         linelabel=linelabel[ii], pos=pos2[*,ii], noerase=(ii gt 0)
    endfor
    
return
end
    
pro qaplot_ediscs_gandalf_specfit, specdata, specfit=specfit, $
  psfile=psfile, solar=solar
; jm10apr28ucsd - build a QAplot from the EDISCS_GANDALF_SPECFIT output

    version = ediscs_version(/ppxf_specfit)    
    specfitpath = ediscs_path(/ppxf)

    if keyword_set(solar) then suffix = 'solar'
    if (n_elements(psfile) eq 0) then psfile = specfitpath+$
      'qaplot_ediscs_'+suffix+'_'+version+'.ps'

    if (n_elements(specdata) eq 0L) then $
      specdata = read_ediscs_gandalf(solar=solar)
    specfit = read_ediscs_gandalf_specfit(specdata,solar=solar)

; build the plot    
    im_plotconfig, 4, pos1, psfile=psfile, $
      yspace=[0.5,0.5], charsize=1.2, ymargin=[0.6,1.0], $
      thick=2.0, height=replicate(2.8,3)

    index = lindgen(n_elements(specdata))
    if (index[0] ne -1) then begin
       for jj = 0, n_elements(index)-1 do begin
          print, format='("Object ",I0,"/",I0,A10,$)', jj+1, $
            n_elements(index), string(13b)
          qaplot_gandalf, specdata[index[jj]], $
            specfit[index[jj]], pos1
       endfor
    endif
       
    im_plotconfig, psfile=psfile, /gzip, /psclose

return
end    
