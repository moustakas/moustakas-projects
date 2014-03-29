pro render_deep2_gandalf_lineplot, wave, flux_noc, linefit, specdata, $
  noerase=noerase, xrange=xrange, scale=scale, pos=pos, $
  linelabel=linelabel, xtickinterval=xtickinterval

    get_element, wave, xrange, xx
    max_norej = max(linefit[xx[0]:xx[1]])
;   max_norej = max(flux_noc[xx[0]:xx[1]])
    max_rej = im_max(flux_noc[xx[0]:xx[1]],sigrej=3.0)
    min_rej = im_min(flux_noc[xx[0]:xx[1]],sigrej=3.0)
    sigma_rej = djsig(flux_noc[xx[0]:xx[1]],sigrej=3.0)

    yrange = [(max_rej*(-0.12))<(-2.5*sigma_rej),$
      (max_norej*1.25)>(3.0*sigma_rej)]
    if max_norej eq 0.0 then yrange =[0<min_rej,1>max_rej]
    if xx[1] eq xx[0] then begin
       ytickname = replicate(' ',10)
       yrange = [-0.05,1.0]
    endif
    
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=xrange, $
      yrange=scale*yrange, position=pos, ytitle='', noerase=noerase, $
      xtickinterval=xtickinterval, ytickname=ytickname
    djs_oplot, wave, scale*flux_noc, ps=10, color='grey', thick=3
    djs_oplot, wave, scale*linefit, ps=10, color='red', thick=4.0
    oplot_deep2_gandalf_linewave, specdata
    im_legend, linelabel, /left, /top, box=0, charsize=1.2

return
end

pro oplot_deep2_gandalf_linewave, specdata
; overplot the wavelength of the line, adjusting for the possible
; differences in the emission/absorption-line redshifts
    for jj = 0, n_elements(specdata.linename)-1 do begin
       line = strupcase(strtrim(specdata.linename[jj],2))
       linewave = specdata.(tag_indx(specdata,line+'_WAVE'))
       linez = (specdata.(tag_indx(specdata,line+'_LINEZ')))[0]
       if (linez eq 0.0) then factor = 1.0 else $
         factor = (1.0+linez)/(1.0+specdata.z)
       djs_oplot, factor*linewave*[1,1], !y.crange, line=5, color='blue'
    endfor
    djs_oplot, !x.crange, [0,0], line=0, color='red'
return
end

pro qaplot_deep2_gandalf, specdata, specfit, pos1, pos2, $
  linear=linear
; make the actual QAplot - see below for the main routine
    
    light = 2.99792458D5 ; speed of light [km/s]

; define some convenient internal variables
    galaxy = string(specdata.mask,format='(I0)')+'/'+$
      string(specdata.objno,format='(I0)')
    
    good = where(specfit.wave gt 0.0)
    if keyword_set(linear) then wave = specfit.wave[good] else $
      wave = exp(specfit.wave[good])
    flux = specfit.flux[good]
    linefit = specfit.linefit[good]
    continuum = specfit.continuum[good]
;   smooth_continuum = specfit.smooth_continuum[good]
    flux_noc = flux - continuum; - smooth_continuum

    xtitle1 = 'Rest Wavelength (\AA)'
    scale = 1.0
;   power = ceil(alog10(abs(median(flux))))
;   scale = 10.0^(-power)
;   ytitle1 = 'Flux (10^{-'+string(power,format='(I0)')+'} '+flam_units()+')'
    ytitle1 = 'Flux (ADU)'
    
; ---------------------------------------------------------------------------    
; page1 - continuum plot

; first third    
    xrange = [min(wave),(max(wave)-min(wave))/3.0+min(wave)]
;   xrange = [min(wave),mean(wave)]
    get_element, wave, xrange, xx
    stats = im_stats(flux[xx[0]:xx[1]],sigrej=5.0)
    yrange = [-0.12,1.2]*stats.maxrej

    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=xrange, $
      yrange=scale*yrange, position=pos1[*,0], ytitle='', $
      xtitle='', title=galaxy
    djs_oplot, wave, scale*flux, ps=10, color='grey'
    djs_oplot, wave, scale*continuum, ps=10, color='red', thick=2.5
;   djs_oplot, wave, scale*(continuum+smooth_continuum), ps=10, $
;     color='blue', thick=2.5

    im_legend, /right, /top, box=0, margin=0, charsize=1.2, $
      [$
;     '\chi^2_{\nu}='+strtrim(string(specdata.continuum_chi2,format='(F12.2)'),2),$
;     'E(B-V)='+strtrim(string(specdata.continuum_ebv,format='(F12.3)'),2),$;+'\pm'+$
;     strtrim(string(specdata.continuum_ebv_err,format='(F12.3)'),2),$
;     '\sigma='+strtrim(string(specdata.vdisp,format='(F12.1)'),2)+' km s^{-1}',$
;     '\sigma='+strtrim(string(specdata.vdisp,format='(F12.1)'),2)+'\pm'+$
;     strtrim(string(specdata.vdisp_err,format='(F12.1)'),2)+' km s^{-1}']
      'S/N='+strtrim(string(specdata.continuum_snr,format='(F12.1)'),2)]

    im_legend, /left, /top, box=0, margin=0, charsize=1.2, $
      ['z='+strtrim(string(specdata.z,format='(F12.5)'),2),$
      'z_{Balmer}='+strtrim(string(specdata.zline_balmer,format='(F12.5)'),2),$
      'z_{forbid}='+strtrim(string(specdata.zline_forbidden,format='(F12.5)'),2)]

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
;   djs_oplot, wave, scale*(continuum+smooth_continuum), ps=10, $
;     color='blue', thick=2.5
    
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
;   djs_oplot, wave, scale*(continuum+smooth_continuum), ps=10, $
;     color='blue', thick=2.5
    
; ---------------------------------------------------------------------------    
; page2 - emission-line plot
    npanel = 8
    xrange = fltarr(2,npanel)
    linelabel = ['Mg II','[Ne V]','[O II]','[Ne III]','H\gamma',$
      'H\beta','[O III]','H\alpha+[N II]']
    xrange[*,0] = 2800.0
    xrange[*,1] = 3426.0
    xrange[*,2] = 3727.0
    xrange[*,3] = 3869.0
    xrange[*,4] = 4340.0
    xrange[*,5] = 4861.0
    xrange[*,6] = [4959.0,5007.0]
    xrange[*,7] = [6548.0,6584.0]
    xtickint = [10.0,10.0,10.0,10.0,10.0,10.0,25.0,25.0]
    xrange = xrange+rebin(15.0*[-1,1],2,npanel)

    for ii = 0, npanel-1 do begin
       render_deep2_gandalf_lineplot, wave, flux_noc, linefit, $
         specdata, xrange=xrange[*,ii], scale=scale, $
         linelabel=linelabel[ii], pos=pos2[*,ii], noerase=(ii gt 0), $
         xtickinterval=xtickint[ii]
    endfor

return
end
    
pro qaplot_deep2_gandalf_specfit_dr4, specdata, specfit, $
  linear=linear, psfile=psfile
; jm13dec20siena - build a QAplot from the READ_DEEP2_GANDALF_SPECFIT
;   output; NOTE! /LINEAR must be set if it was set when calling
;   READ_DEEP2_GANDALF_SPECFIT(), otherwise DOOM!! 
;
; alldata = read_deep2_gandalf(/ppxf)
; data = alldata[where(alldata.z gt 0.5 and alldata.z lt 0.5003)]
; specfit = read_deep2_gandalf_specfit(data,/linear)
; qaplot_deep2_gandalf_specfit, data, specfit, /linear

    nobj = n_elements(specdata)
    if (nobj eq 0L) then begin
       doc_library, 'qaplot_deep2_gandalf_specfit'
       return
    endif

    fixoii = specdata[0].fixoii ; assume all are the same
    if (n_elements(specfit) eq 0L) then begin
       fixoii = specdata[uniq(specdata.fixoii,sort(specdata.fixoii))].fixoii
       if n_elements(fixoii) eq 1 then begin
          specfit = read_deep2_gandalf_specfit_dr4(specdata,linear=linear,fixoii=fixoii[0])
       endif else begin
          zero = where(specdata.fixoii eq 0,comp=one)
          specfit0 = read_deep2_gandalf_specfit_dr4(specdata[zero],linear=linear,fixoii=0)
          specfit1 = read_deep2_gandalf_specfit_dr4(specdata[one],linear=linear,fixoii=1)
          specfit = im_empty_structure(specfit0[0],ncopies=nobj)
          specfit[zero] = temporary(specfit0)
          specfit[one] = temporary(specfit1)
       endelse
    endif

    if (n_elements(psfile) eq 0) then psfile = $
      'qaplot_deep2_gandalf_specfit.ps'
    
    im_plotconfig, 19, pos2, xspace=0.5, yspace=[0.5,0.5,0.5], $
       width=[3.25,3.25], height=[2.0,2.0,2.0,2.0], ymargin=[0.6,0.9]
    im_plotconfig, 4, pos1, psfile=psfile, yspace=[0.5,0.5], $
      charsize=1.2, ymargin=[0.6,1.0], thick=2.0, $
      height=replicate(2.8,3)
       
    for jj = 0, nobj-1 do begin
       print, format='("Object ",I0,"/",I0,A10,$)', jj, nobj-1, $
         string(13b)
       qaplot_deep2_gandalf, specdata[jj], specfit[jj], $
         pos1, pos2, linear=linear
    endfor

    im_plotconfig, psfile=psfile, /psclose, /pdf
;   spawn, 'rsync -auv '+psfile+'.gz ~', /sh

return
end    
