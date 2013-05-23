pro mgfe_oplot, specdata, specfit, pos1, pos2
; driver routine to plot a spectrum and the best fit; see below for
; the main routine
    
; define some convenient internal variables
    galaxy = strtrim(specdata.pass,2)+'/'+strtrim(specdata.aper,2)
    
    good = where(specfit.wave gt 0.0)
    wave = exp(specfit.wave[good])
    flux = specfit.flux[good]
    ferr = specfit.ferr[good]
    linefit = specfit.linefit[good]
    continuum = specfit.continuum[good]
    smooth_continuum = specfit.smooth_continuum[good]

    xtitle1 = 'Rest Wavelength (\AA)'
    power = ceil(abs(alog10(abs(median(flux)))))
    scale = 10.0^power
    ytitle1 = 'Flux (10^{-'+string(power,format='(I0)')+'} '+flam_units()+')'
    
; -------------------------
; top panel -- full spectral range
    xrange = minmax(wave)
    get_element, wave, xrange, xx
    stats = im_stats(flux[xx[0]:xx[1]],sigrej=5.0)
    yrange = [-0.05,1.0]*stats.maxrej

    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=xrange, $
      yrange=scale*yrange, position=pos1, ytitle=ytitle, $
      xtitle=xtitle, title=galaxy
    djs_oplot, wave, scale*flux, ps=10, color='grey', thick=4
    djs_oplot, wave, scale*continuum, ps=10, thick=4, color='black'
    djs_oplot, wave, scale*(flux-continuum), ps=10, $
      color=fsc_color('medium grey',101), thick=6
    djs_oplot, wave, scale*(smooth_continuum+linefit), $
      ps=10, color='red', thick=4
    
;   djs_oplot, wave, scale*(continuum+smooth_continuum), ps=10, $
;     color='blue', thick=2.5
;   djs_oplot, wave, scale*linefit, ps=10, thick=6, color='red'
;   djs_oplot, wave, scale*smooth_continuum, ps=10, color='dark green', thick=6

    xyouts, pos1[0]-0.05, pos1[1], textoidl(ytitle1), align=0.5, $
      orientation=90, /norm
    
    im_legend, /right, /bottom, box=0, margin=0, charsize=1.2, $
      ['\chi^2_{\nu}='+strtrim(string(specdata.continuum_chi2,format='(F12.2)'),2),$
      'E(B-V)='+strtrim(string(specdata.continuum_ebv,format='(F12.3)'),2),$;+'\pm'+$
      '\sigma='+strtrim(string(specdata.vdisp,format='(I0)'),2)+' km/s',$
;     '\sigma='+strtrim(string(specdata.vdisp,format='(I0)'),2)+' km s^{-1}',$
      'S/N='+strtrim(string(specdata.continuum_snr,format='(F12.1)'),2)], /clear

;    if (specdata.zline_balmer le 0.0) then dz_balmer = 'N/A' else $
;      dz_balmer = strtrim(string(im_light()*(specdata.zline_balmer-specdata.zabs),format='(I0)'),2)
;    if (specdata.zline_forbidden le 0.0) then dz_forbidden = 'N/A' else $
;      dz_forbidden = strtrim(string(im_light()*(specdata.zline_forbidden-specdata.zabs),format='(I0)'),2)
        
    im_legend, /left, /top, box=0, margin=0, charsize=1.2, $
      ['Z='+strtrim(string(specdata.continuum_Z,format='(F5.3)'),2),$
      'z_{AGES}='+strtrim(string(specdata.z,format='(F12.5)'),2),$
      'z_{abs}='+strtrim(string(specdata.zabs,format='(F12.5)'),2)]
;     '\delta'+'z_{Balmer}='+dz_balmer,$
;     '\delta'+'z_{forbidden}='+dz_forbidden], /clear

; -------------------------
; bottom panels -- zoom around Mgb, Fe5270, and Fe5335
    xrange = [[5130,5220],[5225,5330],[5295,5370]]
    label = ['Mgb','Fe5270','Fe5335']
    for ii = 0, 2 do begin
       get_element, wave, xrange[*,ii], xx
       stats = im_stats(flux[xx[0]:xx[1]],sigrej=3.0)
       yrange = [0.95*stats.minrej,stats.maxrej]
;      yrange = [-0.04,1.05]*stats.maxrej

;      if (ii eq 1) then xtitle = xtitle1 else delvarx, xtitle
       djs_plot, [0], [0], /nodata, /noerase, /xsty, /ysty, xrange=xrange[*,ii], $
         yrange=scale*yrange, position=pos2[*,ii], ytitle=ytitle, $
         xtitle=xtitle1, xtickinterval=100
; overplot the continuum windows
       case ii of
          0: begin ; Mgb
             wblue = [5142.625,5161.375]
             wred  = [5191.375,5206.375]
             windex = [5160.125,5192.625]
          end
          1: begin ; Fe5270
             wblue  = [5233.150,5248.150]
             wred   = [5285.650,5318.150]
             windex = [5245.650,5285.650]
          end
          2: begin ; Fe5335
             wblue  = [5304.625,5315.875] 
             wred   = [5353.375,5363.375]
             windex = [5312.125,5352.125]
          end
          else: 
       endcase
; blue continuum       
       polyfill, [wblue[0],wblue[1],wblue[1],wblue[0]], $ 
         [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
         /line_fill, orientation=135, color=djs_icolor('blue'), $
         spacing=0.05, linestyle=1, noclip=0
; red continuum       
       polyfill, [wred[0],wred[1],wred[1],wred[0]], $ 
         [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
         /line_fill, orientation=135, color=djs_icolor('red'), $
         spacing=0.05, linestyle=1, noclip=0
; index continuum
       if (windex[0] ne -1) then begin
          polyfill, [windex[0],windex[1],windex[1],windex[0]], $ 
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=135, color=fsc_color('tan',101), $
            spacing=0.05, linestyle=1, noclip=0
       endif       
; now plot the data       
       djs_oplot, wave, scale*flux, ps=10, color='grey', thick=6
       djs_oplot, wave, scale*(continuum+smooth_continuum+linefit), $
         ps=10, color='red', thick=6
       djs_oplot, wave, scale*continuum, ps=10, thick=6, color='black'
;      djs_oplot, wave, scale*(flux-linefit), ps=10, $
;        color=fsc_color('medium grey',101), thick=6
       im_legend, label[ii], /left, /top, box=0, margin=0, charsize=1.7
    endfor

;; -------------------------
;; page 2 - zoom into the emission lines
;    linename = ['oii_3727','neiii_3869','h5','h_epsilon',$
;      'h_delta','h_gamma','h_beta','oiii_5007']
;    nicename = ['[O II]','[Ne III]',$
;      'H5','H\epsilon','H\delta','H\gamma','H\beta','[O III]']
;    nline = n_elements(linename)
;    im_plotconfig, 19, linepos, xspace=0.5, yspace=[0.5,0.5,0.5], $
;      width=[3.25,3.25], height=2.0*[1,1,1,1]
;    xrange = fltarr(2,nline)
;    xrange[*,0] = 3727.0
;    xrange[*,1] = 3869.69
;    xrange[*,2] = 3889.05
;    xrange[*,3] = 3970.07
;    xrange[*,4] = 4101.0
;    xrange[*,5] = 4340.0
;    xrange[*,6] = 4861.0
;    xrange[*,7] = [4959.0,5007.0]
;    xrange = xrange+rebin(30.0*[-1,1],2,nline)
;    flux_noc = flux-continuum-smooth_continuum
;    for ii = 0, nline-1 do begin
;       get_element, wave, xrange[*,ii], xx
;       max_norej = max(linefit[xx[0]:xx[1]])
;       max_rej = im_max(flux_noc[xx[0]:xx[1]],sigrej=3.0)
;       sigma_rej = djsig(flux_noc[xx[0]:xx[1]],sigrej=3.0)
;    
;       yrange = [(max_rej*(-0.12))<(-2.5*sigma_rej),$
;         (max_norej*1.25)>(3.0*sigma_rej)]
;    
;       djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=xrange[*,ii], $
;         yrange=scale*yrange, position=linepos[*,ii], ytitle='', noerase=(ii gt 0), $
;         xtickinterval=25
;       djs_oplot, wave, scale*flux_noc, ps=10, color='grey', thick=3
;       djs_oplot, wave, scale*linefit, ps=10, color='red', thick=2.5
;
;; overplot the position of the line and make a legend       
;       line = strupcase(linename[ii])
;       linewave = specdata.(tag_indx(specdata,line+'_WAVE'))
;       linez = (specdata.(tag_indx(specdata,line+'_LINEZ')))[0]
;       linesigma = (specdata.(tag_indx(specdata,line+'_SIGMA')))[0]
;       linesnr = (specdata.(tag_indx(specdata,line)))[0]/$
;         (specdata.(tag_indx(specdata,line)))[1]
;       if (linez eq 0.0) then factor = 1.0 else $
;         factor = (1.0+linez)/(1.0+specdata.zabs)
;       djs_oplot, factor*linewave*[1,1], !y.crange, line=5, color='blue'
;       djs_oplot, !x.crange, [0,0], line=0, color='red'
;; legend info
;       if (linez le 0) then linezlabel = '' else $
;         linezlabel = '\delta'+'z='+strtrim(string(im_light()*(linez-specdata.zabs),format='(I0)'),2)
;       if (linesigma le 0) then linesigmalabel = '' else $
;         linesigmalabel = '\sigma='+strtrim(string(linesigma,format='(I0)'),2)
;       if (linesnr le 0) then linesnrlabel = '' else $
;         linesnrlabel = 'S/N='+strtrim(string(linesnr,format='(I0)'),2) 
;       im_legend, [nicename[ii],linezlabel,linesigmalabel,linesnrlabel], $
;         /left, /top, box=0, charsize=1.0, margin=0, /clear
;    endfor 
    
return
end

pro qaplots_ages_mgfe
; jm13apr14siena - make some spectroscopic QAplots

    mgfepath = getenv('IM_PROJECTS_DIR')+'/mgfe/'
    ppxf = mrdfits(mgfepath+'ages_mgfe_ppxf.fits.gz',1)
    specfit = mrdfits(mgfepath+'ages_mgfe_specfit.fits.gz',1)
;   kcorr = mrdfits(mgfepath+'ages_mgfe_kcorr.fits.gz',1)
    ngal = n_elements(ppxf)
    
; read Pascale's table of non emission-line galaxies
;   readcol, mgfepath+'greg/table3.dat', galaxy, format='A', $
;     comment='#', /silent

; cluster    
;   match, strtrim(cluster.galaxy,2), galaxy, m1, m2
;   qaplot_mgfe_gandalf_specfit, cluster[m1], psfile=mgfepath+'qaplot_cluster_noemission.ps'
    
;; field
;    field = read_mgfe_sfh_sample(/field)
;    match, strtrim(field.galaxy,2), galaxy, m1, m2
;    qaplot_mgfe_gandalf_specfit, field[m1], psfile=mgfepath+'qaplot_field_noemission.ps'

; --------------------------------------------------
; QAplot of the spectral fits zoomed into the Mgb, Fe5270, and Fe5335
; indices; group many of the passes together to keep the file size
; manageable 
    allpass = fix(ppxf.pass/10.0)
    uindx = uniq(allpass,sort(allpass))
    pass = allpass[uindx]

    for ii = 0, n_elements(pass)-1 do begin
       psfile = mgfepath+'qaplots/qa_specfit_'+string(ii+1,format='(I2.2)')+'.ps'
;      psfile = mgfepath+'qaplots/qa_specfit_'+string(pass[ii],format='(I2.2)')+'.ps'
       indx = where(allpass eq pass[ii],nindx)

; reverse-sort by S/N
       snr = ppxf[indx].continuum_snr
       indx = indx[reverse(sort(snr))]

       im_plotconfig, 13, pos2, xspace=[0.5D,0.5D], width=3.0D*[1,1,1], $
         height=3.5D, ymargin=[5.0,1.1], xmargin=[1.1,0.4]
       im_plotconfig, 0, pos1, psfile=psfile, height=3.5, $
         width=10D, ymargin=[0.6,5.0], xmargin=[1.1,0.4], $
         charsize=1.5

       djs_plot, [0], [0], /nodata, xsty=4, ysty=4, $
         position=pos2[*,1], title='AGES MgFe Sample: Pass '+$
         string(min(ppxf[indx].pass),format='(I0)')+'-'+$
         string(max(ppxf[indx].pass),format='(I0)')
       spec = specfit[indx]
;      for jj = 0, nindx-1 do mgfe_oplot, ppxf[indx[jj]], $
       for jj = 0, 10-1 do mgfe_oplot, ppxf[indx[jj]], $
         spec[jj], pos1, pos2
       im_plotconfig, psfile=psfile, /psclose;, /pdf
    endfor

stop    
    
return
end
    
