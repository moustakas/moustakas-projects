pro ediscs_sfh_oplot, specdata, specfit, pos1, pos2
; driver routine to plot a spectrum and the best fit; see below for
; the main routine
    
; define some convenient internal variables
    galaxy = strtrim(specdata.galaxy,2)
    
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
    
    im_legend, /right, /top, box=0, margin=0, charsize=1.2, $
      ['\chi^2_{\nu}='+strtrim(string(specdata.continuum_chi2,format='(F12.2)'),2),$
      'E(B-V)='+strtrim(string(specdata.continuum_ebv,format='(F12.3)'),2),$;+'\pm'+$
      '\sigma='+strtrim(string(specdata.vdisp,format='(I0)'),2),$
;     '\sigma='+strtrim(string(specdata.vdisp,format='(I0)'),2)+' km s^{-1}',$
      'S/N='+strtrim(string(specdata.continuum_snr,format='(F12.1)'),2)], /clear

    if (specdata.zline_balmer le 0.0) then dz_balmer = 'N/A' else $
      dz_balmer = strtrim(string(im_light()*(specdata.zline_balmer-specdata.zabs),format='(I0)'),2)
    if (specdata.zline_forbidden le 0.0) then dz_forbidden = 'N/A' else $
      dz_forbidden = strtrim(string(im_light()*(specdata.zline_forbidden-specdata.zabs),format='(I0)'),2)
        
    im_legend, /left, /top, box=0, margin=0, charsize=1.2, $
      ['Z='+strtrim(string(specdata.continuum_Z,format='(F5.3)'),2),$
      'z_{EDisCS}='+strtrim(string(specdata.z,format='(F12.5)'),2),$
      'z_{abs}='+strtrim(string(specdata.zabs,format='(F12.5)'),2),$
      '\delta'+'z_{Balmer}='+dz_balmer,$
      '\delta'+'z_{forbidden}='+dz_forbidden], /clear

; -------------------------
; bottom panels -- zoom around 4000-A break, Hd_A, and Hg_A
;   xrange = [[3840,4110],[4000,4200],[4250,4450]]
    xrange = [[3825,4125],[4025,4175],[4275,4425]]
    label = ['D_{n}(4000)','H\delta_{A}','H\gamma_{A}']
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
          0: begin ; D(4000)
             wblue = [3850,3950]
             wred  = [4000,4100]
             windex = -1
          end
          1: begin ; Hd_A
             wblue  = [4041.6,4079.75]
             wred   = [4128.5,4161.0]
             windex = [4083.5,4122.25]
          end
          2: begin ; Hg_A
             wblue  = [4283.5,4319.75]
             wred   = [4367.25,4419.75]
             windex = [4319.75,4363.5]
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

; -------------------------
; page 2 - zoom into the emission lines
    linename = ['oii_3727','neiii_3869','h5','h_epsilon',$
      'h_delta','h_gamma','h_beta','oiii_5007']
    nicename = ['[O II]','[Ne III]',$
      'H5','H\epsilon','H\delta','H\gamma','H\beta','[O III]']
    nline = n_elements(linename)
    im_plotconfig, 19, linepos, xspace=0.5, yspace=[0.5,0.5,0.5], $
      width=[3.25,3.25], height=2.0*[1,1,1,1]
    xrange = fltarr(2,nline)
    xrange[*,0] = 3727.0
    xrange[*,1] = 3869.69
    xrange[*,2] = 3889.05
    xrange[*,3] = 3970.07
    xrange[*,4] = 4101.0
    xrange[*,5] = 4340.0
    xrange[*,6] = 4861.0
    xrange[*,7] = [4959.0,5007.0]
    xrange = xrange+rebin(30.0*[-1,1],2,nline)
    flux_noc = flux-continuum-smooth_continuum
    for ii = 0, nline-1 do begin
       get_element, wave, xrange[*,ii], xx
       max_norej = max(linefit[xx[0]:xx[1]])
       max_rej = im_max(flux_noc[xx[0]:xx[1]],sigrej=3.0)
       sigma_rej = djsig(flux_noc[xx[0]:xx[1]],sigrej=3.0)
    
       yrange = [(max_rej*(-0.12))<(-2.5*sigma_rej),$
         (max_norej*1.25)>(3.0*sigma_rej)]
    
       djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=xrange[*,ii], $
         yrange=scale*yrange, position=linepos[*,ii], ytitle='', noerase=(ii gt 0), $
         xtickinterval=25
       djs_oplot, wave, scale*flux_noc, ps=10, color='grey', thick=3
       djs_oplot, wave, scale*linefit, ps=10, color='red', thick=2.5

; overplot the position of the line and make a legend       
       line = strupcase(linename[ii])
       linewave = specdata.(tag_indx(specdata,line+'_WAVE'))
       linez = (specdata.(tag_indx(specdata,line+'_LINEZ')))[0]
       linesigma = (specdata.(tag_indx(specdata,line+'_SIGMA')))[0]
       linesnr = (specdata.(tag_indx(specdata,line)))[0]/$
         (specdata.(tag_indx(specdata,line)))[1]
       if (linez eq 0.0) then factor = 1.0 else $
         factor = (1.0+linez)/(1.0+specdata.zabs)
       djs_oplot, factor*linewave*[1,1], !y.crange, line=5, color='blue'
       djs_oplot, !x.crange, [0,0], line=0, color='red'
; legend info
       if (linez le 0) then linezlabel = '' else $
         linezlabel = '\delta'+'z='+strtrim(string(im_light()*(linez-specdata.zabs),format='(I0)'),2)
       if (linesigma le 0) then linesigmalabel = '' else $
         linesigmalabel = '\sigma='+strtrim(string(linesigma,format='(I0)'),2)
       if (linesnr le 0) then linesnrlabel = '' else $
         linesnrlabel = 'S/N='+strtrim(string(linesnr,format='(I0)'),2) 
       im_legend, [nicename[ii],linezlabel,linesigmalabel,linesnrlabel], $
         /left, /top, box=0, charsize=1.0, margin=0, /clear
    endfor 
    
return
end

pro qaplots_ediscs_sfh
; jm10jun01ucsd - make some spectroscopic QAplots

    sfhpath = ediscs_path(/projects)+'sfh/'
    cluster = read_ediscs_sfh_sample(/cluster,/all)
    field = read_ediscs_sfh_sample(/field,/all)
    ppxf = read_ediscs(/ppxf)
    kcorr = read_ediscs(/kcorr)
    info = read_ediscs(/spec1d)

; read Pascale's table of non emission-line galaxies
;   readcol, sfhpath+'greg/table3.dat', galaxy, format='A', $
;     comment='#', /silent

; cluster    
;   match, strtrim(cluster.galaxy,2), galaxy, m1, m2
;   qaplot_ediscs_gandalf_specfit, cluster[m1], psfile=sfhpath+'qaplot_cluster_noemission.ps'
    
;; field
;    field = read_ediscs_sfh_sample(/field)
;    match, strtrim(field.galaxy,2), galaxy, m1, m2
;    qaplot_ediscs_gandalf_specfit, field[m1], psfile=sfhpath+'qaplot_field_noemission.ps'

; --------------------------------------------------
; QAplot of the spectral fits zoomed into the high-order Balmer lines;
; make separate plots for each targeted cluster, for convenience

    allcl = strtrim(cluster.cluster,2)
    uindx = uniq(allcl,sort(allcl))
    cl = allcl[uindx]

;   for ii = 6, n_elements(cl)-1 do begin
    for ii = 0, n_elements(cl)-1 do begin
       psfile = sfhpath+'qaplots/qa_specfit_'+cl[ii]+'.ps'
       iscl = where(strtrim(cluster.cluster,2) eq cl[ii],niscl)
       isfield = where(strtrim(field.cluster,2) eq cl[ii],nisfield)
       splog, cl[ii], niscl, nisfield

; reverse-sort by S/N
       cl_snr = cluster[iscl].continuum_snr
       iscl = iscl[reverse(sort(cl_snr))]
       if (nisfield ne 0) then begin
          field_snr = field[isfield].continuum_snr
          isfield = isfield[reverse(sort(field_snr))]
       endif
;      iscl = iscl[1:3] & isfield = isfield[1:3]
;      niscl = 3 & nisfield = 3

       im_plotconfig, 13, pos2, xspace=[0.5D,0.5D], width=3.0D*[1,1,1], $
         height=3.5D, ymargin=[5.0,1.1], xmargin=[1.1,0.4]
       im_plotconfig, 0, pos1, psfile=psfile, height=3.5, $
         width=10D, ymargin=[0.6,5.0], xmargin=[1.1,0.4], $
         charsize=1.5
; cluster galaxies
       djs_plot, [0], [0], /nodata, xsty=4, ysty=4, $
         position=pos2[*,1], title=cl[ii]+': Cluster Sample'
       cspec = read_ediscs_gandalf_specfit(cluster[iscl])
       for jj = 0, niscl-1 do ediscs_sfh_oplot, cluster[iscl[jj]], $
         cspec[jj], pos1, pos2
; field galaxies
       djs_plot, [0], [0], /nodata, xsty=4, ysty=4, $
         position=pos2[*,1], title=cl[ii]+': Field Sample'
       if (nisfield ne 0) then begin
          fspec = read_ediscs_gandalf_specfit(field[isfield])
          for jj = 0, nisfield-1 do ediscs_sfh_oplot, field[isfield[jj]], $
            fspec[jj], pos1, pos2
       endif
       im_plotconfig, psfile=psfile, /psclose, /gzip
    endfor

stop
    
; --------------------------------------------------
; S/N vs I-band magnitude and M_V for each cluster 

    psfile = sfhpath+'qaplots/qaplot_snr_vs_mag.ps'
    im_plotconfig, 4, pos, psfile=psfile, yspace=[0.0,1.0], $
      height=[2.5,2.5,2.5], charsize=1.7, ymargin=[0.8,1.1]

    xrange = [0.2,100]
    yrange1 = [17.5,25]
    yrange2 = [-16,-24.5]

    xtitle = 'Continuum S/N (pixel^{-1})'
    ytitle1 = 'I_{tot} (AB mag)'
    ytitle2 = 'M_{V} (AB mag)'

; full set of plots for each cluster
    allcl = strtrim(info.cluster_fullname,2)
    uindx = uniq(allcl,sort(allcl))
    cl = allcl[uindx]
    target = strtrim(info[uindx].cluster,2) ; targeted cluster
    good = where(strcompress(cl,/remove) ne '',ncl)
    cl = cl[good]
    target = target[good]

    for ii = -1, ncl-1 do begin
       if (ii eq -1) then begin
          iscl = where((kcorr.maggies[3] gt 0.0) and $
            (strmatch(info.memberflag,'*1*')),niscl)
          isfield = -1
          title = 'All clusters'
          symsize1 = 0.7
       endif else begin
          iscl = where((allcl eq cl[ii]) and (kcorr.maggies[3] gt 0.0) and $
            (strmatch(info.memberflag,'*1*')),niscl)
          isfield = where((strtrim(info.cluster,2) eq target[ii]) and $
            (kcorr.maggies[3] gt 0.0) and $
            (strmatch(info.memberflag,'*1*') eq 0) and $
            (abs(info.z-info[iscl[0]].cluster_z) lt 0.1))
          title = cl[ii]
          symsize1 = 1.4
       endelse

; vs I-mag    
       djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
         xtickname=replicate(' ',10), xtitle='', ytitle=ytitle1, $
         xrange=xrange, yrange=yrange1, /xlog, title=title
       if (isfield[0] ne -1) then djs_oplot, ppxf[isfield].continuum_snr, $
         -2.5*alog10(kcorr[isfield].maggies[3]), psym=symcat(15), $
         color=fsc_color('dodger blue',101), symsize=symsize1
       djs_oplot, ppxf[iscl].continuum_snr, -2.5*alog10(kcorr[iscl].maggies[3]), $
         psym=symcat(16), color=fsc_color('firebrick',101), symsize=symsize1
       djs_oplot, 3.0*[1,1], !y.crange, line=0, thick=4
       djs_oplot, 10^!x.crange, 23*[1,1], line=5, thick=6

; vs M_V
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
         xtitle=xtitle, ytitle=ytitle2, xrange=xrange, yrange=yrange2, /xlog
       if (isfield[0] ne -1) then djs_oplot, ppxf[isfield].continuum_snr, $
         kcorr[isfield].ubvrijhk_absmag_00[2], psym=symcat(15), $
         color=fsc_color('dodger blue',101), symsize=symsize1
       djs_oplot, ppxf[iscl].continuum_snr, kcorr[iscl].ubvrijhk_absmag_00[2], $
         psym=symcat(16), color=fsc_color('firebrick',101), symsize=symsize1
;      djs_oplot, allfield.continuum_snr, allfield.ubvrijhk_absmag_00[2], $
;        psym=symcat(15), color=fsc_color('dodger blue',101), symsize=symsize1
       djs_oplot, 3.0*[1,1], !y.crange, line=0, thick=4
       djs_oplot, 10^!x.crange, -19.0*[1,1], line=5, thick=6

if ii eq 5 then stop
; plot absolute magnitude vs I-band magnitude
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,2], xsty=1, ysty=1, $
         xtitle=ytitle1, ytitle=ytitle2, xrange=yrange1, yrange=yrange2
       if (isfield[0] ne -1) then djs_oplot, -2.5*alog10(kcorr[isfield].maggies[3]), $
         kcorr[isfield].ubvrijhk_absmag_00[2], psym=symcat(15), $
         color=fsc_color('dodger blue',101), symsize=symsize1
       djs_oplot, -2.5*alog10(kcorr[iscl].maggies[3]), kcorr[iscl].ubvrijhk_absmag_00[2], $
         psym=symcat(16), color=fsc_color('firebrick',101), symsize=symsize1
;      djs_oplot, -2.5*alog10(allfield.maggies[3]), allfield.ubvrijhk_absmag_00[2], $
;        psym=symcat(15), color=fsc_color('dodger blue',101), symsize=symsize1
       djs_oplot, !x.crange, -19.0*[1,1], line=0, thick=4
       djs_oplot, 23*[1,1], !y.crange, line=0, thick=4
    endfor
       
    im_plotconfig, psfile=psfile, /psclose, /gzip

    
stop    
    
return
end
    
