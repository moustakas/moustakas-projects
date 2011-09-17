pro talk_11sep_preston, keynote=keynote, noevol=noevol
; jm11sep04ucsd - miscellaneous plots for my 2011 Sep IAU talk in
;   Preston, UK

    mfpath = mf_path()
    mfspath = mf_path(/mfs)
    isedpath = mf_path(/isedfit)
    mzpath = mz_path()
    
    rootpath = getenv('IM_RESEARCH_DIR')+'/talks/2011/11sep_preston/'
    if keyword_set(keynote) then talkpath = rootpath+'keynote/' else $
      talkpath = rootpath+'figures/'

; --------------------------------------------------
; metallicity vs redshift in bins of mass - all three calibrations
    massbins = mz_massbins(nmassbins)
    psfile = talkpath+'z_vs_oh_bymass.ps'
    im_plotconfig, 5, pos, psfile=psfile, height=2.6*[1,1], $
      charsize=1.6, keynote=keynote

    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'

    calib = ['kk04','t04','m91']
    calibcolor = ['forest green','orange','red']
;   calibcolor = ['black','black','black']
    polycolor = ['forest green','orange','red']
    calibpsym = [6,5,9]
    sdss_calibpsym = [15,17,16]
;   calibpsym = [15,17,16]
    calibsymsize = [1.5,1.8,1.5]*1.1
    calibline = [0,5,3]
    ncalib = n_elements(calib)

    xrange = [0.0,0.72]
    yrange = [8.67,9.22]
    xtitle1 = 'Redshift'
    ytitle1 = mzplot_ohtitle(/mean)

    zaxis = range(xrange[0]+0.01,xrange[1]-0.01,50)
    
    for mm = 0, n_elements(massbins)-2 do begin
       if (mm lt 2) then begin
          xtitle = ''
          xtickname = replicate(' ',10)
       endif else begin
          xtitle = xtitle1
          delvarx, xtickname
       endelse
       if odd(mm) then begin
          ytitle = ''
          ytickname = replicate(' ',10)
       endif else begin
          ytitle = ytitle1
          delvarx, ytickname
       endelse
       
       plot, [0], [0], /nodata, noerase=(mm gt 0), ysty=1, xsty=1, xrange=xrange, $
         yrange=yrange, xtitle=xtitle, ytitle=ytitle, position=pos[*,mm], $
         xtickname=xtickname, ytickname=ytickname, yminor=2, xtickinterval=0.2, $
         color=im_color(keycolor)
       im_legend, massbins[mm].label, /right, /top, box=0, margin=0, $
         charsize=1.3, textcolor=im_color(keycolor)
       if (mm eq 0) then begin
          im_legend, strupcase(calib), /left, /bottom, box=0, margin=0, $
            charsize=1.2, pspacing=1.9, line=calibline, thick=8, $
            psym=-calibpsym, color=calibcolor, symsize=calibsymsize*0.8, $
            symthick=4.0, textcolor=im_color(keycolor)
       endif
; and now the data
       for ii = 0, ncalib-1 do begin
          mzevol = mrdfits(mzpath+'mzevol_'+calib[ii]+'.fits.gz',1,/silent)
          good = where(mzevol.ohmean_bymass[mm,*] gt -900.0)
; AGES          
          oploterror, mzevol.medz_bymass[mm,good], mzevol.ohmean_bymass[mm,good], $
            mzevol.ohmean_bymass_err[mm,good], $
            psym=symcat(calibpsym[ii],thick=5), symsize=calibsymsize[ii], $
            color=im_color(calibcolor[ii],101), errcolor=im_color(calibcolor[ii],101)
; SDSS
          oploterror, mzevol.sdss_medz_bymass[mm], mzevol.sdss_ohmean_bymass[mm], $
            mzevol.sdss_ohmean_bymass_err[mm], $
            psym=symcat(sdss_calibpsym[ii],thick=5), symsize=calibsymsize[ii], $
            color=im_color(calibcolor[ii],101), errcolor=im_color(calibcolor[ii],101)          
          
          djs_oplot, zaxis, poly(zaxis-mzevol.qz0,mzevol.coeffs_bymass[*,mm]), $
            line=calibline[ii], color=im_color(calibcolor[ii],101)
       endfor
    endfor

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

stop    
    
; ---------------------------------------------------------------------------
; make a legend for the NUV-r vs r-J color-color plot
    psfile = talkpath+'nuvmr_vs_rmj_legend.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote

    djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5, $
      xrange=[0,1], yrange=[0,1]
    label = ['-1.5<log (\psi/M)','-3.0<log (\psi/M)<-1.5',$
      '-3.0>log (\psi/M)']
    im_legend, label, /left, /top, box=0, charsize=1.5, $
      textcolor=['sky blue','spring green1','tomato']
    
    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

stop    
    
; ---------------------------------------------------------------------------
; <SFR/M> vs redshift and stellar mass
    zbins = mf_zbins()
    fit = mrdfits(mfpath+'sfrm_vs_redshift.fits.gz',1)

    psfile = talkpath+'sfrm_mean_vs_mass_redshift.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=5.8, charsize=1.7, keynote=keynote

    if keyword_set(keynote) then begin
       keycolor = 'white'
       color = ['plum','tan','turquoise','sky blue','orange','tomato']
    endif else begin
       keycolor = 'black'
       color = ['purple','forest green','midnight blue','sky blue','orange','dark red']
    endelse

    psym = [16,15,18,17,2,14]
    zrange = [0.12,1.05]
    sfrmrange = [-3.9,1.1]
;   sfrmrange = [-6.25,1.5]
    maxis = range(8,12,50)
    zaxis = range(0.01,1.5,50)
    ageaxis = getage(0.0)-getage(zaxis) ; [Gyr]

    allerrors = 1
    
    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      yrange=sfrmrange, xrange=zrange, xtitle='Redshift', $
      ytitle=textoidl('<log (\psi/M)> (Gyr^{-1})'), color=im_color(keycolor)
    masslabel = string(fit[0].mass,format='(F4.1)')
    im_legend, ['log (M/M'+sunsymbol()+')='+masslabel[0],$
      masslabel[1:n_elements(fit[0].mass)-1]], $
      psym=psym, color=color, /right, /bottom, box=0, margin=0, charsize=1.4, $
      textcolor=im_color(keycolor), position=[zrange[0]+0.31,sfrmrange[0]+0.2], /data
    im_legend, 't_{lookback}^{-1} [Gyr^{-1}]', /left, /top, box=0, $
      line=0, pspacing=1.6, thick=6, charsize=1.6, textcolor=im_color(keycolor), $
      color=im_color(keycolor)
    
    for mm = 0, n_elements(fit[0].mass)-1 do begin
       gd = where(fit.sf_sfrm[mm] gt -900,ngd)
       oploterror, fit[gd].sf_medz, fit[gd].sf_sfrm[mm], $
         fit[gd].sf_sfrm_cv_err[mm], psym=symcat(psym[mm]), $
         symsize=2.0, color=im_color(color[mm]), errcolor=im_color(color[mm])
       if allerrors then begin
; models
          for ii = 0, ngd-1 do begin
             xwid = fit[gd[ii]].sf_sigz*2
             yhei = fit[gd[ii]].sf_sfrm_models_max[mm]-fit[gd[ii]].sf_sfrm_models_min[mm]
             xcen = fit[gd[ii]].sf_medz
             ycen = yhei/2.0+fit[gd[ii]].sf_sfrm_models_min[mm]
             im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
               /data, color=im_color(color[mm]), thick=5
          endfor
; priors
          for ii = 0, ngd-1 do begin
             xwid = fit[gd[ii]].sf_sigz*1.5
             yhei = fit[gd[ii]].sf_sfrm_priors_max[mm]-fit[gd[ii]].sf_sfrm_priors_min[mm]
             xcen = fit[gd[ii]].sf_medz
             ycen = yhei/2.0+fit[gd[ii]].sf_sfrm_priors_min[mm]
             im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
               /data, color=im_color(color[mm]), thick=5, line=1
          endfor
       endif
    endfor
;   djs_oplot, zaxis, -alog10(ageaxis), line=5, thick=5, color=im_color(keycolor)
    djs_oplot, zaxis, -alog10(0.5*ageaxis), line=0, thick=5, color=im_color(keycolor)

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

stop    
    
; --------------------------------------------------
; *relative* number and mass-density evolution for quiescent, and
; star-forming galaxies in bins of mass 
    all = mrdfits(mfpath+'rhoden_numden_vs_redshift_all.fits.gz',1)
    qq = mrdfits(mfpath+'rhoden_numden_vs_redshift_quiescent.fits.gz',1)
    sf = mrdfits(mfpath+'rhoden_numden_vs_redshift_active.fits.gz',1)
    tot = [[qq],[sf],[all]]
    color = ['tomato','sky blue','dim grey']
    symcolor = ['dark red','midnight blue','dim grey']

    xrange = [0.001,1.01]
    rhorange = [-1.1,+0.65]
    numrange = [-1.1,+0.65]
    tmax = 1
    
    psfile = talkpath+'relative_rhoden_numden_bymass.ps'
    im_plotconfig, 17, pos, psfile=psfile, xmargin=[1.2,0.3], ymargin=[0.5,0.8], $
      width=2.9*[1,1,1], height=2.9*[1,1], charsize=1.5, yspace=0.1, xspace=[0.1,0.1], keynote=keynote

    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'

; #########################
; 10-10.5
    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
      ytitle=textoidl('log (n / n_{0})'), color=im_color(keycolor)
    im_xyouts_title, title='log (M/M'+sunsymbol()+')=10-10.5', color=im_color(keycolor)
    for ii = 0, tmax do begin
       gd = where(tot[*,ii].num_10_105 gt -900)
       zz = tot[gd,ii].zbin
       num = tot[gd,ii].num_10_105-tot[0,ii].num_105_11
       err = tot[gd,ii].numerr_cv_10_105
       polyfill, [zz,reverse(zz)],[num-err,reverse(num+err)], $
         /data, color=im_color(color[ii]), noclip=0, /fill
    endfor

    plot, [0], [0], /nodata, /noerase, position=pos[*,3], xsty=1, ysty=1, $
      yrange=rhorange, xrange=xrange, xtitle='Redshift', $
      ytitle=textoidl('log (\rho_{*} / \rho_{*,0})'), color=im_color(keycolor)
    for ii = 0, tmax do begin
       gd = where(tot[*,ii].rho_10_105 gt -900)
       zz = tot[gd,ii].zbin
       rho = tot[gd,ii].rho_10_105-tot[0,ii].rho_105_11
       err = tot[gd,ii].rhoerr_cv_10_105
       polyfill, [zz,reverse(zz)],[rho-err,reverse(rho+err)], $
         /data, color=im_color(color[ii]), noclip=0, /fill
    endfor
    
; #########################
; 10.5-11
    plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), color=im_color(keycolor)
    im_xyouts_title, title='log (M/M'+sunsymbol()+')=10.5-11', color=im_color(keycolor)
    for ii = 0, tmax do begin
       gd = where(tot[*,ii].num_105_11 gt -900)
       zz = tot[gd,ii].zbin
       num = tot[gd,ii].num_105_11-tot[0,ii].num_105_11
       err = tot[gd,ii].numerr_cv_105_11
       polyfill, [zz,reverse(zz)],[num-err,reverse(num+err)], $
         /data, color=im_color(color[ii]), noclip=0, /fill
    endfor

    plot, [0], [0], /nodata, /noerase, position=pos[*,4], xsty=1, ysty=1, $
      yrange=rhorange, xrange=xrange, xtitle='Redshift', ytitle='', $
      ytickname=replicate(' ',10), color=im_color(keycolor)
    for ii = 0, tmax do begin
       gd = where(tot[*,ii].rho_105_11 gt -900)
       zz = tot[gd,ii].zbin
       rho = tot[gd,ii].rho_105_11-tot[0,ii].rho_105_11
       err = tot[gd,ii].rhoerr_cv_105_11
       polyfill, [zz,reverse(zz)],[rho-err,reverse(rho+err)], $
         /data, color=im_color(color[ii]), noclip=0, /fill
    endfor
    
; #########################
; 11-11.5
    plot, [0], [0], /nodata, /noerase, position=pos[*,2], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), color=im_color(keycolor)
    im_xyouts_title, title='log (M/M'+sunsymbol()+')=11-11.5', color=im_color(keycolor)
    for ii = 0, tmax do begin
       gd = where(tot[*,ii].num_11_115 gt -900)
       zz = tot[gd,ii].zbin
       num = tot[gd,ii].num_11_115-tot[0,ii].num_105_11
       err = tot[gd,ii].numerr_cv_11_115
       polyfill, [zz,reverse(zz)],[num-err,reverse(num+err)], $
         /data, color=im_color(color[ii]), noclip=0, /fill
    endfor

    plot, [0], [0], /nodata, /noerase, position=pos[*,5], xsty=1, ysty=1, $
      yrange=rhorange, xrange=xrange, xtitle='Redshift', ytitle='', $
      ytickname=replicate(' ',10), color=im_color(keycolor)
    for ii = 0, tmax do begin
       gd = where(tot[*,ii].rho_11_115 gt -900)
       zz = tot[gd,ii].zbin
       rho = tot[gd,ii].rho_11_115-tot[0,ii].rho_105_11
       err = tot[gd,ii].rhoerr_cv_11_115
       polyfill, [zz,reverse(zz)],[rho-err,reverse(rho+err)], $
         /data, color=im_color(color[ii]), noclip=0, /fill
    endfor
    
    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

stop

; --------------------------------------------------
; primus - sfrm vs mass coded by galaxy type
    fit = mrdfits(mfpath+'sfrm_vs_redshift.fits.gz',1)
    all = mrdfits(mfspath+'mfdata_all_supergrid01.fits.gz',1) ; fiducial supergrid
    qq = mrdfits(mfspath+'mfdata_quiescent_supergrid01.fits.gz',1)
    sf = mrdfits(mfspath+'mfdata_active_supergrid01.fits.gz',1)
    zbins = mf_zbins(nzbins)

    psfile = talkpath+'sfrm_vs_mass.ps'
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.6, $
      height=2.8*[1,1,1], width=3.5*[1,1], keynote=keynote

    if keyword_set(keynote) then begin
       keycolor = 'white'
       qqcolor = 'tomato'
       sfcolor = 'sky blue'
    endif else begin
       keycolor = 'black'
       qqcolor = 'dark red'
       sfcolor = 'steel blue'
    endelse

    xrange = [8.2,12.2]
    yrange = [-4.7,+1.8]
    xtitle1 = mfplot_masstitle()
    ytitle1 = mfplot_sfrmtitle()
    maxis = range(xrange[0]+0.05,xrange[1]-0.05,50)
    
    for iz = 0, nzbins-1 do begin
       allindx = where((all.z ge zbins[iz].zlo) and (all.z lt zbins[iz].zup),nall)
       qqindx = where((qq.z ge zbins[iz].zlo) and (qq.z lt zbins[iz].zup))
       sfindx = where((sf.z ge zbins[iz].zlo) and (sf.z lt zbins[iz].zup))

       sf_mass = sf[sfindx].mass
       sf_sfrm = sf[sfindx].sfrm100>(-4)
       sf_weight = sf[sfindx].weight

       qq_mass = qq[qqindx].mass
       qq_sfrm = qq[qqindx].sfrm100>(-4)
       qq_weight = qq[qqindx].weight
       
       all_mass = all[allindx].mass
       all_sfrm = all[allindx].sfrm100>(-4)
       all_weight = all[allindx].weight
       
       if odd(iz) then begin
          ytitle = '' & ytickname = replicate(' ',10)
       endif else begin
          ytitle = ytitle1 & delvarx, ytickname
       endelse
       if (iz lt 4) then xtickname = replicate(' ',10) else delvarx, xtickname

       mfplot_scatterplot, sf_mass, sf_sfrm, weight=sf_weight, noerase=(iz gt 0), $
         position=pos[*,iz], xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
         xtitle='', ytitle=ytitle, ccolor=im_color(sfcolor,10), $
         /nogrey, outcolor=im_color(sfcolor,11), xtickname=xtickname, $
         ytickname=ytickname, npix=16, color=im_color(keycolor,5)
       mfplot_scatterplot, qq_mass, qq_sfrm, weight=qq_weight, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
         xtitle='', ytitle='', /noerase, position=pos[*,iz], ccolor=im_color(qqcolor,12), $
         /nogrey, outcolor=im_color(qqcolor,12), npix=16, color=im_color(keycolor,5)
       djs_oplot, sf[sfindx[0]].masslimit*[1,1], !y.crange, line=2, color=im_color(keycolor,5)
       djs_oplot, sf[sfindx[0]].masslimit_50*[1,1], !y.crange, line=1, color=im_color(keycolor,5)
;      djs_oplot, maxis, -maxis+9, line=5 ; =1 Msun/yr
       djs_oplot, maxis, poly(maxis-fit[iz].mpivot,fit[iz].coeff), line=0, thick=6, color=im_color(keycolor,5)
;      djs_oplot, maxis, poly(maxis-fit[iz].mpivot,fit[iz].coeff)-4.0*fit[iz].sigma, line=0, thick=6
;      djs_oplot, maxis, poly(maxis-fit[iz].mpivot,fit[iz].coeff)+2.0*fit[iz].sigma, line=0, thick=6
       legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[iz].zup,format='(F12.2)'),2), $
         /right, /top, box=0, margin=0, charsize=1.4, textcolor=im_color(keycolor,5)

;      djs_oplot, maxis, poly(maxis,[-6.19,+0.67-1])+9, line=3, thick=5 ; Noeske+07
       djs_oplot, maxis, poly(maxis,[-6.33,-0.35])+9, line=5, thick=6, color=im_color(keycolor,5) ; Salim+07
;      djs_oplot, maxis, poly(maxis,[-6.33,-0.35])+9+alog10(5.0), line=5, thick=5    ; factor of 3

    endfor 
    xyouts, pos[0,5], pos[1,5]-0.07, align=0.5, /norm, xtitle1, color=im_color(keycolor,5)

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

; --------------------------------------------------
; number and mass-density evolution for all, quiescent, and
; star-forming galaxies in bins of mass
    all = mrdfits(mfpath+'rhoden_numden_vs_redshift_all.fits.gz',1)
    qq = mrdfits(mfpath+'rhoden_numden_vs_redshift_quiescent.fits.gz',1)
    sf = mrdfits(mfpath+'rhoden_numden_vs_redshift_active.fits.gz',1)

    tot = [[all],[qq],[sf]]
    psym = [15,14,16]
    if keyword_set(keynote) then $
      color = ['white','tomato','sky blue'] else $
        color = ['black','dark red','steel blue']
    symsize = [2.0,2.2,2.0]
    off = [0.0,+0.01,-0.01]

    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'

    allerrors = 0
    if allerrors then $
      psfile = talkpath+'rhoden_numden_bymass_allerrors.ps' else $
        psfile = talkpath+'rhoden_numden_bymass.ps'
    im_plotconfig, 17, pos, psfile=psfile, xmargin=[1.2,0.3], ymargin=[0.7,0.85], $
      width=2.9*[1,1,1], height=[3.0,3.0], charsize=1.6, yspace=0.1, xspace=[0.1,0.1], keynote=keynote

    xrange = [0.001,1.01]
    rhorange = [6.6,8.3]
    numrange = [-4.3,-2.35]

; #########################
; 10-10.5
    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
      ytitle=mfplot_numdentitle(), color=im_color(keycolor)
    im_xyouts_title, title='log (M/M'+sunsymbol()+')=10-10.5', color=im_color(keycolor)
    im_legend, ['All','Quiescent','Star-Forming'], /left, /bottom, box=0, $
      charsize=1.2, psym=psym, symsize=symsize*0.9, color=color, textcolor=im_color(keycolor)

    for ii = 0, 2 do begin
       gd = where(tot[*,ii].num_10_105 gt -999.0,ngd)
       if (ngd ne 0) then begin
          djs_oploterr, tot[gd,ii].zbin+off[ii], tot[gd,ii].num_10_105, xerr=tot[gd,ii].zsigma*0, $
            yerr=tot[gd,ii].numerr_cv_10_105, $
            psym=symcat(psym[ii],color=im_color(keycolor)), symsize=symsize[ii], $
            color=im_color(color[ii]), errcolor=im_color(color[ii]), errstyle=0, errthick=6, /cap
          if allerrors then begin
; models             
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.9
                yhei = tot[gd[jj],ii].num_models_max_10_105-tot[gd[jj],ii].num_models_min_10_105
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].num_models_min_10_105
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=3
             endfor
; priors
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.5
                yhei = tot[gd[jj],ii].num_priors_max_10_105-tot[gd[jj],ii].num_priors_min_10_105
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].num_priors_min_10_105
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=5, line=1
             endfor
          endif
       endif
    endfor 

    plot, [0], [0], /nodata, /noeras, position=pos[*,3], xsty=1, ysty=1, $
      yrange=rhorange, xrange=xrange, xtitle='Redshift', $
      ytitle=mfplot_rhodentitle(), color=im_color(keycolor)
    for ii = 0, 2 do begin
       gd = where(tot[*,ii].rho_10_105 gt -999.0,ngd)
       if (ngd ne 0) then begin
          djs_oploterr, tot[gd,ii].zbin+off[ii], tot[gd,ii].rho_10_105, xerr=tot[gd,ii].zsigma*0, yerr=tot[gd,ii].rhoerr_cv_10_105, $
            psym=symcat(psym[ii],color=im_color(keycolor)), symsize=symsize[ii], color=im_color(color[ii]), $
            errcolor=im_color(color[ii]), errstyle=0, errthick=6, /cap
          if allerrors then begin
; models             
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.9
                yhei = tot[gd[jj],ii].rho_models_max_10_105-tot[gd[jj],ii].rho_models_min_10_105
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].rho_models_min_10_105
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=3
             endfor
; priors
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.5
                yhei = tot[gd[jj],ii].rho_priors_max_10_105-tot[gd[jj],ii].rho_priors_min_10_105
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].rho_priors_min_10_105
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=5, line=1
             endfor
          endif
       endif
    endfor 
    
; #########################
; 10.5-11
    plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), ytitle='', color=im_color(keycolor)
    im_xyouts_title, title='log (M/M'+sunsymbol()+')=10.5-11', color=im_color(keycolor)
;   djs_oplot, combo.z, combo.num_105_11, psym=6, color='green', symsize=2.0
;   djs_oplot, ilbert_qq.z, ilbert_qq.num_105_11, psym=6, color='orange', symsize=2.0
    for ii = 0, 2 do begin
       gd = where(tot[*,ii].num_105_11 gt -999.0,ngd)
       if (ngd ne 0) then begin
          djs_oploterr, tot[gd,ii].zbin+off[ii], tot[gd,ii].num_105_11, xerr=tot[gd,ii].zsigma*0, yerr=tot[gd,ii].numerr_cv_105_11, $
            psym=symcat(psym[ii],color=im_color(keycolor)), symsize=symsize[ii], color=im_color(color[ii]), $
            errcolor=im_color(color[ii]), errstyle=0, errthick=6, /cap
          if allerrors then begin
; models             
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.9
                yhei = tot[gd[jj],ii].num_models_max_105_11-tot[gd[jj],ii].num_models_min_105_11
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].num_models_min_105_11
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=3
             endfor
; priors
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.5
                yhei = tot[gd[jj],ii].num_priors_max_105_11-tot[gd[jj],ii].num_priors_min_105_11
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].num_priors_min_105_11
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=5, line=1
             endfor
          endif
       endif
    endfor 

    plot, [0], [0], /nodata, /noerase, position=pos[*,4], xsty=1, ysty=1, $
      yrange=rhorange, xrange=xrange, xtitle='Redshift', ytickname=replicate(' ',10), ytitle='', $
      color=im_color(keycolor)
;   djs_oplot, combo.z, combo.rho_105_11, psym=6, color='green', symsize=2.0
;   djs_oplot, ilbert_qq.z, ilbert_qq.rho_105_11, psym=6, color='orange', symsize=2.0
    for ii = 0, 2 do begin
       gd = where(tot[*,ii].rho_105_11 gt -999.0,ngd)
       if (ngd ne 0) then begin
          djs_oploterr, tot[gd,ii].zbin+off[ii], tot[gd,ii].rho_105_11, xerr=tot[gd,ii].zsigma*0, yerr=tot[gd,ii].rhoerr_cv_105_11, $
            psym=symcat(psym[ii],color=im_color(keycolor)), symsize=symsize[ii], color=im_color(color[ii]), $
            errcolor=im_color(color[ii]), errstyle=0, errthick=6, /cap
          if allerrors then begin
; models             
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.9
                yhei = tot[gd[jj],ii].rho_models_max_105_11-tot[gd[jj],ii].rho_models_min_105_11
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].rho_models_min_105_11
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=3
             endfor
; priors
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.5
                yhei = tot[gd[jj],ii].rho_priors_max_105_11-tot[gd[jj],ii].rho_priors_min_105_11
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].rho_priors_min_105_11
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=5, line=1
             endfor
          endif
       endif
    endfor 

; #########################
; 11-11.5
    plot, [0], [0], /nodata, /noerase, position=pos[*,2], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), ytitle='', color=im_color(keycolor)
    im_xyouts_title, title='log (M/M'+sunsymbol()+')=11-11.5', color=im_color(keycolor)
;   djs_oplot, combo.z, combo.num_11_115, psym=6, color='green', symsize=2.0
;   djs_oplot, ilbert_qq.z, ilbert_qq.num_11_115, psym=6, color='orange', symsize=2.0
    for ii = 0, 2 do begin
       gd = where(tot[*,ii].num_11_115 gt -999.0,ngd)
       if (ngd ne 0) then begin
          djs_oploterr, tot[gd,ii].zbin+off[ii], tot[gd,ii].num_11_115, xerr=tot[gd,ii].zsigma*0, yerr=tot[gd,ii].numerr_cv_11_115, $
            psym=symcat(psym[ii],color=im_color(keycolor)), symsize=symsize[ii], $
            color=im_color(color[ii]), errcolor=im_color(color[ii]), errstyle=0, errthick=6, /cap
          if allerrors then begin
; models             
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.9
                yhei = tot[gd[jj],ii].num_models_max_11_115-tot[gd[jj],ii].num_models_min_11_115
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].num_models_min_11_115
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=3
             endfor
; priors
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.5
                yhei = tot[gd[jj],ii].num_priors_max_11_115-tot[gd[jj],ii].num_priors_min_11_115
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].num_priors_min_11_115
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=5, line=1
             endfor
          endif
       endif
    endfor 

    plot, [0], [0], /nodata, /noerase, position=pos[*,5], xsty=1, ysty=1, $
      yrange=rhorange, xrange=xrange, xtitle='Redshift', ytickname=replicate(' ',10), ytitle='', $
      color=im_color(keycolor)
;   djs_oplot, combo.z, combo.rho_11_115, psym=6, color='green', symsize=2.0
;   djs_oplot, ilbert_qq.z, ilbert_qq.rho_11_115, psym=6, color='orange', symsize=2.0
    for ii = 0, 2 do begin
       gd = where(tot[*,ii].rho_11_115 gt -999.0,ngd)
       if (ngd ne 0) then begin
          djs_oploterr, tot[gd,ii].zbin+off[ii], tot[gd,ii].rho_11_115, xerr=tot[gd,ii].zsigma*0, yerr=tot[gd,ii].rhoerr_cv_11_115, $
            psym=symcat(psym[ii],color=im_color(keycolor)), symsize=symsize[ii], $
            color=im_color(color[ii]), errcolor=im_color(color[ii]), errstyle=0, errthick=6, /cap
          if allerrors then begin
; models             
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.9
                yhei = tot[gd[jj],ii].rho_models_max_11_115-tot[gd[jj],ii].rho_models_min_11_115
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].rho_models_min_11_115
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=3
             endfor
; priors
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.5
                yhei = tot[gd[jj],ii].rho_priors_max_11_115-tot[gd[jj],ii].rho_priors_min_11_115
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].rho_priors_min_11_115
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=5, line=1
             endfor
          endif
       endif
    endfor 

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

; ---------------------------------------------------------------------------
; 6x3 panel plot showing the field-averaged MFs for 'all',
; 'quiescent', and 'star-forming galaxies'
    supergrid = 1
    xrange = [8.05,12.5]
;   xrange = [8.3,12.5]
    yrange = [-5.2,-1.1]
    showerr = 0
    zbins = mf_zbins(nzbins)
    maxis1 = range(8.8,12.5,100)

    nlimitshow = 3
    pcolor = 'black' & pline = 0
    scolor = 'black' & sline = 5

    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'
    
    psfile = talkpath+'mf_vmax_everything.ps'
    im_plotconfig, 20, pos, psfile=psfile, charsize=1.4, $
      height=1.9*[1,1,1], keynote=keynote

; 'all' galaxies
    mfdata = read_vmax_mf(noevol=noevol,/avgfield,supergrid=supergrid,/log)
    mffit = read_vmax_mf(/bestfit,noevol=noevol,/avgfield,supergrid=supergrid)
    mffit_sdss = read_vmax_mf('sdss',/bestfit,noevol=noevol,supergrid=supergrid)
    for iz = 0, nzbins-1 do begin
       if (iz gt 0) then begin
          ytitle = ''
          ytickname = replicate(' ',10)
       endif else begin
          ytitle = '' ; mfplot_phititle()
          delvarx, ytickname
       endelse
       plot, [0], [0], /nodata, noerase=(iz gt 0), position=pos[*,iz+nzbins*0], xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, xtickname=replicate(' ',10), ytickname=ytickname, xtitle='', $
         ytitle=ytitle, xtickinterval=1, xminor=4, color=im_color(keycolor)
       if (iz eq 0) then legend, 'All', /left, /bottom, box=0, charsize=1.0, margin=0, textcolor=im_color(keycolor)
       legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[iz].zup,format='(F12.2)'),2), $
         /right, /top, box=0, margin=0, charsize=1.0, textcolor=im_color(keycolor)

       good = where(mfdata[iz].limit eq 1)
       limit = where(mfdata[iz].limit eq 0 and mfdata[iz].number gt 0 and mfdata[iz].mass gt xrange[0]+0.05,nlimit)
       if nlimit gt 0 then limit = (reverse(limit))[0:nlimitshow<nlimit-1]
       if keyword_set(showerr) then begin
          oploterror, mfdata[iz].mass[good], mfdata[iz].phi[good], $
            mfdata[iz].phierr[good], psym=symcat(8,thick=3), symsize=0.75, $
            errthick=!p.thick
       endif else begin
          phimass = mfdata[iz].mass[good]
; Poisson+CV errors
          phimin = mfdata[iz].phi[good]-sqrt(mfdata[iz].phierr_lower[good]^2+mfdata[iz].phierr_cv[good]^2)
          phimax = mfdata[iz].phi[good]+sqrt(mfdata[iz].phierr_upper[good]^2+mfdata[iz].phierr_cv[good]^2)
          polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
            /data, color=im_color('dim grey'), noclip=0, /fill
; Poisson errors
          phimin = mfdata[iz].phi[good]-mfdata[iz].phierr_lower[good]
          phimax = mfdata[iz].phi[good]+mfdata[iz].phierr_upper[good]
          polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
            /data, color=im_color('light gray'), noclip=0, /fill

          if (nlimit gt 0) then begin
             oploterror, mfdata[iz].mass[limit], mfdata[iz].phi[limit], mfdata[iz].phierr_upper[limit], $
               psym=symcat(9,thick=3), symsize=1, /nohat, /hibar, $
               color=im_color('gray40'), errcolor=im_color('gray40')
             oploterror, mfdata[iz].mass[limit], mfdata[iz].phi[limit], mfdata[iz].phierr_lower[limit], $
               psym=3, symsize=1, /nohat, /lobar, color=im_color('gray40'), errcolor=im_color('gray40')
          endif
       endelse
       djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=mffit_sdss)), $
         line=sline, thick=4
;      djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=mffit[iz])), $
;        line=pline, thick=4, color=im_color(pcolor)
    endfor 

; 'active' galaxies
    mfdata = read_vmax_mf(noevol=noevol,/avgfield,supergrid=supergrid,/active,/log)
    mffit = read_vmax_mf(/bestfit,noevol=noevol,/avgfield,supergrid=supergrid,/active)
    mffit_sdss = read_vmax_mf('sdss',/bestfit,noevol=noevol,supergrid=supergrid,/active)
    for iz = 0, nzbins-1 do begin
       if (iz gt 0) then begin
          ytitle = ''
          ytickname = replicate(' ',10)
       endif else begin
          ytitle = mfplot_phititle()
          delvarx, ytickname
       endelse
       plot, [0], [0], /nodata, /noerase, position=pos[*,iz+nzbins*1], xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, xtickname=replicate(' ',10), ytickname=ytickname, xtitle='', $
         ytitle=ytitle, xtickinterval=1, xminor=4, color=im_color(keycolor)
       if (iz eq 0) then legend, 'Star-Forming', /left, /bottom, box=0, $
         charsize=1.0, margin=0, textcolor=im_color(keycolor)

       good = where(mfdata[iz].limit eq 1)
       limit = where(mfdata[iz].limit eq 0 and mfdata[iz].number gt 0 and mfdata[iz].mass gt xrange[0]+0.05,nlimit)
       if nlimit gt 0 then limit = (reverse(limit))[0:nlimitshow<nlimit-1]
       if keyword_set(showerr) then begin
          oploterror, mfdata[iz].mass[good], mfdata[iz].phi[good], $
            mfdata[iz].phierr[good], psym=symcat(8,thick=3), symsize=0.75, $
            errthick=!p.thick
       endif else begin
;         djs_oplot, mfdata[iz].mass[good], mfdata[iz].phi[good], $
;           psym=symcat(8,thick=3), symsize=0.75
          phimass = mfdata[iz].mass[good]
; Poisson+CV errors
          phimin = mfdata[iz].phi[good]-sqrt(mfdata[iz].phierr_lower[good]^2+mfdata[iz].phierr_cv[good]^2)
          phimax = mfdata[iz].phi[good]+sqrt(mfdata[iz].phierr_upper[good]^2+mfdata[iz].phierr_cv[good]^2)
          polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
            /data, color=im_color('blue'), noclip=0, /fill
; Poisson errors
          phimin = mfdata[iz].phi[good]-mfdata[iz].phierr_lower[good]
          phimax = mfdata[iz].phi[good]+mfdata[iz].phierr_upper[good]
          polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
            /data, color=im_color('sky blue'), noclip=0, /fill

          if (nlimit gt 0) then begin
             oploterror, mfdata[iz].mass[limit], mfdata[iz].phi[limit], mfdata[iz].phierr_upper[limit], $
               psym=symcat(9,thick=3), symsize=1, /nohat, /hibar, $
               color=im_color('steel blue'), errcolor=im_color('steel blue')
             oploterror, mfdata[iz].mass[limit], mfdata[iz].phi[limit], mfdata[iz].phierr_lower[limit], $
               psym=3, symsize=1, /nohat, /lobar, color=im_color('steel blue'), errcolor=im_color('steel blue')
          endif
       endelse
       djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=mffit_sdss)), $
         line=sline, thick=4
;      djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=mffit[iz])), $
;        line=pline, thick=4, color=im_color(pcolor)
    endfor 
    
; 'quiescent' galaxies
    mfdata = read_vmax_mf(noevol=noevol,/avgfield,supergrid=supergrid,/quiescent,/log)
    mffit = read_vmax_mf(/bestfit,noevol=noevol,/avgfield,supergrid=supergrid,/quiescent)
    mffit_sdss = read_vmax_mf('sdss',/bestfit,noevol=noevol,supergrid=supergrid,/quiescent)
    for iz = 0, nzbins-1 do begin
       if (iz gt 0) then begin
          ytitle = ''
          ytickname = replicate(' ',10)
       endif else begin
          ytitle = '' ; mfplot_phititle()
          delvarx, ytickname
       endelse
       plot, [0], [0], /nodata, /noerase, position=pos[*,iz+nzbins*2], xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, ytickname=ytickname, $
         xtitle='', ytitle=ytitle, xtickinterval=1, xminor=4, color=im_color(keycolor)
       if (iz eq 0) then legend, 'Quiescent', /left, /bottom, box=0, charsize=1.0, $
         margin=0, textcolor=im_color(keycolor)

       good = where(mfdata[iz].limit eq 1)
       limit = where(mfdata[iz].limit eq 0 and mfdata[iz].number gt 0 and mfdata[iz].mass gt xrange[0]+0.05,nlimit)
       if nlimit gt 0 then limit = (reverse(limit))[0:nlimitshow<nlimit-1]
       if keyword_set(showerr) then begin
          oploterror, mfdata[iz].mass[good], mfdata[iz].phi[good], $
            mfdata[iz].phierr[good], psym=symcat(8,thick=3), symsize=0.75, $
            errthick=!p.thick
       endif else begin
;         djs_oplot, mass, phi, psym=symcat(8,thick=3), symsize=0.75
          phimass = mfdata[iz].mass[good]
; Poisson+CV errors
          phimin = mfdata[iz].phi[good]-sqrt(mfdata[iz].phierr_lower[good]^2+mfdata[iz].phierr_cv[good]^2)
          phimax = mfdata[iz].phi[good]+sqrt(mfdata[iz].phierr_upper[good]^2+mfdata[iz].phierr_cv[good]^2)
          polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
            /data, color=im_color('dark red'), noclip=0, /fill
; Poisson errors
          phimin = mfdata[iz].phi[good]-mfdata[iz].phierr_lower[good]
          phimax = mfdata[iz].phi[good]+mfdata[iz].phierr_upper[good]
          polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
            /data, color=im_color('tomato'), noclip=0, /fill

          if (nlimit gt 0) then begin
             oploterror, mfdata[iz].mass[limit], mfdata[iz].phi[limit], mfdata[iz].phierr_upper[limit], $
               psym=symcat(9,thick=3), symsize=1, /nohat, /hibar, $
               color=im_color('orange red'), errcolor=im_color('orange red')
             oploterror, mfdata[iz].mass[limit], mfdata[iz].phi[limit], mfdata[iz].phierr_lower[limit], $
               psym=3, symsize=1, /nohat, /lobar, color=im_color('orange red'), errcolor=im_color('orange red')
          endif
       endelse
       djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=mffit_sdss)), $
         line=sline, thick=4
;      djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=mffit[iz])), $
;        line=pline, thick=4, color=im_color(pcolor)
    endfor 

; x-title
    xyouts, pos[0,15], pos[1,15]-0.08, mfplot_masstitle(), align=0.5, $
      /normal, color=im_color(keycolor)
    
    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

; ---------------------------------------------------------------------------
; sdss - compare all, quiescent, and active samples
    supergrid = 1
    psfile = talkpath+'mf_vmax_sdss.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=6.5, charsize=2, keynote=keynote

    if keyword_set(keynote) then keycolor = im_color('white',10) else $
      keycolor = im_color('black',10)

    xrange = [9,12.1]
    yrange = [-6.2,-1.7]
    maxis1 = range(xrange[0]+0.15,12.5,100)

    subsample = ['active','quiescent','all']
    if keyword_set(keynote) then begin
       fitcolor = ['medium blue','red','white']
       mfcolor = ['dodger blue','tomato','antique white']
    endif else begin
       fitcolor = ['medium blue','red','black']
       mfcolor = ['dodger blue','tomato','dark grey']
    endelse

    fitcolor_double = 'black'
    mfline_double = 1
    
    psymsize = [1.2,1.3,1.1]*1.3
    psymsize_limit = [0.8,0.9,0.8]*1.2
;   psymsize = [2.0,2.5,2.0]*0.7
;   mfpsym = [9,4,6]
    mfpsym = [16,15,14]
    mfpsym_limit = [9,6,4]
    mfline = [5,3,0]
    nlimitshow = 2

    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytitle=mfplot_phititle(), color=keycolor
    im_legend, ['All','Quiescent','Star-Forming'], /left, /bottom, box=0, $
      color=reverse(mfcolor), psym=reverse(mfpsym), symsize=[1.7,2.1,1.9], $
      symthick=8, spacing=2.4, margin=0, charsize=1.7, textcolor=keycolor
    
; print the "quenching mass"
    qq = alog10(mf_schechter(maxis1,schechter=read_vmax_mf('sdss',/quie,/best)))
    aa = alog10(mf_schechter(maxis1,schechter=read_vmax_mf('sdss',/act,/best)) )
    splog, 'quenching mass ', interpol(maxis1,qq-aa,0.0)

    mffitall = maxis1*0.0 ; add the active and quiescent Schechter fits
    for jj = 0, n_elements(subsample)-1 do begin
       case jj of
          0: begin
             mfdata = read_vmax_mf('sdss',supergrid=supergrid,/active,/log)
             mffit = read_vmax_mf('sdss',supergrid=supergrid,/active,/bestfit)
             model = mf_schechter(maxis1,schechter=mffit)
             mffitall = model
          end
          1: begin
             mfdata = read_vmax_mf('sdss',supergrid=supergrid,/quiescent,/log)
             mffit = read_vmax_mf('sdss',supergrid=supergrid,/quiescent,/bestfit)
             model = mf_schechter(maxis1,schechter=mffit)
             mffitall += model
          end
          2: begin
             mfdata = read_vmax_mf('sdss',supergrid=supergrid,/log)
             mffit = read_vmax_mf('sdss',supergrid=supergrid,/bestfit)
;            mffit_double = read_vmax_mf('sdss',supergrid=supergrid,/bestfit,/double)
             model = mf_schechter(maxis1,schechter=mffit)
;            model_double = mf_schechter_plus(maxis1,schechter=mffit_double)
          end
       endcase
       
;      good = where(mfdata.number gt 0 and mfdata.mass gt xrange[0]+0.08,nthese)
       good = where(mfdata.limit eq 1,nthese)
       limit = where(mfdata.number gt 0 and mfdata.limit eq 0 and mfdata.mass gt xrange[0]+0.08,nlimit)

       oploterror, mfdata.mass[good], mfdata.phi[good], mfdata.phierr_upper[good], $
         psym=symcat(mfpsym[jj],thick=5), symsize=psymsize[jj], errthick=5, $
         color=im_color(mfcolor[jj],101), errcolor=im_color(mfcolor[jj],101), /hibar, /nohat
       oploterror, mfdata.mass[good], mfdata.phi[good], mfdata.phierr_lower[good], psym=3, $
         symsize=psymsize[jj], errthick=5, color=im_color(mfcolor[jj],101), $
         errcolor=im_color(mfcolor[jj],101), /lobar, /nohat

       if nlimit gt 0 then limit = (reverse(limit))[0:nlimitshow<nlimit-1]
;      if (nlimit gt 0) then begin
;         oploterror, mfdata.mass[limit], mfdata.phi[limit], mfdata.phierr_upper[limit], $
;           psym=symcat(mfpsym_limit[jj],thick=2), symsize=psymsize_limit[jj], errthick=5, $
;           color=im_color(mfcolor[jj]), errcolor=im_color(mfcolor[jj]), /hibar, /nohat
;         oploterror, mfdata.mass[limit], mfdata.phi[limit], mfdata.phierr_upper[limit], psym=3, $
;           symsize=psymsize[jj], errthick=5, color=im_color(mfcolor[jj]), $
;           errcolor=im_color(mfcolor[jj]), /lobar, /nohat
;      endif
       
;       if (nlimit gt 0) then begin
;          limit = reverse(these[0]-(lindgen(nlimit)+1))
;          plotsym, 2, thick=6
;          djs_oplot, mfdata.mass[limit], alog10(mfdata.phi[limit]), psym=8, $
;            color=im_color(mfcolor[jj],100+jj), symsize=psymsize[jj]
;       endif
       
;       if subsample[jj] eq 'all' then begin
;;         djs_oplot, maxis1, alog10(mffitall), line=mfline[jj], $
;;           thick=7, color=im_color('dark green')
;          djs_oplot, maxis1, alog10(model), line=mfline[jj], thick=7, color=im_color(fitcolor[jj])
;;         djs_oplot, maxis1, alog10(model_double), line=mfline_double, thick=7, color=im_color(fitcolor_double)
;       endif else begin
;          djs_oplot, maxis1, alog10(model), line=mfline[jj], thick=7, color=im_color(fitcolor[jj])
;       endelse
    endfor
;   mfoplot_lit, maxis1, color=im_color('dark green'), line=5, thick=8, /baldry, /salp, /log

;; add an inset showing the M*-alpha covariance
;    scale = 100
;    nmonte = 500
;    inset_pos = [pos[0]+0.15,pos[1]+0.15,pos[0]+0.45,pos[1]+0.4]
;    djs_plot, [0], [0], /nodata, /noerase, position=inset_pos, $
;      xsty=1, ysty=1, xrange=scale*0.025*[-1,1], yrange=scale*0.025*[-1,1], $
;      xtitle='10^{-2} \Delta'+'log (M_{*}/M'+sunsymbol()+')', $
;      ytitle='10^{-2} \Delta\alpha', charsize=1.4, $
;      xtickinterval=1, ytickinterval=1
;    for jj = 0, n_elements(subsample)-1 do begin
;       mffit = read_vmax_mf('sdss',quiescent=(subsample[jj] eq 'quiescent'),$
;         active=(subsample[jj] eq 'active'),supergrid=supergrid,/bestfit)
;       covar = mffit.covar[1:2,1:2]*sqrt(mffit.chi2_dof)
;
;       rand = mrandomn(seed,covar,nmonte)
;       rand[*,0] = rand[*,0];+mffit.logmstar
;       rand[*,1] = rand[*,1];+mffit.alpha
;       ell1 = covar2ellipse(covar,nsigma=1.0)
;       ell3 = covar2ellipse(covar,nsigma=3.0)
;       tvellipse, scale*ell1.major, scale*ell1.minor, 0.0, 0.0, ell1.angle, $
;         thick=7, color=im_color(fitcolor[jj],101), /data, line=mfline[jj]
;;      tvellipse, scale*ell3.major, scale*ell3.minor, 0.0, 0.0, ell1.angle, $
;;        line=5, thick=4, color=im_color(fitcolor[jj],101), /data
;    endfor

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

; --------------------------------------------------
; primus - NUV-r vs r-J
    zbins = mf_zbins(nzbins)
    all = mrdfits(mfspath+'mfdata_all_supergrid01.fits.gz',1)
    sfrm = all.sfr100+9-all.mass ; [Gyr^-1]
    
    psfile = talkpath+'nuvmr_vs_rmj.ps'
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.6, $
      height=2.8*[1,1,1], width=3.5*[1,1], keynote=keynote

    if keyword_set(keynote) then keycolor = im_color('white',10) else $
      keycolor = im_color('black',10)
    
    xrange = [-0.2,1.9]
    yrange = [0.01,7.99]
    xtitle1 = textoidl('^{0.1}(r - J)')
    ytitle1 = textoidl('^{0.1}(NUV - r)')

    rmjaxis = range(xrange[0],xrange[1],50)
    junk = get_mf_quiescent(box=box)

    for iz = 0, nzbins-1 do begin
       if odd(iz) then begin
          ytitle = '' & ytickname = replicate(' ',10)
       endif else begin
          ytitle = ytitle1 & delvarx, ytickname
       endelse
       if (iz lt 4) then xtickname = replicate(' ',10) else delvarx, xtickname

       indx = where((all.z ge zbins[iz].zlo) and (all.z lt zbins[iz].zup),nindx)
       mass = all[indx].mass
       sfrm = (all[indx].sfr100+9-mass)>(-4)
;      qq = where(sfrm lt (poly(mass-fit[iz].mpivot,fit[iz].coeff)-4.0*fit[iz].sigma))

       isqq = where(all[indx].nuvmr_01 gt poly(all[indx].rmj_01-box.rmj_turn,$
         [box.nuvmr_const,box.rmj_slope])>box.nuvmr_const,nisqq)
       splog, nisqq, nindx, nisqq/float(nindx)
       
       low = where(sfrm lt -3.0,nlow)
       high = where(sfrm gt -1.5,nhigh)
       interm = where(sfrm gt -3.0 and sfrm lt -1.5,ninterm)
       
       mfplot_scatterplot, all[indx[high]].rmj_01, all[indx[high]].nuvmr_01, $
         weight=all[indx[high]].weight, noerase=(iz gt 0), /nogrey, $
         position=pos[*,iz], xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
         xtitle='', ytitle=ytitle, xtickname=xtickname, $
         ytickname=ytickname, npix=16, ccolor=im_color('sky blue'), $
         outcolor=im_color('sky blue'), color=keycolor
       mfplot_scatterplot, all[indx[interm]].rmj_01, all[indx[interm]].nuvmr_01, $
         weight=all[indx[interm]].weight, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
         xtitle='', ytitle='', /noerase, position=pos[*,iz], ccolor=im_color('spring green1'), $
         /nogrey, outcolor=im_color('spring green1'), npix=16
       mfplot_scatterplot, all[indx[low]].rmj_01, all[indx[low]].nuvmr_01, $
         weight=all[indx[low]].weight, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
         xtitle='', ytitle='', /noerase, position=pos[*,iz], ccolor=im_color('tomato'), $
         /nogrey, outcolor=im_color('tomato'), npix=16

       legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[iz].zup,format='(F12.2)'),2), $
         /left, /top, box=0, margin=0, charsize=1.4, textcolor=keycolor

       djs_oplot, rmjaxis, poly(rmjaxis-box.rmj_turn,[box.nuvmr_const,box.rmj_slope])>box.nuvmr_const, $
         color=keycolor, line=5
;      djs_oplot, [pars.x1,pars.x2], [pars.y1,pars.y1], $
;        line=line, color=color, thick=thick
;      djs_oplot, [pars.x2,pars.x3], [pars.y1,pars.y2], $
;        line=line, color=color, thick=thick
;      djs_oplot, [pars.x3,pars.x3], [pars.y2,pars.y3], $
;        line=line, color=color, thick=thick
    endfor
    xyouts, pos[0,5], pos[1,5]-0.07, align=0.5, /norm, xtitle1, color=keycolor

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

return
end
