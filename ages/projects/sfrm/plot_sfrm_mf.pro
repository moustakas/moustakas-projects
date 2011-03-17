;+
; NAME:
;   PLOT_SFRM_MF
;
; PURPOSE:
;   Plot the mass function results from SFRM_MF.
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 05, UCSD
;
; Copyright (C) 2010, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro plot_mf, mf_fit, mf_data, zbins, sdss_mf_fit=sdss_mf_fit, $
  psfile=psfile
; simple wrapper to make plotting the MFs of the various subsamples
; (e.g., all, quiescent, active) easier 

    im_plotconfig, 7, pos, psfile=psfile, charsize=1.6, $
      height=3.0*[1,1,1]

    maxis1 = im_array(7.5,12.5,0.01)
    xrange = [8.0,12.5]
    yrange = [5E-7,0.2]
    
    for ii = 0, n_elements(mf_fit)-1 do begin
       if odd(ii) then begin
          ytitle = ''
          ytickname = replicate(' ',10)
       endif else begin
          if (ii eq 2) then $
            ytitle = sfrm_phititle() else delvarx, ytitle
          delvarx, ytickname 
       endelse

       if (ii lt 4) then begin
          xtitle = ''
          xtickname = replicate(' ',10)
       endif else begin
          delvarx, xtickname 
       endelse
       
       djs_plot, [0], [0], /nodata, noerase=(ii gt 0), $
         position=pos[*,ii], xsty=1, ysty=1, /ylog, $
         yrange=yrange, xrange=xrange, xtickname=xtickname, $
         ytickname=ytickname, xtitle=xtitle, ytitle=ytitle

; generate a shaded histogram showing the effects of cosmic variance
; and model assumptions       
       splog, 'Fix this!!'
       these = where((mf_data[ii].phi gt 0.0),these)
;      these = where((mf_data[ii].phi gt 0.0) and $
;        (mf_data[ii].phierr_cv gt -999.0) and $
;        (mf_data[ii].phierr_model gt -999.0),nthese)
       fullbin = mf_data[ii].fullbin[these]
       mass = mf_data[ii].mass[these]
       phi = mf_data[ii].phi[these]
       phierr = mf_data[ii].phierr[these]
       splog, 'Fix this!!'
       phierr_cv = 0.0
       phierr_model = 0.0
;      phierr_cv = mf_data[ii].phierr_cv[these]
;      phierr_model = mf_data[ii].phierr_model[these]

;      niceprint, mass, phi, phierr, phierr_cv, phierr_model
       allerr = sqrt(phierr^2+phierr_cv^2+phierr_model^2)
       phimin = phi-allerr
       phimax = phi+allerr

       polyfill, [mass,reverse(mass)],[phimin,reverse(phimax)], $
         /data, /fill, color=fsc_color('tan',100), $
         noclip=0
       
; plot the mf_data, distinguishing between objects that are above and
; below the 75% completeness limit
       these = where((mf_data[ii].phi gt 0.0),nthese)
       fullbin = mf_data[ii].fullbin[these]
       mass = mf_data[ii].mass[these]
       phi = mf_data[ii].phi[these]
       phierr = mf_data[ii].phierr[these]
       phierr_cv = mf_data[ii].phierr_cv[these]
       phierr_model = mf_data[ii].phierr_model[these]
       
       above = where(fullbin,comp=below)
       oploterror, mass[above], phi[above], phierr[above], $
         psym=symcat(9,thick=!p.thick), symsize=1.2
       oploterror, mass[below], phi[below], phierr[below], $
         psym=symcat(9,thick=!p.thick), symsize=1.2, color=djs_icolor('grey'), $
         errcolor=djs_icolor('grey')

; overplot the Schechter fit and the SDSS results for reference
       djs_oplot, maxis1, mf_schechter(10^maxis1,mf_fit[ii]), $
         line=0, color='red'
       djs_oplot, maxis1, mf_schechter(10^maxis1,sdss_mf_fit), $
         line=5, color=fsc_color('dodger blue',100)
       legend, 'z='+strtrim(string(zbins[ii].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[ii].zup,format='(F12.2)'),2), $
         /right, /top, box=0, margin=0, charsize=1.5
    endfor
    xyouts, pos[0,5], pos[1,5]-0.07, align=0.5, /norm, sfrm_masstitle()

    im_plotconfig, /psclose

return
end

function get_mstar, all=all, quiescent=quiescent, $
  active=active, mstar_err=mstar_err
    sfrmpath = ages_path(/projects)+'sfrm/'
    if keyword_set(all) then suffix = 'all'
    if keyword_set(quiescent) then suffix = 'quiescent'
    if keyword_set(active) then suffix = 'active'
    mf_fit = mrdfits(sfrmpath+'mf_fit_'+suffix+'.fits.gz',1,/silent)
    mstar = mf_fit.mstar
    splog, 'Fix this!'
    mstar_err = sqrt(mf_fit.mstar_err^2+mf_fit.mstar_cv_err^2+mf_fit.mstar_model_err^2)
return, mstar
end

function get_sdss_mstar, all=all, quiescent=quiescent, $
  active=active, mstar_err=mstar_err
    sfrmpath = ages_path(/projects)+'sfrm/'
    if keyword_set(all) then suffix = 'all'
    if keyword_set(quiescent) then suffix = 'quiescent'
    if keyword_set(active) then suffix = 'active'
    mf_fit = mrdfits(sfrmpath+'sdss_mf_fit_'+suffix+'.fits.gz',1,/silent)
    mstar = mf_fit.mstar
    mstar_err = mf_fit.mstar_err
return, mstar
end

function get_phistar, phiscale=phiscale, all=all, quiescent=quiescent, $
  active=active, phistar_err=phistar_err
    sfrmpath = ages_path(/projects)+'sfrm/'
    if keyword_set(all) then suffix = 'all'
    if keyword_set(quiescent) then suffix = 'quiescent'
    if keyword_set(active) then suffix = 'active'
    mf_fit = mrdfits(sfrmpath+'mf_fit_'+suffix+'.fits.gz',1,/silent)
    phistar = mf_fit.phistar
    phistar_err = sqrt(mf_fit.phistar_err^2+mf_fit.phistar_cv_err^2+$
      mf_fit.phistar_model_err^2)
    phistar = phistar*phiscale
    phistar_err = phistar_err*phiscale
return, phistar
end

function get_sdss_phistar, phiscale=phiscale, all=all, $
  quiescent=quiescent, active=active, phistar_err=phistar_err
    sfrmpath = ages_path(/projects)+'sfrm/'
    if keyword_set(all) then suffix = 'all'
    if keyword_set(quiescent) then suffix = 'quiescent'
    if keyword_set(active) then suffix = 'active'
    mf_fit = mrdfits(sfrmpath+'sdss_mf_fit_'+suffix+'.fits.gz',1,/silent)
    phistar = mf_fit.phistar
    phistar_err = mf_fit.phistar_err
    phistar = phistar*phiscale
    phistar_err = phistar_err*phiscale
return, phistar
end

pro plot_sfrm_mf, ps=ps
; jm10feb05ucsd - plot the stellar mass functions in AGES

    sfrmpath = ages_path(/projects)+'sfrm/'
    paperpath = ages_path(/papers)+'sfrm/'
    zbins = sfrm_zbins(nzbins)

    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'
    subsample = ['all','quiescent','active']

; --------------------------------------------------
; mass functions in six redshift bins for all/quiescent/active
; subsamples

    for jj = 0, n_elements(subsample)-1 do begin
       mf_data = mrdfits(sfrmpath+'mf_data_'+subsample[jj]+'.fits.gz',1,/silent)
       mf_fit = mrdfits(sfrmpath+'mf_fit_'+subsample[jj]+'.fits.gz',1,/silent)
       sdss_mf_fit = mrdfits(sfrmpath+'sdss_mf_fit_'+subsample[jj]+'.fits.gz',1,/silent)
       psfile = paperpath+'mf_'+subsample[jj]+suffix
       plot_mf, mf_fit, mf_data, zbins, sdss_mf_fit=sdss_mf_fit, psfile=psfile    
    endfor

; --------------------------------------------------
; M*, Phi* versus redshift
    psfile = paperpath+'mstar_phistar_vs_redshift'+suffix
    im_plotconfig, 6, pos, psfile=psfile, xmargin=[1.4,0.4], $
      width=6.7, height=[4.0,4.0]

    ztitle = 'Redshift'
    ytitle1 = sfrm_masstitle(/mstar)
    ytitle2 = sfrm_phititle(/phistar)

    phiscale = 1E3
    zrange = [-0.02,0.75]
    yrange1 = [10.5,11.1]
    yrange2 = [0,0.004]*phiscale

; -------------------------
; M* vs redshift    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=yrange1, xrange=zrange, xtitle='', ytitle=ytitle1, $
      ytickinterval=0.2, xtickname=replicate(' ',10)

; all galaxies
    mstar = get_mstar(/all,mstar_err=mstar_err)
    oploterror, zbins.zbin, mstar, mstar_err, $
      psym=-symcat(6,thick=8), symsize=2.5, errthick=8, line=0
    sdss_mstar = get_sdss_mstar(/all,mstar_err=sdss_mstar_err)
    oploterror, 0.025, sdss_mstar, sdss_mstar_err, $
      psym=-symcat(15,thick=8), symsize=2.5, errthick=8, line=0
; quiescent galaxies
    mstar = get_mstar(/quiescent,mstar_err=mstar_err)
    oploterror, zbins.zbin-0.01, mstar, mstar_err, $
      psym=-symcat(9,thick=8), symsize=2.5, errthick=8, line=3, $
      color=fsc_color('firebrick',100), errcolor=fsc_color('firebrick',100)
    sdss_mstar = get_sdss_mstar(/quiescent,mstar_err=sdss_mstar_err)
    oploterror, 0.025, sdss_mstar, sdss_mstar_err, $
      psym=-symcat(16,thick=8), symsize=2.5, errthick=8, line=3, $
      color=fsc_color('firebrick',100), errcolor=fsc_color('firebrick',100)
; star-forming galaxies
    mstar = get_mstar(/active,mstar_err=mstar_err)
    oploterror, zbins.zbin+0.01, mstar, mstar_err, $
      psym=-symcat(4,thick=8), symsize=3.5, errthick=8, line=5, $
      color=fsc_color('royal blue',100), errcolor=fsc_color('royal blue',100)
    sdss_mstar = get_mstar(/active,mstar_err=sdss_mstar_err)
    oploterror, 0.025, sdss_mstar, sdss_mstar_err, $
      psym=-symcat(14,thick=8), symsize=3.5, errthick=8, line=5, $
      color=fsc_color('royal blue',100), errcolor=fsc_color('royal blue',100)

; SDSS/Baldry & Glazebrook+08
;   plots, 0.025, 10.648, psym=symcat(15), symsize=3.0
;   oploterror, 0.025, 10.648, 0.013, psym=symcat(46), $
;     color=djs_icolor('orange'), errcolor=djs_icolor('orange'), $
;     symsize=5.0

    label = ['All','Quiescent','Star-Forming']
    im_legend, label, /right, /bottom, box=0, margin=0, $
      charsize=1.7, psym=[6,9,4], line=[0,3,5], $
      thick=8, symsize=[2.0,2.0,2.5], color=['','firebrick','royal blue'], $
      symthick=8
    
; -------------------------
; Phi* vs redshift    
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      yrange=yrange2, xrange=zrange, xtitle=ztitle, ytitle=ytitle2

; all galaxies
    phistar = get_phistar(/all,phiscale=phiscale,phistar_err=phistar_err)
    oploterror, zbins.zbin, phistar, phistar_err, $
      psym=-symcat(6,thick=8), symsize=2.5, errthick=8, line=0
    sdss_phistar = get_sdss_phistar(/all,phiscale=phiscale,phistar_err=sdss_phistar_err)
    oploterror, 0.025, sdss_phistar, sdss_phistar_err, $
      psym=-symcat(15,thick=8), symsize=2.5, errthick=8, line=0
; quiescent galaxies
    phistar = get_phistar(/quiescent,phiscale=phiscale,phistar_err=phistar_err)
    oploterror, zbins.zbin-0.01, phistar, phistar_err, $
      psym=-symcat(9,thick=8), symsize=2.5, errthick=8, line=3, $
      color=fsc_color('firebrick',100), errcolor=fsc_color('firebrick',100)
    sdss_phistar = get_sdss_phistar(/quiescent,phiscale=phiscale,phistar_err=sdss_phistar_err)
    oploterror, 0.025, sdss_phistar, sdss_phistar_err, $
      psym=-symcat(16,thick=8), symsize=2.5, errthick=8, line=3, $
      color=fsc_color('firebrick',100), errcolor=fsc_color('firebrick',100)
; all galaxies
    phistar = get_phistar(/active,phiscale=phiscale,phistar_err=phistar_err)
    oploterror, zbins.zbin+0.01, phistar, phistar_err, $
      psym=-symcat(4,thick=8), symsize=3.5, errthick=8, line=5, $
      color=fsc_color('royal blue',100), errcolor=fsc_color('royal blue',100)
    sdss_phistar = get_phistar(/active,phiscale=phiscale,phistar_err=sdss_phistar_err)
    oploterror, 0.025, sdss_phistar, sdss_phistar_err, $
      psym=-symcat(14,thick=8), symsize=3.5, errthick=8, line=5, $
      color=fsc_color('royal blue',100), errcolor=fsc_color('royal blue',100)

; SDSS/Baldry & Glazebrook+08
;   plots, 0.025, phiscale*4.26E-3, psym=symcat(15), symsize=3.0

    im_plotconfig, /psclose

;; --------------------------------------------------
;; M* versus redshift
;    psfile = paperpath+'mstar_vs_redshift'+suffix
;    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.4,0.4], $
;      width=6.7, height=5.0
;
;    zrange = [-0.02,0.75]
;    ztitle = 'Redshift'
;    yrange = [10.45,11.1]
;    ytitle = sfrm_masstitle(/mstar)
;    
;    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;      yrange=yrange, xrange=zrange, xtitle=ztitle, ytitle=ytitle, $
;      ytickinterval=0.2
;
;; all galaxies
;    mstar = get_mstar(/all,mstar_err=mstar_err)
;    oploterror, zbins.zbin, mstar, mstar_err, $
;      psym=-symcat(6,thick=8), symsize=2.5, errthick=8, line=0
;; quiescent galaxies
;    mstar = get_mstar(/quiescent,mstar_err=mstar_err)
;    oploterror, zbins.zbin-0.01, mstar, mstar_err, $
;      psym=-symcat(9,thick=8), symsize=2.5, errthick=8, line=3, $
;      color=fsc_color('firebrick',100), errcolor=fsc_color('firebrick',100)
;; all galaxies
;    mstar = get_mstar(/active,mstar_err=mstar_err)
;    oploterror, zbins.zbin+0.01, mstar, mstar_err, $
;      psym=-symcat(4,thick=8), symsize=3.5, errthick=8, line=5, $
;      color=fsc_color('royal blue',100), errcolor=fsc_color('royal blue',100)
;
;; SDSS/Baldry & Glazebrook+08
;    plots, 0.025, 10.648, psym=symcat(15), symsize=3.0
;;   oploterror, 0.025, 10.648, 0.013, psym=symcat(46), $
;;     color=djs_icolor('orange'), errcolor=djs_icolor('orange'), $
;;     symsize=5.0
;
;    label = ['All','Quiescent','Star-Forming']
;    im_legend, label, /right, /bottom, box=0, margin=0, $
;      charsize=1.7, psym=[6,9,4], line=[0,3,5], $
;      thick=8, symsize=[2.0,2.0,2.5], color=['','firebrick','royal blue'], $
;      symthick=8
;    
;    im_plotconfig, /psclose
;
;; --------------------------------------------------
;; Phi* versus redshift
;
;    psfile = paperpath+'phistar_vs_redshift'+suffix
;    im_plotconfig, 0, pos, psfile=psfile, $ ; xmargin=[1.5,0.3], $
;      width=6.7, height=5.0
;
;    scale = 1E3
;    yrange = [0,0.007]*scale
;    ytitle = textoidl('\Phi^{*}(!8M!6) (10^{-3} h_{70}^3 Mpc^{-3} dex^{-1})')
;
;    phistar = mf_fit.phistar
;    phistar_err = sqrt(mf_fit.phistar_err^2+mf_fit.phistar_cv_err^2+$
;      mf_fit.phistar_model_err^2)
;    
;    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;      yrange=yrange, xrange=zrange, xtitle=ztitle, ytitle=ytitle
;    oploterror, zbins.zbin, scale*phistar, scale*phistar_err, $
;      psym=symcat(6,thick=8), symsize=2.5, errthick=8
;; SDSS/Baldry & Glazebrook+08
;;   oploterror, 0.025, 10^(-0.46)*1D-3, 0.013, psym=symcat(46), $
;;     color=djs_icolor('orange'), errcolor=djs_icolor('orange'), $
;;     symsize=5.0    
;    plots, 0.025, scale*4.26E-3, psym=symcat(46), $
;      color=djs_icolor('orange'), symsize=5.0    
;    
;    im_plotconfig, /psclose

return
end
    
