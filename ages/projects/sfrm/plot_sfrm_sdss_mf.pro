;+
; NAME:
;   PLOT_SFRM_SDSS_MF
;
; PURPOSE:
;   Plot the SDSS (local) mass function results from SFRM_SDSS_MF. 
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 17, UCSD
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

pro plot_sfrm_sdss_mf, ps=ps

    sfrmpath = ages_path(/projects)+'sfrm/'
    paperpath = ages_path(/papers)+'sfrm/'

    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

; --------------------------------------------------
; stellar mass functions

    maxis1 = im_array(8.0,12.5,0.01)
    xrange = [7.6,12.0]
    yrange = [5E-7,0.1]

    subsample = ['all','quiescent','active']
    fit_color = ['black','red','navy blue']
    mf_color = ['grey','tomato','dodger blue']
    psym = [9,6,4]
    line = [0,0,0]
    
    psfile = paperpath+'sdss_mf'+suffix
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], width=6.8

    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      /ylog, yrange=yrange, xrange=xrange, xtitle=sfrm_masstitle(), $
      ytitle=sfrm_phititle()
    im_legend, ['All','Quiescent','Star-Forming'], /left, /bottom, box=0, $
      color=mf_color, psym=psym, line=[0,0,0], thick=8.0, symsize=2.0, $
      symthick=8
    legend, 'Baldry+08', /right, /top, box=0, line=3, pspacing=1.5, $
      color=djs_icolor('dark green'), charsize=1.6
    
    for jj = 0, n_elements(subsample)-1 do begin
       mf_data = mrdfits(sfrmpath+'sdss_mf_data_'+subsample[jj]+'.fits.gz',1,/silent)
       mf_fit = mrdfits(sfrmpath+'sdss_mf_fit_'+subsample[jj]+'.fits.gz',1,/silent)
       mf_fit_double = mrdfits(sfrmpath+'sdss_mf_fit_'+subsample[jj]+'_double.fits.gz',1,/silent)
       
       these = where((mf_data.phi gt 0.0),nthese)
       fullbin = mf_data.fullbin[these]
       mass = mf_data.mass[these]
       phi = mf_data.phi[these]
       phierr = mf_data.phierr[these]
       phierr_cv = mf_data.phierr_cv[these]
       phierr_model = mf_data.phierr_model[these]
       
       above = where(fullbin,comp=below)
       oploterror, mass[above], phi[above], phierr[above], $
         psym=symcat(psym[jj],thick=8), symsize=1.6, $
         color=fsc_color(mf_color[jj],100+jj), errcolor=fsc_color(mf_color[jj],100+jj)
;      oploterror, mass[below], phi[below], phierr[below], $
;        psym=symcat(psym[jj],thick=!p.thick), symsize=1.2, color=djs_icolor('grey'), $
;        errcolor=djs_icolor('grey')

; overplot the Schechter fits
;      if strmatch(subsample[jj],'**') eq 0
       djs_oplot, maxis1, mf_schechter_plus(10^maxis1,mf_fit_double), $
         line=line[jj], color=fsc_color(fit_color[jj],10+jj), thick=8.0
       oplot_bgd08_mf, maxis1, color='dark green', line=3
    endfor

    im_plotconfig, /psclose
    
return
end
    
