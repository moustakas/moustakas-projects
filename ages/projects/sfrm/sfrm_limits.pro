;+
; NAME:
;   SFRM_LIMITS
;
; PURPOSE:
;   Compute the limiting stellar mass as a function of redshift. 
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;    If I_lim=19.95 is the apparent magnitude limit of the survey, I
;    is the apparent magnitude of any source, and M_g is its absolute
;    magnitude (computed using K-correct), then the limiting absolute
;    g-band magnitude for that source is given simply by: 
;       M_g_lim = M_g - (I-I_lim) 
;
;    The limiting stellar mass, log_M_lim, is similarly given by:
;    log_M_lim = log_M + 0.4(I-I_lim) 
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

pro sfrm_limits, quiescent=quiescent, active=active

    sfrmpath = ages_path(/projects)+'sfrm/'
    parent = read_sfrm_sample()
    zbins = sfrm_zbins(nzbins)

    suffix = 'all'
    if keyword_set(quiescent) then suffix = 'quiescent'
    if keyword_set(active) then suffix = 'active'
    outfile = sfrmpath+'sfrm_limits_'+suffix+'.fits'

; specify the binning parameters    
    I_lim = 19.95
    jv2ab = k_vega2ab(filterlist='twomass_J.par',/kurucz,/silent)
    itag = tag_indx(parent,'i_tot')
    
    zbin = 0.02
    zaxis = im_array(0.05,0.75,zbin)
    nzaxis = n_elements(zaxis)
    quant = [0.5,0.75,0.95]

    limits = {zaxis: zaxis, $
      minmr_50: fltarr(nzaxis), minmr_75: fltarr(nzaxis), minmr_95: fltarr(nzaxis), $
      minmg_50: fltarr(nzaxis), minmg_75: fltarr(nzaxis), minmg_95: fltarr(nzaxis), $
      minmass_50: fltarr(nzaxis), minmass_75: fltarr(nzaxis), minmass_95: fltarr(nzaxis), $
      minmr_poly: fltarr(nzaxis), minmg_poly: fltarr(nzaxis), minmass_poly: fltarr(nzaxis), $
      zbin: zbins.zbin, minmr: fltarr(nzbins), minmg: fltarr(nzbins), minmass: fltarr(nzbins)}

; select the sample    
    nuvmr = parent.galex_absmag[1]-parent.ugriz_absmag[2]
    rmj = parent.ugriz_absmag[2]-(parent.ubvrijhk_absmag[5]+jv2ab)
    qq = select_quiescent(nuvmr,rmj,active=aa)
    if keyword_set(quiescent) then parent = parent[qq]
    if keyword_set(active) then parent = parent[aa]
    weight = parent.final_weight ; spectroscopic weight

; 0.1r absolute magnitude
    minmr = parent.ugriz_absmag[2] - (parent.(itag)-I_lim)
    for ii = 0, nzaxis-2 do begin
       these = where((parent.z ge zaxis[ii]) and (parent.z lt zaxis[ii+1]),nthese)
       weight1 = weight[these]
       minmr1 = minmr[these]
       stats = im_quantile(minmr1,weight1,quant=quant)
       limits.minmr_50[ii] = stats[0]
       limits.minmr_75[ii] = stats[1]
       limits.minmr_95[ii] = stats[2]
    endfor

; 0.1g absolute magnitude
    minmg = parent.ugriz_absmag[1] - (parent.(itag)-I_lim)
    for ii = 0, nzaxis-2 do begin
       these = where((parent.z ge zaxis[ii]) and (parent.z lt zaxis[ii+1]),nthese)
       weight1 = weight[these]
       minmg1 = minmg[these]
       stats = im_quantile(minmg1,weight1,quant=quant)
       limits.minmg_50[ii] = stats[0]
       limits.minmg_75[ii] = stats[1]
       limits.minmg_95[ii] = stats[2]
    endfor

; mass    
    masslim = parent.mass + 0.4*(parent.(itag)-I_lim)
    for ii = 0, nzaxis-2 do begin
       these = where((parent.z ge zaxis[ii]) and (parent.z lt zaxis[ii+1]),nthese)
       weight1 = weight[these]
       masslim1 = masslim[these]
       stats = im_quantile(masslim1,weight1,quant=quant)
       limits.minmass_50[ii] = stats[0]
       limits.minmass_75[ii] = stats[1]
       limits.minmass_95[ii] = stats[2]
    endfor
    
; fit a polynomial model to the 75% limiting magnitude and mass,
; excluding the first and last redshift bins
    good = where((limits.zaxis gt 0.05) and (limits.zaxis lt 0.75) and $
      (limits.minmr_75 ne 0.0),ngood)
    im_poly_iter, limits.zaxis[good], limits.minmr_75[good], 4, coeff=coeff
    limits.minmr_poly = poly(limits.zaxis,coeff)
    
    good = where((limits.zaxis gt 0.05) and (limits.zaxis lt 0.75) and $
      (limits.minmg_75 ne 0.0),ngood)
    im_poly_iter, limits.zaxis[good], limits.minmg_75[good], 4, coeff=coeff
    limits.minmg_poly = poly(limits.zaxis,coeff)
    
    good = where((limits.zaxis gt 0.05) and (limits.zaxis lt 0.75) and $
      (limits.minmass_75 ne 0.0),ngood)
    im_poly_iter, limits.zaxis[good], limits.minmass_75[good], 4, coeff=coeff
    limits.minmass_poly = poly(limits.zaxis,coeff)

; finally get the limiting magnitude in mass at the center of each of
; our six redshift bins by interpolating the polynomial fit
    limits.minmr = interpol(limits.minmr_poly,limits.zaxis,limits.zbin)
    limits.minmg = interpol(limits.minmg_poly,limits.zaxis,limits.zbin)
    limits.minmass = interpol(limits.minmass_poly,limits.zaxis,limits.zbin)
    limits.minmass = fix(limits.minmass*10.0)/10.0
    
;; compute the average 75% magnitude and stellar mass limits in each of
;; our six redshift bins
;    for ii = 0, n_elements(zbins)-1 do begin
;       these = where((limits.zaxis ge zbins[ii].zlo) and (limits.zaxis lt zbins[ii].zup))
;;      limits.minmg[ii] = im_min(limits.minmg_75[these],sigrej=3.0)
;;      limits.minmass[ii] = im_min(limits.minmass_75[these],sigrej=3.0)
;       limits.minmg[ii] = djs_median(limits.minmg_75[these])
;       limits.minmass[ii] = djs_median(limits.minmass_75[these])
;    endfor

; write out    
    im_mwrfits, limits, outfile, /clobber
    
; make a QAplot
    colors = ['red','green','blue']

    psfile = sfrmpath+'qaplots/sfrm_limits_'+suffix+'.ps'
    im_plotconfig, 0, psfile=psfile

    djs_plot, parent.z, parent.ugriz_absmag[2], psym=2, yrange=[-15,-24], $
      xtitle='Redshift', ytitle='M_{0.1r}', sym=0.3, xsty=1, ysty=1, $
      xrange=[0,0.8]
    im_legend, ['50%','75%','95%'], /right, /bottom, $
      box=0, color=colors, line=[0,0,0], pspacing=1.2
    djs_oplot, zaxis, limits.minmr_50, color=colors[0], line=0, psym=-8
    djs_oplot, zaxis, limits.minmr_75, color=colors[1], line=0, psym=-8
    djs_oplot, zaxis, limits.minmr_95, color=colors[2], line=0, psym=-8
    djs_oplot, limits.zbin, limits.minmr, psym=symcat(6,thick=8), color='orange', symsize=3.0
    djs_oplot, zaxis, limits.minmr_poly, line=0

    djs_plot, parent.z, parent.ugriz_absmag[1], psym=2, yrange=[-15,-24], $
      xtitle='Redshift', ytitle='M_{0.1g}', sym=0.3, xsty=1, ysty=1, $
      xrange=[0,0.8]
    im_legend, ['50%','75%','95%'], /right, /bottom, $
      box=0, color=colors, line=[0,0,0], pspacing=1.2
    djs_oplot, zaxis, limits.minmg_50, color=colors[0], line=0, psym=-8
    djs_oplot, zaxis, limits.minmg_75, color=colors[1], line=0, psym=-8
    djs_oplot, zaxis, limits.minmg_95, color=colors[2], line=0, psym=-8
    djs_oplot, limits.zbin, limits.minmg, psym=symcat(6,thick=8), color='orange', symsize=3.0
    djs_oplot, zaxis, limits.minmg_poly, line=0

    djs_plot, parent.z, parent.mass, psym=2, yrange=[8,12], $
      xtitle='Redshift', ytitle='log (M_{*}/M'+sunsymbol()+')', sym=0.3, $
      xsty=1, ysty=1, xrange=[0,0.8]
    im_legend, ['50%','75%','95%'], /right, /bottom, $
      box=0, color=colors, line=[0,0,0], pspacing=1.2
    djs_oplot, zaxis, limits.minmass_50, color=colors[0], line=0, psym=-8
    djs_oplot, zaxis, limits.minmass_75, color=colors[1], line=0, psym=-8
    djs_oplot, zaxis, limits.minmass_95, color=colors[2], line=0, psym=-8
    djs_oplot, limits.zbin, limits.minmass, psym=symcat(6,thick=8), color='orange', symsize=3.0
    djs_oplot, zaxis, limits.minmass_poly, line=0

    im_plotconfig, psfile=psfile, /psclose, /gzip
    
return
end
    
