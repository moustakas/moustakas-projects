pro mz_limits, quiescent=quiescent, active=active
; jm11nov02ucsd - estimate the limiting stellar mass and absolute
; magnitude as a function of redshift
    
;    If I_lim=20 is the apparent magnitude limit of the survey, I is
;    the apparent magnitude of any source, and M_B is its absolute
;    magnitude (computed using K-correct), then the limiting absolute
;    B-band magnitude for that source is given simply by: M_B_lim =
;    M_B - (I-I_lim)
;
;    The limiting stellar mass, log_M_lim, is similarly given by:
;    log_M_lim = log_M + 0.4(I-I_lim) 
    
    mzpath = mz_path()

    parent = read_mz_sample(/mzhii_ancillary)
    mass = read_mz_sample(/mzhii_mass)
    zbins = mz_zbins(nzbins)

    outfile = mzpath+'mz_limits.fits'

; specify the binning parameters    
    ifaint = mz_ifaint(select_filter=ifilter,/vega)
    itag = tag_indx(parent,'i_tot')
    quant = [0.5,0.75,0.95]

    zbinsize = 0.04
    zmin = min(zbins.zlo);+zbinsize/2.0
    zmax = max(zbins.zup);-zbinsize/2.0
    nzaxis = (zmax-zmin)/zbinsize
    zaxis = range(zmin+zbinsize/2.0,zmax-zbinsize/2.0,nzaxis)

    ncoeff = 3
    limits = {zaxis: zaxis, $
      mblimit_50_zaxis: fltarr(nzaxis), mblimit_75_zaxis: fltarr(nzaxis), mblimit_95_zaxis: fltarr(nzaxis), $
      masslimit_50_zaxis: fltarr(nzaxis), masslimit_75_zaxis: fltarr(nzaxis), masslimit_95_zaxis: fltarr(nzaxis), $
      mblimit_50_coeff: fltarr(ncoeff), mblimit_75_coeff: fltarr(ncoeff), mblimit_95_coeff: fltarr(ncoeff), $
      masslimit_50_coeff: fltarr(ncoeff), masslimit_75_coeff: fltarr(ncoeff), masslimit_95_coeff: fltarr(ncoeff), $
; the numbers we want:
      zbin: zbins.zbin, $
      mblimit_50: fltarr(nzbins), mblimit_75: fltarr(nzbins), mblimit_95: fltarr(nzbins), $
      masslimit_50: fltarr(nzbins), masslimit_75: fltarr(nzbins), masslimit_95: fltarr(nzbins)}

    weight = parent.final_weight ; spectroscopic weight

; B-band magnitude
    mblimit = parent.k_ubvrijhk_absmag_00[1] - (parent.(itag)-ifaint)
    for ii = 0, nzaxis-1 do begin
       these = where((parent.z ge zaxis[ii]-zbinsize/2) and (parent.z lt zaxis[ii]+zbinsize/2),nthese)
       weight1 = weight[these]
       mblimit1 = mblimit[these]
       stats = weighted_quantile(mblimit1,weight1,quant=quant)
       limits.mblimit_50_zaxis[ii] = stats[0]
       limits.mblimit_75_zaxis[ii] = stats[1]
       limits.mblimit_95_zaxis[ii] = stats[2]
    endfor

; mass    
    masslim = mass.mass_50 + 0.4*(parent.(itag)-ifaint)
    for ii = 0, nzaxis-1 do begin
       these = where((parent.z ge zaxis[ii]-zbinsize/2) and (parent.z lt zaxis[ii]+zbinsize/2),nthese)
       weight1 = weight[these]
       masslim1 = masslim[these]
       stats = weighted_quantile(masslim1,weight1,quant=quant)
       limits.masslimit_50_zaxis[ii] = stats[0]
       limits.masslimit_75_zaxis[ii] = stats[1]
       limits.masslimit_95_zaxis[ii] = stats[2]
    endfor
    
; fit a polynomial model to the limiting magnitude and mass, excluding
; the first and last redshift bins
    these = where(zaxis lt 0.7)

    coeff_50 = robust_poly_fit(zaxis[these],limits.mblimit_50_zaxis[these],ncoeff-1)
    coeff_75 = robust_poly_fit(zaxis[these],limits.mblimit_75_zaxis[these],ncoeff-1)
    coeff_95 = robust_poly_fit(zaxis[these],limits.mblimit_95_zaxis[these],ncoeff-1)
;   im_poly_iter, zaxis[these], limits.mblimit_50_zaxis[these], ncoeff-1, coeff=coeff_50
;   im_poly_iter, zaxis[these], limits.mblimit_75_zaxis[these], ncoeff-1, coeff=coeff_75
;   im_poly_iter, zaxis[these], limits.mblimit_95_zaxis[these], ncoeff-1, coeff=coeff_95
    limits.mblimit_50_coeff = coeff_50
    limits.mblimit_75_coeff = coeff_75
    limits.mblimit_95_coeff = coeff_95
    
    coeff_50 = robust_poly_fit(zaxis[these],limits.masslimit_50_zaxis[these],ncoeff-1)
    coeff_75 = robust_poly_fit(zaxis[these],limits.masslimit_75_zaxis[these],ncoeff-1)
    coeff_95 = robust_poly_fit(zaxis[these],limits.masslimit_95_zaxis[these],ncoeff-1)
;   im_poly_iter, zaxis[these], limits.masslimit_50_zaxis[these], ncoeff-1, coeff=coeff_50
;   im_poly_iter, zaxis[these], limits.masslimit_75_zaxis[these], ncoeff-1, coeff=coeff_75
;   im_poly_iter, zaxis[these], limits.masslimit_95_zaxis[these], ncoeff-1, coeff=coeff_95
    limits.masslimit_50_coeff = coeff_50
    limits.masslimit_75_coeff = coeff_75
    limits.masslimit_95_coeff = coeff_95
    
; finally get the limiting magnitude in mass at the center of each of
; our six redshift bins by interpolating the polynomial fit
    limits.mblimit_50 = interpol(poly(zaxis,limits.mblimit_50_coeff),zaxis,limits.zbin)
    limits.mblimit_75 = interpol(poly(zaxis,limits.mblimit_75_coeff),zaxis,limits.zbin)
    limits.mblimit_95 = interpol(poly(zaxis,limits.mblimit_95_coeff),zaxis,limits.zbin)

    limits.masslimit_50 = interpol(poly(zaxis,limits.masslimit_50_coeff),zaxis,limits.zbin)
    limits.masslimit_75 = interpol(poly(zaxis,limits.masslimit_75_coeff),zaxis,limits.zbin)
    limits.masslimit_95 = interpol(poly(zaxis,limits.masslimit_95_coeff),zaxis,limits.zbin)

;; compute the average 75% magnitude and stellar mass limits in each of
;; our six redshift bins
;    for ii = 0, n_elements(zbins)-1 do begin
;       these = where((limits.zaxis ge zbins[ii].zlo) and (limits.zaxis lt zbins[ii].zup))
;;      limits.mblimit[ii] = im_min(limits.mblimit_75[these],sigrej=3.0)
;;      limits.masslimit[ii] = im_min(limits.masslimit_75[these],sigrej=3.0)
;       limits.mblimit[ii] = djs_median(limits.mblimit_75[these])
;       limits.masslimit[ii] = djs_median(limits.masslimit_75[these])
;    endfor

; write out    
    im_mwrfits, limits, outfile, /clobber

    niceprint, limits.zbin, limits.masslimit_50, limits.masslimit_75, limits.masslimit_95, $
      limits.mblimit_50, limits.mblimit_75, limits.mblimit_95, $
      limits.masslimit_75-limits.masslimit_50, limits.masslimit_95-limits.masslimit_50, $
      limits.mblimit_75-limits.mblimit_50, limits.mblimit_95-limits.mblimit_50
    splog, mean(limits.masslimit_75-limits.masslimit_50), mean(limits.masslimit_95-limits.masslimit_50), $
      mean(limits.mblimit_75-limits.mblimit_50), mean(limits.mblimit_95-limits.mblimit_50)
    
; make a QAplot
    colors = ['red','green','blue']

    psfile = mzpath+'qaplots/mz_limits.ps'
    im_plotconfig, 6, pos, psfile=psfile

    djs_plot, parent.z, parent.k_ubvrijhk_absmag_00[1], psym=2, yrange=[-15,-24], $
      xtitle='', ytitle=mzplot_mbtitle(), sym=0.3, xsty=1, ysty=1, $
      xrange=[0,0.8], position=pos[*,0], xtickname=replicate(' ',10)
    im_legend, ['50%','75%','95%'], /right, /bottom, $
      box=0, color=colors, line=[0,0,0], pspacing=1.2
    djs_oplot, zaxis, limits.mblimit_50_zaxis, color=colors[0], line=0, psym=-8
    djs_oplot, zaxis, limits.mblimit_75_zaxis, color=colors[1], line=0, psym=-8
    djs_oplot, zaxis, limits.mblimit_95_zaxis, color=colors[2], line=0, psym=-8
    djs_oplot, limits.zbin, limits.mblimit_50, psym=symcat(6,thick=8), color='orange', symsize=3.0
    djs_oplot, zaxis, poly(zaxis,limits.mblimit_50_coeff), line=0, thick=8

    djs_plot, parent.z, mass.mass_50, psym=2, yrange=[8,12], $
      xtitle='Redshift', ytitle=mzplot_masstitle(), sym=0.3, xsty=1, ysty=1, $
      xrange=[0,0.8], position=pos[*,1], /noerase
    im_legend, ['50%','75%','95%'], /right, /bottom, $
      box=0, color=colors, line=[0,0,0], pspacing=1.2
    djs_oplot, zaxis, limits.masslimit_50_zaxis, color=colors[0], line=0, psym=-8
    djs_oplot, zaxis, limits.masslimit_75_zaxis, color=colors[1], line=0, psym=-8
    djs_oplot, zaxis, limits.masslimit_95_zaxis, color=colors[2], line=0, psym=-8
    djs_oplot, limits.zbin, limits.masslimit_50, psym=symcat(6,thick=8), color='orange', symsize=3.0
    djs_oplot, zaxis, poly(zaxis,limits.masslimit_50_coeff), line=0, thick=8

    im_plotconfig, psfile=psfile, /psclose
    
return
end
    
