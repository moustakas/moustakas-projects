pro uband_junk

; JUNKY attempts at things    
    
    rootpath = getenv('CATALOGS_DIR')+'/sfhgrid/' ; parameter I/O path
    modelspath = rootpath+'basemodels/'           ; base models output path

    Zstr = 'm62'
    imfstr = 'salp'
    
    sfhfile = 'sfh_base_models_'+Zstr+'_'+imfstr+'.fits.gz'
    sfh = mrdfits(rootpath+sfhfile,1,/silent)

; ---------------------------------------------------------------------------    
; 
; ---------------------------------------------------------------------------    

; read the continuous SFH models    
    
    tau = [2,3,5,8,20]
    ntau = n_elements(tau)

    match, tau, sfh.tau, matchindx, contindx

    splog, 'Reading the continuous models.'
    for i = 0L, ntau-1L do begin
    
       cont1 = im_read_bc03(isedfile=sfh.contfile[contindx[i]],$
         isedpath=modelspath,/salpeter,bc03_extras=continfo1,$
         minwave=2500,maxwave=7000,/silent)
       if (i eq 0L) then begin
          cont = cont1
          continfo = continfo1
       endif else begin
          cont = [ [cont], [cont1] ]
          continfo = [ [continfo], [continfo1] ]
       endelse

    endfor
    cont = reform(cont)

    get_element, continfo[*,0].logage, alog10(age_universe), ageindx
    
    LU = const*10^(-0.4*continfo[ageindx,*].Umag) ; [erg/s]
    LHa = continfo[ageindx,*].sfr_yr/7.9D-42      ; [erg/s]
    D4000 = continfo[ageindx,*].b4_vn

; read the 500 Myr burst model

    splog, 'Reading the 500 Myr burst model.'
    burst = im_read_bc03(isedfile='burst_m62_salp_500Myr.ised',$
      isedpath=datapath,/salpeter,bc03_extras=burstinfo,$
      minwave=2500,maxwave=7000,/silent)

    burstages = [0.5,1.0,2.0,3.0,4.0,5.0]*1E9 ; [yr]
    get_element, burstinfo.logage, alog10(burstages), burstindx
    
    LU_burst = const*10^(-0.4*burstinfo[burstindx].Umag) ; [erg/s]
    LHa_burst = burstinfo[burstindx].sfr_yr/7.9D-42      ; [erg/s]
    D4000_burst = burstinfo[burstindx].b4_vn

; generate the plot    

    xrange = [0.5,2.1]
    yrange = [0.1,3.9]

    xtitle = 'D_{n}(4000)'
    ytitle = 'log [L(U)/L(H\alpha)] '
    
    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
      xsty=3, ysty=3, charsize=postthick, charthick=postthick, $
      xthick=postthick, ythick=postthick, xtitle=xtitle, ytitle=ytitle

    plotsym, 0, 2, /fill
    djs_oplot, D4000, alog10(LU/LHa), ps=-8, thick=postthick

    plotsym, 8, 2, /fill
    djs_oplot, D4000_burst, alog10(LU_burst/LHa_burst), ps=-8, thick=postthick
    
    stop

       
; ---------------------------------------------------------------------------    

; ---------------------------------------------------------------------------    
; consider an underlying stellar population with tau=4 Gyr (~Sbc
; galaxy; Bruzual & Charlot 1993); show the time evolution of the
; L(U)/L(Ha) ratio for this model, and for the superposition of a 10% 
; bursts lasting 500 Myr
; ---------------------------------------------------------------------------    

; run galaxev_csp, input bc2003_hr_m62_salp_ssp.ised, select "2",
; single burst of finite length, select 0.5 Gyr as the burst duration,
; and finally output filename = burst_m62_salp_500Myr.ised

; run add_bursts, using cont_m62_salp_008.ised and
; burst_m62_salp_500Myr.ised as the inputs; Burst1: beginning time = 0
; Gyr, amplitude = 1; Burst2: beginning time = 8 Gyr, amplitude = 0.1

;   burst = im_read_bc03(isedfile='burst_m62_salp_500Myr.ised',$
;     isedpath=modelspath,/salpeter,bc03_extras=eburst,/silent)

    cont = im_read_bc03(isedfile='cont_m62_salp_008.ised',$
      isedpath=modelspath,/salpeter,bc03_extras=cinfo,/silent)
    junk = im_read_bc03(isedfile='junk.ised',$
      isedpath='./',/salpeter,bc03_extras=info,/silent)

    ckeep = where(cinfo.logage lt alog10(age_universe))
    ctime = 10^cinfo[ckeep].logage
    cLU = const*10^(-0.4*cinfo[ckeep].Umag) ; [erg/s]
    cLHa = cinfo[ckeep].sfr_yr/7.9D-42      ; [erg/s]
    
    time = 10^info.logage
    LU = const*10^(-0.4*info.Umag) ; [erg/s]
    LHa = info.sfr_yr/7.9D-42      ; [erg/s]

    djs_plot, [0], [0], /nodata, xrange=[0,6], yrange=[0.5,2.5], $
      xsty=3, ysty=3, charsize=2.0, charthick=2.0, xthick=2.0, /xlog, $
      ythick=2.0, xtitle='Time [Gyr Ago]', ytitle='log [L(U)/L(H\alpha)]'
    djs_oplot, age_universe-ctime, alog10(cLU/cLHa), thick=2
    djs_oplot, time, alog10(LU/LHa), thick=2, line=2

stop    

; ---------------------------------------------------------------------------    
; this plot shows the time-evolution of the U-band luminosity    
; ---------------------------------------------------------------------------    

    djs_plot, [0], [0], /nodata, xrange=[6.5,10.2], yrange=[16,-3], $
      xsty=3, ysty=3, charsize=2.0, charthick=2.0, xthick=2.0, $
      ythick=2.0, xtitle='log age [yr]', ytitle='U'
    
    tau = sfh.tau
    ntau = n_elements(tau)

    isedfile = [sfh.contfile[[0,40]],sfh.instfile[9]]
    nised = n_elements(isedfile)
    
    for i = 0L, nised-1L do begin
    
       ssp = im_read_bc03(isedfile=isedfile[i],isedpath=modelspath,$
         /salpeter,bc03_extras=bce,/silent)

       LU = Uinfo.weff*Uinfo.vega_flam*10^(-0.4*bce.Umag) ; [erg/s/cm2]
       djs_oplot, bce.logage, bce.Umag
;      cc = get_kbrd(1)

    endfor

stop
    
; ---------------------------------------------------------------------------    

    minage = 0.01 ; [yr]
    maxage = 13.0 ; [yr]

    djs_plot, [0], [0], /nodata, xrange=[minage,maxage], yrange=[0.1,3.5], $
      xsty=3, ysty=3, charsize=2.0, charthick=2.0, xthick=2.0, $
      ythick=2.0, /xlog, xtitle='Time [Gyr]', ytitle='log [L(U)/L(H\alpha)]'
    
    tau = sfh.tau
    ntau = n_elements(tau)

    for i = 1L, ntau-1L do begin
    
       ssp = im_read_bc03(isedfile=sfh.contfile[i],isedpath=modelspath,$
         /salpeter,bc03_extras=bce,/silent)

       get_element, ssp.age/1E9, [minage,maxage], xx

       time = 10^bce[xx[0]:xx[1]].logage/1E9

       const = Uinfo.weff*Uinfo.vega_flam*(4*!dpi*vegadist^2)

       LU = const*10^(-0.4*bce[xx[0]:xx[1]].Umag) ; [erg/s]
       LHa = bce[xx[0]:xx[1]].sfr_yr/7.9D-42      ; [erg/s]

       djs_oplot, time, alog10(LU/LHa)
;      cc = get_kbrd(1)

    endfor

stop
    
; ---------------------------------------------------------------------------    
; what is the time evolution of the L(U)/L(Ha) ratio for continuous
; star formation, truncated after 30-300 Myr, viewed at ages ranging
; from 10-300 Myr?

    tburst = sfh.tburst
    nburst = n_elements(tburst)

;   for i = 9L, nburst-1L do begin
    for i = 0L, nburst-1L do begin
    
       ssp = im_read_bc03(isedfile=sfh.instfile[i],isedpath=modelspath,$
         /salpeter,bc03_extras=bce,/silent)
       maxage = 0.3 < tburst[i] ; [Gyr]

       get_element, ssp.age/1E9, [minage,maxage], xx

       time = 10^bce[xx[0]:xx[1]].logage/1E6

       LU = Uinfo.vega_flam*10^(-0.4*bce[xx[0]:xx[1]].Umag)      ; [erg/s/cm2]
       LHa = bce[xx[0]:xx[1]].sfr_yr/7.9D-42/(4*!dpi*vegadist^2) ; [erg/s/cm2]

       djs_oplot, time, alog10(LU/LHa)

    endfor
    

return
end
    
