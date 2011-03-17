pro mz_pegase_models
; jm10oct20ucsd - build some simple toy models to interpret our
; results using Pegase
    
    mzpath = ages_path(/projects)+'mz/'
    rawinfo = get_pegase_info()
    info = mz_galform_model(rawinfo,alpha=alpha,mass0=mass0,tau0=tau0)
    
stop    
    
; given the stellar mass grid, assume tau is a power-law function of
; mass, and then interpolate the models appropriately
    mass = range(8.0,11.5,20) 
    alpha = -0.9
    mchar = 10.0
    tau = 10^(alpha*mass + mchar)
;   tau = (10^(mass-mchar))^alpha
    keep = where((tau ge min(info.tau)) and (tau le max(info.tau)),ntau)
    mass = mass[keep]
    tau = tau[keep]
    print, tau
    tindx = findex(info.tau,tau)

; ##################################################
; build some QAplots
    psfile = mzpath+'qaplots/qa_pegase_models.ps'
    im_plotconfig, 0, pos, psfile=psfile

; -------------------------
; dependence on zform
    zobs = 0.1
    zform = im_array(1.0,3.0,0.5)
    plot, [0], [0], /nodata, xrange=[8.5,11.5], yrange=[7.5,9.5], $
      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle=mzplot_ohtitle(), $
      position=pos
    legend, 'zform', /left, /top, box=0, margin=0
    for ii = 0, n_elements(zform)-1 do begin
       aindx = findex(info[0].age,getage(zobs)-getage(zform[ii]))
       log12oh = reform(interpolate(info.log12oh,aindx,tindx,/grid,missing=-1.0))
       djs_oplot, mass, log12oh, psym=-6, line=ii
    endfor

; -------------------------
; dependence on the observed redshift
    zobs = [0.1,0.3,0.5,0.7,1.0,1.5]
    zform = 2.0

    aindx = findex(info[0].age,getage(zobs)-getage(zform))
    log12oh = reform(interpolate(info.log12oh,aindx,tindx,/grid,missing=-1.0))
    mstar = reform(interpolate(info.mstar,aindx,tindx,/grid,missing=-1.0))

    newmass = log12oh*0.0 ; adjust the stellar mass of the models to the final desired mass     
    for ii = 0, ntau-1 do newmass[*,ii] = mass[ii]+alog10(mstar[*,ii]/mstar[0,ii])

; at fixed zobs    
    plot, [0], [0], /nodata, xrange=[8.5,11.5], yrange=[7.5,9.5], $
      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle=mzplot_ohtitle(), $
      position=pos
    for jj = 0, n_elements(zobs)-1 do djs_oplot, newmass[jj,*], log12oh[jj,*], psym=-6
;   for ii = 0, ntau-1 do djs_oplot, newmass[*,ii], log12oh[*,ii], psym=-6
    legend, 'fixed zobs', /left, /top, box=0, margin=0

; at fixed tau    
    plot, [0], [0], /nodata, xrange=[8.5,11.5], yrange=[7.5,9.5], $
      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle=mzplot_ohtitle(), $
      position=pos
    for ii = 0, ntau-1 do djs_oplot, newmass[*,ii], log12oh[*,ii], psym=-6
    legend, 'fixed tau', /left, /top, box=0, margin=0

; ssfr 
    
    
    im_plotconfig, psfile=psfile, /psclose, /gzip

    
stop
    ssfr = reform(interpolate(info.zgas,aindx,tindx,/grid))


    zgas = reform(interpolate(info.zgas,aindx,tindx,/grid))
    plot, mass, zgas, psym=-6, ysty=3, xsty=3;, yrange=[7,9.5]

    
    
    
stop
stop
    
; tests with Noeske's models
    
;   zf = 10^(-2.7)*(10^mass)^0.3-1 ; eq 6

    ngal = 20 & nz = 30
    mass = range(9.5,12.0,ngal)
    zf = 2.0
    tau = 10^20.7/10^mass/1D9      ; [Gyr]
    age = getage(0.3)-getage(zf) ; observed at z=0.3

    recy = 0.5
    
;   sfr_zf = 1.0/(1D9*tau*(1-recy)) ; sfr at z=zf
    sfr_zf = 10^mass/(1D9*tau*(1-recy)) ; sfr at z=zf
    sfr = alog10(sfr_zf*exp(-age/tau)) ; sfr at z=0.3
    ssfr = sfr - mass

    plot, mass, ssfr, yrange=[-11,-8], ysty=3
    
    zaxis = range(0.05,1.0,nz)

    age = fltarr(nz,ngal)
    for ii = 0, ngal-1 do age[*,ii] = getage(zaxis)-getage(zf[ii])
    age = getage(zf)

return
end
    
;;
;;    zf = [1.0,2.0,3.0] ; formation redshift
;;    nzf = n_elements(zf)
;;
;;    nage = 100
;;    info = {$
;;      tau:              0.0,$
;;      zform:            0.0,$
;;      age:     fltarr(nage),$
;;      log12oh: fltarr(nage),$
;;      mstar:   fltarr(nage),$
;;      sfr:     fltarr(nage)}
;;    info = replicate(info,ntau,nzf)
;;    info.tau = tau#(fltarr(nzf)+1)
;;    info.zform = transpose(zf#(fltarr(ntau)+1))
;;    
;;    for ii = 0, ntau-1 do begin
;;       peg = mrdfits(sfhfile[ii],1,/silent)
;;       for jj = 0, nzf-1 do begin
;;          age = range(0.1,getage(0.0)-getage(zf[jj]),nage) ; [Gyr]
;;          findx = findex(peg.age/1D9,age)
;;          info[ii,jj].age = age
;;          info[ii,jj].log12oh = alog10(interpolate(peg.zgas,findx)/zsun)+log12ohsun
;;          info[ii,jj].mstar = interpolate(peg.mstar,findx)
;;          info[ii,jj].sfr = interpolate(peg.sfr,findx)
;;       endfor
;;    endfor
