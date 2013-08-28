pro qaplot_mz_pegase_models
; jm10oct20ucsd - build some simple QAplots of the Pegase models 
    
    mzpath = ages_path(/projects)+'mz/'
    peginfo = mzpegase_read_models()

    mzavg = mrdfits(mzpath+'mzevol_avg.fits.gz',1,/silent)
    massaxis = range(8.8,11.0,500)
    
    ohrange1 = [8.0,9.5]
    massrange1 = [8.5,11.5]
    
; build a fiducial galaxy formation model and then make the QAplots 
    model = mz_galform_model(peginfo,alpha=alpha,mass0=mass0,tau0=tau0)
    nmodel = n_elements(model)
    
    psfile = mzpath+'qaplots/qa_pegase_models.ps'
    im_plotconfig, 0, pos, psfile=psfile

; -------------------------
; dependence of the MZ relation on alpha
    zobs = 0.1 & zform = 2.0

    tau0 = 14.0 & mass0 = 10D^10.2
    alpha = range(-0.2,-0.8,3)

    plot, [0], [0], /nodata, xrange=massrange1, yrange=ohrange1, $
      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle=mzplot_ohtitle(), $
      position=pos
    legend, textoidl('\alpha='+string(alpha,format='(F4.1)')), /right, $
      /bottom, box=0, line=lindgen(n_elements(alpha)), pspacing=1.8
    legend, textoidl(['\tau='+string(tau0,format='(F4.1)')+' Gyr',$
      'log (M_{0}/M_{'+sunsymbol()+'})='+string(alog10(mass0),format='(F4.1)')]), $
      /left, /top, box=0
    aindx = findex(model[0].age,getage(zobs)-getage(zform))
    for ii = 0, n_elements(alpha)-1 do begin
       model1 = mz_galform_model(peginfo,alpha=alpha[ii],tau=tau0,mass0=mass0)
       log12oh = reform(interpolate(model1.log12oh,aindx,findgen(nmodel),/grid))
       mstar = reform(interpolate(model1.mstar,aindx,findgen(nmodel),/grid))
       djs_oplot, mstar, log12oh, psym=-6, line=ii
    endfor
    for jj = 0, 2 do djs_oplot, massaxis, mz_brokenpl(massaxis,mzlocal[jj].coeff_bin), color='red'
    djs_oplot, massaxis, tremonti_mz(massaxis=massaxis), color='blue', line=5

; -------------------------
; dependence of the MZ relation on tau0
    zobs = 0.1 & zform = 2.0

    alpha = -1.0 & mass0 = 10D^10.5
    tau0 = range(2.0,12.0,3)

    plot, [0], [0], /nodata, xrange=massrange1, yrange=ohrange1, $
      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle=mzplot_ohtitle(), $
      position=pos
    legend, textoidl('\tau='+string(tau0,format='(F4.1)')+' Gyr'), /right, $
      /bottom, box=0, line=lindgen(n_elements(tau0)), pspacing=1.8
    legend, textoidl(['\alpha='+string(alpha,format='(F4.1)'),$
      'log (M_{0}/M_{'+sunsymbol()+'})='+string(alog10(mass0),format='(F4.1)')]), $
      /left, /top, box=0
    aindx = findex(model[0].age,getage(zobs)-getage(zform))
    for ii = 0, n_elements(tau0)-1 do begin
       model1 = mz_galform_model(peginfo,alpha=alpha,tau=tau0[ii],mass0=mass0)
       log12oh = reform(interpolate(model1.log12oh,aindx,findgen(nmodel),/grid))
       mstar = reform(interpolate(model1.mstar,aindx,findgen(nmodel),/grid))
       djs_oplot, mstar, log12oh, psym=-6, line=ii
    endfor
    for jj = 0, 2 do djs_oplot, massaxis, mz_brokenpl(massaxis,mzlocal[jj].coeff_bin), color='red'
    djs_oplot, massaxis, tremonti_mz(massaxis=massaxis), color='blue', line=5

; -------------------------
; dependence of the MZ relation on mass0
    zobs = 0.1 & zform = 2.0

    alpha = -1.0 & tau0 = 6.0
    mass0 = 10^range(9.0,12.0,3.0)

    plot, [0], [0], /nodata, xrange=massrange1, yrange=ohrange1, $
      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle=mzplot_ohtitle(), $
      position=pos
    legend, textoidl('log (M_{0}/M_{'+sunsymbol()+'})='+string(alog10(mass0),format='(F4.1)')), /right, $
      /bottom, box=0, line=lindgen(n_elements(mass0)), pspacing=1.8
    legend, textoidl(['\alpha='+string(alpha,format='(F4.1)'),$
      '\tau='+string(tau0,format='(F4.1)')+' Gyr']), $
      /left, /top, box=0
    aindx = findex(model[0].age,getage(zobs)-getage(zform))
    for ii = 0, n_elements(mass0)-1 do begin
       model1 = mz_galform_model(peginfo,alpha=alpha,tau=tau0,mass0=mass0[ii])
       log12oh = reform(interpolate(model1.log12oh,aindx,findgen(nmodel),/grid))
       mstar = reform(interpolate(model1.mstar,aindx,findgen(nmodel),/grid))
       djs_oplot, mstar, log12oh, psym=-6, line=ii
    endfor
    for jj = 0, 2 do djs_oplot, massaxis, mz_brokenpl(massaxis,mzlocal[jj].coeff_bin), color='red'
    djs_oplot, massaxis, tremonti_mz(massaxis=massaxis), color='blue', line=5

; -------------------------
; dependence of the MZ relation on zform
    zobs = 0.1
    zform = im_array(1.0,3.0,0.5)
    plot, [0], [0], /nodata, xrange=massrange1, yrange=ohrange1, $
      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle=mzplot_ohtitle(), $
      position=pos
    legend, 'zform', /left, /top, box=0, margin=0
    for ii = 0, n_elements(zform)-1 do begin
       aindx = findex(model[0].age,getage(zobs)-getage(zform[ii]))
       log12oh = reform(interpolate(model.log12oh,aindx,findgen(nmodel),/grid))
       mstar = reform(interpolate(model.mstar,aindx,findgen(nmodel),/grid))
       djs_oplot, mstar, log12oh, psym=-6, line=ii

stop       
       
    endfor

    im_plotconfig, psfile=psfile, /psclose, /gzip

stop
    
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
    plot, [0], [0], /nodata, xrange=massrange1, yrange=ohrange1, $
      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle=mzplot_ohtitle(), $
      position=pos
    for jj = 0, n_elements(zobs)-1 do djs_oplot, newmass[jj,*], log12oh[jj,*], psym=-6
;   for ii = 0, ntau-1 do djs_oplot, newmass[*,ii], log12oh[*,ii], psym=-6
    legend, 'fixed zobs', /left, /top, box=0, margin=0

; at fixed tau    
    plot, [0], [0], /nodata, xrange=massrange1, yrange=ohrange1, $
      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle=mzplot_ohtitle(), $
      position=pos
    for ii = 0, ntau-1 do djs_oplot, newmass[*,ii], log12oh[*,ii], psym=-6
    legend, 'fixed tau', /left, /top, box=0, margin=0

; ssfr 
    
    
    im_plotconfig, psfile=psfile, /psclose, /gzip
    
    
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
