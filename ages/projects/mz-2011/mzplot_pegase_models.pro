pro mzplot_pegase_models, ps=ps
; jm11oct10ucsd - Pegase modeling results
    
    mzpath = mz_path()
    qapath = mzpath+'qaplots/'
    if keyword_set(ps) then begin
       pspath = qapath
       suffix = '.ps'
    endif else begin
       pspath = mz_path(/paper)
       suffix = '.eps'
    endelse

; -------------------------
    psfile = paperpath+'mzpegase_'
    im_plotconfig, 1, pos, psfile=psfile

; dependence of the MZ relation on alpha
    mzavg = mrdfits(mzpath+'mzevol_avg.fits.gz',1,/silent)
    massaxis = range(8.8,11.0,500)
    
; build a fiducial galaxy formation model and then make the QAplots 
    model = mz_galform_model(peginfo,alpha=alpha,mass0=mass0,tau0=tau0)
    nmodel = n_elements(model)
    

    ohrange1 = [8.0,9.5]
    massrange1 = [8.5,11.5]
    
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

    im_plotconfig, psfile=psfile, /psclose

return
end
