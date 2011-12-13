function mzpegase_func, mstar, params, zobs=zobs, normmass=normmass
; params[0] = alpha
; params[1] = beta
; params[2] = tau0
; params[3] = zform0
    model = mzpegase_galform_model(zobs,alpha=params[0],$
      beta=params[1],tau0=params[2],zform0=params[3])
    log12oh = interpol(model.log12oh,model.mstar,mstar)
    log12oh = log12oh-interpol(log12oh,mstar,normmass)
return, log12oh
end


pro mzpegase_plots, ps=ps, optimize=optimize
; jm11nov11ucsd - Pegase modeling for the discussion section
    
    mzpath = mz_path()
    qapath = mzpath+'qaplots/'
    if keyword_set(ps) then begin
       pspath = qapath
       suffix = '.ps'
    endif else begin
       pspath = mz_path(/paper)
       suffix = '.eps'
    endelse

; read the measured MZ relation and construct a "mean" MZ relation
; among the three calibrations normalized at 10^10.3 Msun
    zbins = mz_zbins(nz)

    mzavg = mrdfits(mzpath+'mzevol_avg.fits.gz',1,/silent)
    calib = strtrim(mzavg.calib,2)
    calibcolor = ['firebrick','forest green','dodger blue']
    calibpsym = [16,17,15]
    calibsymsize = [1.5,1.8,1.5]*0.8
    calibline = [0,3,5]

    normmass = 10.5D
    nmass = 50
    mzobsmass = range(9.0,11.2,nmass)
    mzobs = mzobsmass*0.0
    mzerr = mzobsmass*0.0
    zobs = 0.1

    ohmodel = fltarr(nmass,3)
    for jj = 0, 2 do begin
       ohmodel[*,jj] = mzevol_func(mzobsmass,mzavg.mzevol_coeffs_r0zero[*,jj],z=zobs,qz0=mzavg.qz0)
       ohmodel[*,jj] = ohmodel[*,jj]-interpol(ohmodel[*,jj],mzobsmass,normmass)
    endfor
    for kk = 0, nmass-1 do begin
       mzobs[kk] = djs_mean(ohmodel[kk,*])
;      mzerr[kk] = djsig(ohmodel[kk,*])
       mzerr[kk] = (max(ohmodel[kk,*])-min(ohmodel[kk,*])) ;/2.0
    endfor

; ---------------------------------------------------------------------------

    maxis = range(9.0,11.2,50)
    zz = [0.1,0.3,0.5,0.7]
    nzz = n_elements(zz)
    yield = 0.01
    
    djs_plot, [0], [0], /nodata, xrange=[9,11.4], yrange=[8.0,9.5]
    for ii = 0, nzz-1 do begin
       eta = (10^maxis/10D^10.0)^(-2.0/3.0)
       alpha = (0.5-0.1*zz[ii])*(10^maxis/10D^10.0)^0.25
       zism = yield/(1+eta)*(1/(1-alpha))
       djs_oplot, maxis, alog10(zism/0.02)+8.7, line=ii
    endfor
    

stop    
    
; ---------------------------------------------------------------------------
; optimize the parameters of the model
    if keyword_set(optimize) then begin
       parinfo = replicate({value: 1D, limited: [0,0], limits: [0D,0D], fixed: 0},4)
       parinfo.value = [-1.0D,+0.3D,5D,4D]
;      parinfo[0].limited[1]  = 1 & parinfo[0].limits[1] = -0.01 ; alpha
;      parinfo[1].limited[1]  = 1 & parinfo[1].limits[1] = 2.0 ; 0.8         ; beta
;      parinfo[2].limited = 1 & parinfo[2].limits = 10D^[9.5,11.5] ; mass0
       
       functargs = {zobs: zobs, normmass: normmass}
       params = mpfitfun('mzpegase_func',mzobsmass,mzobs,mzerr,$
         parinfo=parinfo,functargs=functargs,perror=perror,dof=dof,$
         covar=covar,status=mpstatus,quiet=quiet,bestnorm=chi2,$
         yfit=yfit)
       djs_plot, mzobsmass, mzobs, xsty=3, ysty=3
       djs_oplot, mzobsmass, mzobs+mzerr/2.0, line=5
       djs_oplot, mzobsmass, mzobs-mzerr/2.0, line=5
       djs_oplot, mzobsmass, yfit, line=0, color='blue'
       return
    endif
    
; ---------------------------------------------------------------------------
; Figure 23 - dlogoh and SFR/M vs redshift for three stellar masses 
    psfile = pspath+'mzpegase_dlogoh_sfrm'+suffix
    im_plotconfig, 6, pos, psfile=psfile, xmargin=[1.5,0.4], width=6.6, $
      height=[3.5,3.5]

    xrange1 = [-0.02,0.8]
    yrange1 = [8.4,9.25]
;   yrange1 = [-0.18,+0.05]
    yrange2 = [-1.8,-0.05]

    xtitle1 = 'Redshift'
    ytitle1 = mzplot_ohtitle()
;   ytitle1 = '\Delta'+' [12 + log (O/H)]'
    ytitle2 = textoidl('log (sSFR / Gyr^{-1})') ; mzplot_sfrmtitle()

    mstar = [10.0,10.5,11.0]
    mstarlabel = ['log(M/M'+sunsymbol()+')=','','']+string(mstar,format='(F4.1)')
    mstarcolor = ['navy','tan','black']
    modelcolor = ['navy','tan','black']
    mstarpsym = [16,17,15]
    mstarsymsize1 = [1.2,1.2,1.2]*2.0
    mstarsymsize2 = [1.2,1.2,1.2]*1.8
    mstarline = [0,5,3]
    nmstar = n_elements(mstar)

    nz = 20
    zaxis = range(0.02,xrange1[1]-0.02,nz)
    zaxis2 = range(0.1,0.7,7)

; build the pegase models    
    zform0 = 3.0
    tau0 = 8.0
    alpha = 1.0
    beta = +0.8

    modelsfrm = fltarr(nmstar,nz)
    modeloh = fltarr(nmstar,nz)
;   refmodel = mzpegase_galform_model(0.1,alpha=alpha,beta=beta,$
;     zform0=zform0,tau0=tau0)
    for iz = 0, nz-1 do begin
       model = mzpegase_galform_model(zaxis[iz],alpha=alpha,beta=beta,$
         zform0=zform0,tau0=tau0)
       modelsfrm[*,iz] = interpol(model.sfrm,model.mstar,mstar)
       modeloh[*,iz] = interpol(model.log12oh,model.mstar,mstar);-$
;        interpol(refmodel.log12oh,refmodel.mstar,mstar)
    endfor
    
; ####################
; metallicity evolution    
    djs_plot, [0], [0], /nodata, ysty=1, xsty=1, position=pos[*,0], $
      xrange=xrange1, yrange=yrange1, xtitle='', ytitle=ytitle1, $
      xtickname=replicate(' ',10), ytickinterval=0.2
    zoff = [0.0,0.0,+0.005]
    for mm = nmstar-1, 0, -1 do begin
       mzbest = mzevol_func(mstar[mm],mzavg.mzevol_coeffs_avg,qz0=mzavg.qz0,z=zaxis2)
       mzhi = mzevol_func(mstar[mm],mzavg.mzevol_coeffs_avg+mzavg.mzevol_coeffs_avg_err/2.0,qz0=mzavg.qz0,z=zaxis2)
       mzlo = mzevol_func(mstar[mm],mzavg.mzevol_coeffs_avg-mzavg.mzevol_coeffs_avg_err/2.0,qz0=mzavg.qz0,z=zaxis2)
       
       oploterror, zaxis2+zoff[mm], mzbest, mzhi-mzbest, color=im_color(mstarcolor[mm]), psym=symcat(mstarpsym[mm]), $
         symsize=mstarsymsize1[mm], errcolor=im_color(mstarcolor[mm]), /hibar
       oploterror, zaxis2+zoff[mm], mzbest, mzbest-mzlo, color=im_color(mstarcolor[mm]), $
         errcolor=im_color(mstarcolor[mm]), psym=3, /lobar

;      djs_oplot, zaxis2, mzevol_func(mstar[mm],mzavg.mzevol_coeffs_avg,qz0=mzavg.qz0,z=zaxis2), $
;        color=im_color(mstarcolor[mm]), psym=symcat(mstarpsym[mm]), $
;        symsize=mstarsymsize1[mm]*0.8; line=mstarline[mm], thick=8

;      slope = poly(mstar[mm]-mzavg.dlogohdz_normmass,mzavg.dlogohdz_coeff)
;      djs_oplot, zaxis2, poly(zaxis2-mzavg.qz0,[0.0,slope]), $
;        color=im_color(mstarcolor[mm]), psym=symcat(mstarpsym[mm]), $
;        symsize=mstarsymsize1[mm]*0.8; line=mstarline[mm], thick=8
    endfor
; now overlay the pegase model results
    for mm = 0, nmstar-1 do djs_oplot, zaxis, modeloh[mm,*], thick=10, $
      line=mstarline[mm], color=im_color(modelcolor[mm])

; ####################
; sfrm evolution
    djs_plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, position=pos[*,1], $
      xrange=xrange1, yrange=yrange2, xtitle=xtitle1, ytitle=ytitle2, $
      ytickinterval=0.5
    im_legend, mstarlabel, /right, /bottom, charsize=1.6, line=mstarline, $
      pspacing=1.8, thick=8, color=mstarcolor, psym=-mstarpsym, box=0, margin=0
    
; salim+07    
    psym_salim = 5
    symsize_salim = 3.0
    for mm = 0, nmstar-1 do djs_oplot, [0.1], [poly(mstar[mm],[-6.33,-0.35])+9], $
      color=im_color(mstarcolor[mm]), $
      psym=symcat(mstarpsym[mm],thick=6), symsize=mstarsymsize2[mm]
;     psym=symcat(psym_salim,thick=6), 

; elbaz+07
    psym_elbaz = 4
    symsize_elbaz = 3.0
    for mm = 0, nmstar-1 do djs_oplot, [0.07], $
      [poly(mstar[mm],[alog10(1.29)-10.0*0.77,0.77-1.0])+9], $
      color=im_color(mstarcolor[mm]), $
      psym=symcat(mstarpsym[mm],thick=6), symsize=mstarsymsize2[mm]
;     psym=symcat(psym_elbaz,thick=6), 

; Noeske+07, from Dunne+10
    psym_noeske = 9
    symsize_noeske = 2.0
    zz = [0.36,0.58]
;   zz = [0.36,0.58,0.78]
    for mm = 0, nmstar-1 do begin
       sfrm = alog10([2.3,4.0]*(10D^mstar[mm]/1D10)^0.67)-mstar[mm]
;      sfrm = alog10([2.3,4.0,5.5]*(10D^mstar[mm]/1D10)^0.67)-mstar[mm]
       djs_oplot, zz, sfrm+9, color=im_color(mstarcolor[mm]), $
         psym=symcat(mstarpsym[mm],thick=6), symsize=mstarsymsize2[mm]
;        psym=symcat(psym_noeske,thick=6), 
    endfor
    
; Oliver+10
    psym_oliver = 9
    symsize_oliver = 2.0
    zz = [0.1,0.25,0.35,0.45,0.55,0.7]
    norm = [-1.08,-0.93,-0.96,-0.84,-0.75,-0.61]
    bb = [0.0,0.01,-0.16,-0.12,-0.20,-0.26]
    for mm = 0, nmstar-1 do begin
       for iz = 0, n_elements(zz)-1 do begin
          sfrm = poly(mstar[mm]-11,[norm[iz],bb[iz]])
          djs_oplot, [zz[iz]], [sfrm], color=im_color(mstarcolor[mm]), $
            psym=symcat(mstarpsym[mm],thick=6), symsize=mstarsymsize2[mm]
;           psym=symcat(psym_oliver,thick=6),
       endfor
    endfor
    
; Karim+11
    psym_karim = 6
    symsize_karim = 1.8

    zz = [0.3,0.5,0.7]
    mm = 0
    mass = 9.99
    sfr = [1.5,3.4,5.4]
    sfrerr = [0.2,0.5,1.0]
    djs_oplot, zz, alog10(sfr)-mass+9, $ ; sfrerr/sfr/alog(10), $ ; 10^10
      color=im_color(mstarcolor[mm]), $;, errcolor=im_color(mstarcolor[mm]), errthick=!p.thick
      psym=symcat(mstarpsym[mm],thick=6), symsize=mstarsymsize2[mm]
;     psym=-symcat(psym_karim,thick=6), 

    mm = 1
    mass = 10.37
    sfr = [2.4,6.2,9.4]
    sfrerr = [0.3,1.0,0.8]
    djs_oplot, zz, alog10(sfr)-mass+9, $ ; sfrerr/sfr/alog(10), $ ; 10^10.5
      color=im_color(mstarcolor[mm]), $
;     psym=-symcat(psym_karim), symsize=symsize_karim
      psym=symcat(mstarpsym[mm],thick=6), symsize=mstarsymsize2[mm]
    
    mm = 2
    mass = 11.1
    sfr = [5.9,14.1,19.5]
    sfrerr = [1.7,2.0,3.5]
    djs_oplot, zz, alog10(sfr)-mass+9, $ ; sfrerr/sfr/alog(10), $ ; 10^11
      color=im_color(mstarcolor[mm]), $
;     psym=-symcat(psym_karim), symsize=symsize_karim
      psym=symcat(mstarpsym[mm],thick=6), symsize=mstarsymsize2[mm]
    
; now overlay the pegase model results
    for mm = 0, nmstar-1 do djs_oplot, zaxis, modelsfrm[mm,*], thick=8, $
      line=mstarline[mm], color=im_color(modelcolor[mm])

    im_plotconfig, psfile=psfile, /psclose

; ---------------------------------------------------------------------------
; Figure 22 - dependence of the MZ and SFR/M relations on our model
; parameters 
    psfile = pspath+'mzpegase_params'+suffix
    im_plotconfig, 2, pos, psfile=psfile, xspace=0.2, charsize=1.5, $
      height=[3.0,3.0], width=4.5*[1,1], xmargin=[1.1,1.1]

    zobs = 0.1
    zform0 = 3.0
    tau0 = 8.0

    ohrange = [8.05,9.35]
    massrange = [8.8,11.3]
    sfrmrange = [-1.8,-0.3]
    maxis1 = range(massrange[0]+0.05,massrange[1]-0.1,50)
    maxis2 = range(massrange[0]+0.05,11.0,50)

; ###############
; dependence of the MZ relation on alpha
    alpha = [0.5,1,2]
    line = lindgen(n_elements(alpha))
    beta = 0.0

; oh vs mass    
    plot, [0], [0], /nodata, xrange=massrange, yrange=ohrange, $
      xsty=1, ysty=1, xtitle='', ytitle=mzplot_ohtitle(), $
      position=pos[*,0], xtickinterval=1, xtickname=replicate(' ',10)
;   im_legend, '\beta = '+string(beta,format='(F4.1)'), $
;     /left, /top, box=0, charsize=1.4, margin=0
    im_legend, '\beta = '+string(beta,format='(F3.1)'), $
      /right, /bottom, box=0, charsize=1.4, position=[10.88,8.45], /data
    im_legend, ['\alpha = ','','','']+string(alpha,format='(F3.1)'), /bottom, $
      /right, box=0, line=line, pspacing=1.8, thick=6, charsize=1.4, margin=0
    im_legend, strupcase(calib), /left, /top, charsize=1.1, line=calibline, $
      pspacing=1.8, thick=6, color=calibcolor;, position=[10.0,8.8], /data

    for jj = 0, 2 do begin
       ohmodel = mzevol_func(maxis1,mzavg.mzevol_coeffs_r0zero[*,jj],z=zobs,qz0=mzavg.qz0)
       ww = where(ohmodel gt ohrange[0]+0.02)
       djs_oplot, maxis1[ww], ohmodel[ww], color=im_color(calibcolor[jj]), thick=10, $
         line=calibline[jj];, $
;        psym=symcat(calibpsym[jj]), symsize=calibsymsize[jj]*0.5
    endfor

    for ii = 0, n_elements(alpha)-1 do begin
       model = mzpegase_galform_model(zobs,alpha=alpha[ii],beta=beta,$
         zform0=zform0,tau0=tau0)
       ww = where(model.log12oh gt ohrange[0]+0.02 and model.mstar gt massrange[0]+0.05 and $
         model.mstar lt massrange[1]-0.05)
       djs_oplot, model[ww].mstar, model[ww].log12oh, line=line[ii], thick=6
    endfor       

; sfrm vs mass    
    plot, [0], [0], /nodata, /noerase, xrange=massrange, yrange=sfrmrange, $
      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle=textoidl('log (sSFR / Gyr^{-1})'), $
      position=pos[*,2], xtickinterval=1, ytickinterval=0.5
    djs_oplot, maxis2, poly(maxis2,[alog10(1.29)-10.0*0.77,0.77-1.0])+9, line=0, $
      thick=6, color=im_color('dark violet')
    djs_oplot, maxis2, poly(maxis2,[-6.33,-0.35])+9, line=0, thick=6, color=im_color('orange')
;   im_legend, ['Salim+07','Elbaz+07'], /left, /bottom, line=0, thick=6, $
;     color=['orange','dark violet'], pspacing=1.8, charsize=1.4, position=[9.1,-1.8]

    for ii = 0, n_elements(alpha)-1 do begin
       model = mzpegase_galform_model(zobs,alpha=alpha[ii],beta=beta,$
         zform0=zform0,tau0=tau0)
       djs_oplot, model.mstar, model.sfrm, line=line[ii], thick=6
    endfor       
    
; ###############
; dependence of the MZ relation on beta
    beta = [0.0,0.5,1.0]
    line = lindgen(n_elements(beta))
    alpha = 1.0

; oh vs mass
    plot, [0], [0], /nodata, /noerase, xrange=massrange, yrange=ohrange, $
      xsty=1, ysty=1, xtitle='', ytitle='', $
      position=pos[*,1], ytickname=replicate(' ',10), xtickinterval=1, xtickname=replicate(' ',10)
    axis, /yaxis, yrange=ohrange, ysty=1, ytitle=mzplot_ohtitle()
    im_legend, '\alpha = '+string(alpha,format='(F3.1)'), /left, $
      /top, box=0, charsize=1.4, margin=0, position=[10.22,8.55], /data
    im_legend, ['\beta = ','','','']+string(beta,format='(F3.1)'), /right, $
      /bottom, box=0, line=line, pspacing=1.8, thick=6, charsize=1.4, margin=0
    
    for jj = 0, 2 do begin
       ohmodel = mzevol_func(maxis1,mzavg.mzevol_coeffs_r0zero[*,jj],z=zobs,qz0=mzavg.qz0)
       ww = where(ohmodel gt ohrange[0]+0.02)
       djs_oplot, maxis1[ww], ohmodel[ww], color=im_color(calibcolor[jj]), thick=10, $
         line=calibline[jj];, $
;        psym=symcat(calibpsym[jj])
    endfor

    for ii = 0, n_elements(beta)-1 do begin
       model = mzpegase_galform_model(zobs,alpha=alpha,beta=beta[ii],$
         zform0=zform0,tau0=tau0)
       ww = where(model.log12oh gt ohrange[0]+0.02 and model.mstar gt massrange[0]+0.05 and $
         model.mstar lt massrange[1]-0.05)
       djs_oplot, model[ww].mstar, model[ww].log12oh, line=line[ii], thick=6
    endfor       

; sfrm vs mass
    plot, [0], [0], /nodata, /noerase, xrange=massrange, yrange=sfrmrange, $
      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle='', $
      position=pos[*,3], xtickinterval=1, ytickname=replicate(' ',10), $
      ytickinterval=0.5
    axis, /yaxis, yrange=sfrmrange, ysty=1, ytitle=textoidl('log (sSFR / Gyr^{-1})'), $
      ytickinterval=0.5
    djs_oplot, maxis2, poly(maxis2,[alog10(1.29)-10.0*0.77,0.77-1.0])+9, line=0, $
      thick=6, color=im_color('dark violet')
    djs_oplot, maxis2, poly(maxis2,[-6.33,-0.35])+9, line=0, thick=6, color=im_color('orange')
    im_legend, ['Salim+07','Elbaz+07'], /left, /bottom, line=0, thick=6, $
      color=['orange','dark violet'], pspacing=1.8, charsize=1.4, position=[9.2,-1.6]

    for ii = 0, n_elements(beta)-1 do begin
       model = mzpegase_galform_model(zobs,alpha=alpha,beta=beta[ii],$
         zform0=zform0,tau0=tau0)
       djs_oplot, model.mstar, model.sfrm, line=line[ii], thick=6
    endfor       
    
    im_plotconfig, psfile=psfile, /psclose

; ---------------------------------------------------------------------------
; *normalized* model and observed MZ relation as a function of
; redshift
    psfile = pspath+'mzpegase_mzevol'+suffix
    im_plotconfig, 0, pos, psfile=psfile, height=4.7, width=6.5, $
      xmargin=[1.6,0.4]

    xrange = [8.9,11.2]
    yrange = [-0.5,+0.1]
    xtitle1 = mzplot_masstitle()
;   ytitle1 = '12 + log (O/H) relative to 10^{10.5} M'+sunsymbol()
    ytitle1 = '\Delta'+'[12 + log (O/H)]'
    maxis1 = range(xrange[0]+0.02,xrange[1]-0.1,50)

;   zval = zbins.zbin
    zval = [0.1,0.2,0.3,0.4,0.5,0.7]
    zlabel = 'z='+string(zval,format='(F4.2)')
    zline = [0,3,5,1,4,2]
    zcolor = ['black','firebrick','dodger blue','tan','forest green','orange']

;   zthese = lindgen(nz)
    zthese = [0,2,4,5]
    nzthese = n_elements(zthese)

    alpha = 1.5
    beta = +0.0
    
    djs_plot, [0], [0], /nodata, ysty=1, xsty=1, position=pos, $
      xrange=xrange, yrange=yrange, xtitle=xtitle1, ytitle=ytitle1
    im_legend, zlabel[zthese], /right, /bottom, box=0, line=zline[zthese], $
      color=zcolor[zthese], pspacing=2, thick=6, margin=0, $
      charsize=1.4

    for iz = 0, nzthese-1 do begin
       mzdata = interpol(mzobs+mzavg.p0avg*(zval[iz]-mzavg.qz0),mzobsmass,maxis1)
       ww = where(mzdata gt yrange[0]+0.01)
       djs_oplot, maxis1[ww], mzdata[ww], line=zline[zthese[iz]], $
         color=im_color(zcolor[zthese[iz]]), thick=2

       model = mzpegase_galform_model(zval[iz],alpha=alpha,beta=beta,$
         mass0=mass0,zform0=zform0,tau0=tau0)
       mzmodel = interpol(model.log12oh,model.mstar,maxis1)
       if (iz eq 0) then mzmodelnorm = interpol(mzmodel,maxis1,normmass)
       mzmodel = mzmodel-mzmodelnorm
       
       ww = where(mzmodel gt yrange[0]+0.01)
       djs_oplot, maxis1[ww], mzmodel[ww], line=zline[zthese[iz]], $
         color=im_color(zcolor[zthese[iz]]), thick=10
    endfor
    
    im_plotconfig, psfile=psfile, /psclose

return
end



;; ---------------------------------------------------------------------------
;; Figure 22 - dependence of the MZ relation on our model parameters
;    psfile = pspath+'mzpegase_params'+suffix
;    im_plotconfig, 3, pos, psfile=psfile, xspace=0.0, charsize=1.5
;
;    zobs = 0.1
;    zform0 = 2.0
;    tau0 = 3.0
;
;    ohrange = [8.4,9.35]
;    massrange = [8.8,11.3]
;    maxis = range(massrange[0]+0.02,massrange[1]-0.1,50)
;
;; ###############
;; dependence of the MZ relation on alpha
;    alpha = [-0.1,-0.5,-1,-2]
;    line = lindgen(n_elements(alpha))
;    beta = 0.0
;    m0 = 10.5
;    mass0 = 10.0^m0
;
;    plot, [0], [0], /nodata, xrange=massrange, yrange=ohrange, $
;      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle=mzplot_ohtitle(), $
;      position=pos[*,0], xtickinterval=1
;    legend, textoidl('\alpha='+string(alpha,format='(F4.1)')), /right, $
;      /bottom, box=0, line=line, pspacing=1.8, thick=8, charsize=1.1, margin=0
;    im_legend, ['\beta='+string(beta,format='(F4.1)'),$
;      'log(M_{0}/M'+sunsymbol()+')='+string(m0,format='(F4.1)')], $
;      /left, /top, box=0, charsize=1.1, margin=0
;
;    for jj = 0, 2 do begin
;       ohmodel = mzevol_func(maxis,mzavg.mzevol_coeffs_r0zero[*,jj],z=zobs,qz0=mzavg.qz0)
;       ww = where(ohmodel gt ohrange[0]+0.02)
;       djs_oplot, maxis[ww], ohmodel[ww], color=im_color(calibcolor[jj]), thick=10, $
;         line=calibline[jj];, $
;;        psym=symcat(calibpsym[jj])
;    endfor
;
;    for ii = 0, n_elements(alpha)-1 do begin
;       model = mzpegase_galform_model(zobs,alpha=alpha[ii],beta=beta,$
;         mass0=mass0,zform0=zform0,tau0=tau0)
;       ww = where(model.log12oh gt ohrange[0]+0.02 and model.mstar gt massrange[0]+0.05 and $
;         model.mstar lt massrange[1]-0.05)
;       djs_oplot, model[ww].mstar, model[ww].log12oh, line=line[ii], thick=6
;    endfor       
;
;; ###############
;; dependence of the MZ relation on beta
;    beta = [0.0,0.2,0.5,0.8]
;    line = lindgen(n_elements(beta))
;    alpha = -1.0
;    m0 = 10.5
;    mass0 = 10.0^m0
;
;    plot, [0], [0], /nodata, /noerase, xrange=massrange, yrange=ohrange, $
;      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle='', $
;      position=pos[*,1], ytickname=replicate(' ',10), xtickinterval=1
;    legend, textoidl('\beta='+string(beta,format='(F3.1)')), /right, $
;      /bottom, box=0, line=line, pspacing=1.8, thick=8, charsize=1.1, margin=0
;    im_legend, ['\alpha='+string(alpha,format='(F4.1)'),$
;      'log(M_{0}/M'+sunsymbol()+')='+string(m0,format='(F4.1)')], /left, $
;      /top, box=0, charsize=1.1, margin=0
;
;    for jj = 0, 2 do begin
;       ohmodel = mzevol_func(maxis,mzavg.mzevol_coeffs_r0zero[*,jj],z=zobs,qz0=mzavg.qz0)
;       ww = where(ohmodel gt ohrange[0]+0.02)
;       djs_oplot, maxis[ww], ohmodel[ww], color=im_color(calibcolor[jj]), thick=10, $
;         line=calibline[jj];, $
;;        psym=symcat(calibpsym[jj])
;    endfor
;
;    for ii = 0, n_elements(beta)-1 do begin
;       model = mzpegase_galform_model(zobs,alpha=alpha,beta=beta[ii],$
;         mass0=mass0,zform0=zform0,tau0=tau0)
;       ww = where(model.log12oh gt ohrange[0]+0.02 and model.mstar gt massrange[0]+0.05 and $
;         model.mstar lt massrange[1]-0.05)
;       djs_oplot, model[ww].mstar, model[ww].log12oh, line=line[ii], thick=6
;    endfor       
;
;; ###############
;; dependence of the MZ relation on M0
;    m0 = [9.5,10.0,10.5,11.0]
;    mass0 = 10.0^m0
;    line = lindgen(n_elements(mass0))
;    alpha = -1.0
;    beta = +0.3
;
;    plot, [0], [0], /nodata, /noerase, xrange=massrange, yrange=ohrange, $
;      xsty=1, ysty=1, xtitle=mzplot_masstitle(), ytitle='', $
;      position=pos[*,2], ytickname=replicate(' ',10), xtickinterval=1
;    legend, textoidl('log(M_{0}/M'+sunsymbol()+')='+string(m0,format='(F4.1)')), /right, $
;      /bottom, box=0, line=line, pspacing=1.8, thick=8, charsize=1.1, margin=0
;    im_legend, ['\alpha='+string(alpha,format='(F4.1)'),$
;      '\beta='+string(beta,format='(F4.1)')], /left, /top, box=0, charsize=1.1, margin=0
;
;    for jj = 0, 2 do begin
;       ohmodel = mzevol_func(maxis,mzavg.mzevol_coeffs_r0zero[*,jj],z=zobs,qz0=mzavg.qz0)
;       ww = where(ohmodel gt ohrange[0]+0.02)
;       djs_oplot, maxis[ww], ohmodel[ww], color=im_color(calibcolor[jj]), thick=10, $
;         line=calibline[jj];, $
;;        psym=symcat(calibpsym[jj])
;    endfor
;
;    for ii = 0, n_elements(mass0)-1 do begin
;       model = mzpegase_galform_model(zobs,alpha=alpha,beta=beta,$
;         mass0=mass0[ii],zform0=zform0,tau0=tau0)
;       ww = where(model.log12oh gt ohrange[0]+0.02 and model.mstar gt massrange[0]+0.05 and $
;         model.mstar lt massrange[1]-0.05)
;       djs_oplot, model[ww].mstar, model[ww].log12oh, line=line[ii], thick=6
;    endfor       
;
;    im_plotconfig, psfile=psfile, /psclose
;
