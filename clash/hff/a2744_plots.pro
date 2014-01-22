pro a2744_plots

    isedfit_dir = getenv('CLASH_PROJECTS')+'/hff/a2744/'
    montegrids_dir = isedfit_dir+'montegrids/'

    suffix = 'jan03'

    filt = hff_filterlist(pivotwave=weff,width=hwhm,/useirac)
    weff = weff/1D4
    hwhm = hwhm/1D4
    
; --------------------------------------------------
; QAplot: compare the BPZ & iSEDfit best redshifts
    isedfit_paramfile = isedfit_dir+'a2744_photoz_paramfile.par'

    cat = read_a2744(/photoz,bpz_dz=bpz_dz,bpz_redshift=bpz_redshift)
    rr = read_isedfit(isedfit_paramfile,params=pp,isedfit_dir=isedfit_dir)
;   keep = where(cat.group lt 3)
;   cat = cat[keep]
;   rr = rr[keep]
    ngal = n_elements(rr)

; print out some statistics    
    niceprint, cat.id, cat.group, cat.zb1, rr.z, rr.z-cat.zb1
    diff = rr.z-cat.zb1
    good = where(abs(diff) lt 1.0,comp=bad)
    splog, djs_mean(diff[good]), djsig(diff[good])
    niceprint, cat[bad].id, cat[bad].group, cat[bad].zb1, rr[bad].z, rr[bad].mstar_50, $
      rr[bad].tau_50, rr[bad].av_50, rr[bad].age_50, rr[bad].sfr_50
    niceprint, cat[bad].id, cat[bad].group, cat[bad].zb1, rr[bad].z

stop    
    
    psfile = isedfit_dir+'bpz_vs_isedfit.ps'

    xrange = [0.0,11.0]
    yrange = xrange
    
    im_plotconfig, 0, pos, psfile=psfile, height=5.0
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='BPZ Redshift', ytitle='iSEDfit Redshift', $
      xrange=xrange, yrange=yrange
    djs_oplot, !x.crange, !y.crange, line=0, color=cgcolor('grey')
    oploterror, cat.zb1, rr.z, rr.z_err, psym=symcat(15)

    ww = where(cat.zb1 gt 6 and rr.z lt 4.0,nww)
    djs_oplot, cat[ww].zb1, max(rr[ww].photoz,dim=1), psym=symcat(6)
    for ii = 0, nww-1 do djs_oplot, [cat[ww[ii]].zb1,cat[ww[ii]].zb1], $
      [rr[ww[ii]].photoz[0],max(rr[ww[ii]].photoz)], color=cgcolor('grey'), line=0
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

; --------------------------------------------------
; paper plot: multi-panel SED plot adopting the BPZ redshifts
    isedfit_paramfile = isedfit_dir+'a2744_paramfile.par'

; print out some properties for the full sample of 10 objects for the
; paper... 
    cat = read_a2744()
    these = where(cat.group eq 1 and (cat.irac_ch1_flux gt 0.0 or $
      cat.irac_ch2_flux gt 0.0),nobj)
    cat = cat[these]
    rr = read_isedfit(isedfit_paramfile,index=these,isedfit_post=post,$
      isedfit_dir=isedfit_dir)
    
    srt = reverse(sort(cat.zb1))
    cat = cat[srt]
    rr = rr[srt]
    post = post[srt]

    logmu = alog10(cat.mu)

    rr.mstar_50 = rr.mstar_50-logmu
    rr.sfr_50 = 10^rr.sfr_50-logmu
    rr.sfr_err = rr.sfr_err*rr.sfr_50*alog(10)
    sfrage95 = fltarr(nobj)
    for ii = 0, nobj-1 do sfrage95[ii] = weighted_quantile(post[ii].sfrage,quant=95)

    zform = getredshift(getage(cat.zb1)-sfrage95)
    
    splog, 'ID, mass, SFR, sSFR, age, age95, t(z), zform'
    niceprint, cat.id, rr.mstar_50, rr.sfr_50, 1D9*rr.sfr_50/10.0^rr.mstar, $
      rr.sfrage_50*1E3, sfrage95*1E3, getage(cat.zb1)*1E3, zform
    print

    splog, 'Median, minmax mu', djs_median(logmu), minmax(logmu)
    splog, 'Median, minmax mass, median err ', djs_median(10^rr.mstar_50), $
      minmax(10^rr.mstar_50), djs_median(rr.mstar_err*10^rr.mstar*alog(10))
    splog, 'Median, minmax SFR, median err ', djs_median(rr.sfr_50), $
      minmax(rr.sfr_50), djs_median(rr.sfr_err)
    splog, 'Median, minmax zform', djs_median(zform), minmax(zform)
    splog, 'Median, minmax sfrage', djs_median(rr.sfrage), minmax(rr.sfrage)
    splog, 'Median, minmax sfrage95', djs_median(sfrage95), minmax(sfrage95)
    splog, 'Doubling time ', 2*djs_median(10^rr.mstar_50)/djs_median(rr.sfr_50)/1D6
    splog, 'Median redshift, median age of Universe', djs_median(cat.zb1), $
      getage(djs_median(cat.zb1))
    splog, 'Median formation redshift ', getredshift(getage(djs_median(cat.zb1))-$
      djs_median(sfrage95) )
    djs_plot, 1E3*getage(cat.zb1), 2*10^rr.mstar_50/rr.sfr_50/1D6, $
      psym=8, xr=[0,1000], yr=[0,1000]
    djs_oplot, !x.crange, !y.crange

stop    
    
; ...but only make a plot for four representative objects
    cat = read_a2744()
    these = where(cat.id eq 406 or cat.id eq 270 or $
      cat.id eq 292 or cat.id eq 3903,nobj)
    cat = cat[these]
    rr = read_isedfit(isedfit_paramfile,index=these,/getmodels)
    
    srt = reverse(sort(cat.zb1))
    cat = cat[srt]
    rr = rr[srt]

    col = 'dodger blue'
    psfile = isedfit_dir+'a2744_seds.ps'
    xrange = [0.3,6.0]
    im_plotconfig, 5, pos, psfile=psfile, xmargin=[1.1,0.2], $
      height=[2.6,2.6], xspace=0.05, yspace=0.05
    for ii = 0, nobj-1 do begin
       if odd(ii) then ytickname = replicate(' ',10) else delvarx, ytickname
       if ii le 1 then xtickname = replicate(' ',10) else delvarx, xtickname
       ticks = loglevels(xrange)
       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=[31,23], $
         xsty=1, ysty=1, /xlog, position=pos[*,ii], noerase=ii gt 0, $
         ytickname=ytickname, xtickname=xtickname, xtickv=ticks, $
         xticks=n_elements(ticks)-1, ytickinterval=2
;      djs_oplot, rr[ii].wave/1D4, rr[ii].flux
       mask = (rr[ii].wave gt 1214.5*(1+rr[ii].z) and $
         rr[ii].wave lt 1215*1.0025*(1+rr[ii].z)) eq 1
       djs_oplot, rr[ii].wave/1D4, djs_maskinterp(rr[ii].flux,mask,rr[ii].wave/1D4)
; model photometry
       notzero = where(rr[ii].bestmaggies gt 0.0)
       bestmab = -2.5*alog10(rr[ii].bestmaggies[notzero])
       djs_oplot, weff[notzero], bestmab, psym=symcat(6,thick=6), symsize=2.5, $
         color=cgcolor('grey')
; overplot the data; distinguish between three different cases, based
; on the input photometry
       mab = maggies2mag(rr[ii].maggies,nsigma=2.0,$
         ivar=rr[ii].ivarmaggies,magerr=maberr,$
         lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
       used = where(mab gt -90.0,nused)
       upper = where(mab lt -90.0 and mabupper gt -90,nupper)
       if (nused ne 0L) then begin
          oploterror, weff[used], mab[used], hwhm[used], $
            mabhierr[used], psym=symcat(16), $
            symsize=1.7, color=im_color(col), /hibar, $
            errcolor=im_color(col), errthick=!p.thick
          oploterror, weff[used], mab[used], hwhm[used], $
            mabloerr[used], psym=3, color=im_color('dodger blue'), /lobar, $
            errcolor=im_color(col), errthick=!p.thick
       endif
       if (nupper ne 0) then begin
          djs_oplot, weff[upper], mabupper[upper], $
            psym=symcat(11,thick=6), symsize=3.0, color=im_color('forest green')
       endif
       im_legend, /left, /top, box=0, spacing=1.5, charsize=1.5, margin=0, $
         ['ID'+strtrim(cat[ii].id,2),'z = '+$
         strtrim(string(cat[ii].zb1,format='(F12.1)'),2)]
    endfor 
    
    xyouts, pos[0,0]-0.08, pos[1,0], 'AB Magnitude', align=0.5, orientation=90, /normal
    xyouts, pos[0,3], pos[1,3]-0.12, textoidl('Observed-Frame Wavelength (\mu'+'m)'), $
      align=0.5, /normal

    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
; --------------------------------------------------
; QAplot: multi-panel SED plot adopting the BPZ redshifts
    isedfit_paramfile = isedfit_dir+'a2744_paramfile.par'

    cat = read_a2744()
    these = where(cat.group eq 1 and (cat.irac_ch1_flux gt 0.0 or $
      cat.irac_ch2_flux gt 0.0),nobj)
;   these = where(cat.id eq 383 or cat.id eq cat.id eq cat.id eq 
    cat = cat[these]
    rr = read_isedfit(isedfit_paramfile,index=these,/getmodels)
    
    srt = reverse(sort(cat.zb1))
    cat = cat[srt]
    rr = rr[srt]

    col = 'dodger blue'
    psfile = isedfit_dir+'a2744_seds_'+suffix+'.ps'
    xrange = [0.3,6.0]
;   im_plotconfig, 5, pos, psfile=psfile, xmargin=[1.1,0.2], $
;     height=[2.6,2.6], xspace=0.05, yspace=0.05
    im_plotconfig, 23, pos, psfile=psfile
    for ii = 0, nobj-1 do begin
       if odd(ii) then ytickname = replicate(' ',10) else delvarx, ytickname
       if ii le 7 then xtickname = replicate(' ',10) else delvarx, xtickname
       ticks = loglevels(xrange)
       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=[31,23], $
         xsty=1, ysty=1, /xlog, position=pos[*,ii], noerase=ii gt 0, $
         ytickname=ytickname, xtickname=xtickname, xtickv=ticks, $
         xticks=n_elements(ticks)-1, ytickinterval=2
;      djs_oplot, rr[ii].wave/1D4, rr[ii].flux
       mask = (rr[ii].wave gt 1214.5*(1+rr[ii].z) and $
         rr[ii].wave lt 1215*1.0025*(1+rr[ii].z)) eq 1
       djs_oplot, rr[ii].wave/1D4, djs_maskinterp(rr[ii].flux,mask,rr[ii].wave/1D4)
; model photometry
       notzero = where(rr[ii].bestmaggies gt 0.0)
       bestmab = -2.5*alog10(rr[ii].bestmaggies[notzero])
       djs_oplot, weff[notzero], bestmab, psym=symcat(6,thick=6), symsize=2.5, $
         color=cgcolor('grey')
; overplot the data; distinguish between three different cases, based
; on the input photometry
       mab = maggies2mag(rr[ii].maggies,nsigma=2.0,$
         ivar=rr[ii].ivarmaggies,magerr=maberr,$
         lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
       used = where(mab gt -90.0,nused)
       upper = where(mab lt -90.0 and mabupper gt -90,nupper)
       if (nused ne 0L) then begin
          oploterror, weff[used], mab[used], hwhm[used], $
            mabhierr[used], psym=symcat(16), $
            symsize=1.7, color=im_color(col), /hibar, $
            errcolor=im_color(col), errthick=!p.thick
          oploterror, weff[used], mab[used], hwhm[used], $
            mabloerr[used], psym=3, color=im_color('dodger blue'), /lobar, $
            errcolor=im_color(col), errthick=!p.thick
       endif
       if (nupper ne 0) then begin
          djs_oplot, weff[upper], mabupper[upper], $
            psym=symcat(11,thick=6), symsize=3.0, color=im_color('forest green')
       endif
       im_legend, /left, /top, box=0, spacing=1.5, charsize=1.5, margin=0, $
         ['ID'+strtrim(cat[ii].id,2),'z = '+$
         strtrim(string(cat[ii].zb1,format='(F12.2)'),2)]
    endfor 
    
    xyouts, pos[0,0]-0.08, (pos[3,4]-pos[1,4])/2.0+pos[1,4], 'AB Magnitude', $
      align=0.5, orientation=90, /normal
    xyouts, pos[2,8], pos[1,8]-0.08, textoidl('Observed-Frame Wavelength (\mu'+'m)'), $
      align=0.5, /normal

    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
; --------------------------------------------------
; paper plot: P(z)'s
    isedfit_paramfile = isedfit_dir+'a2744_photoz_paramfile.par'

; just show groups 1&2
    cat = read_a2744(bpz_dz=bpz_dz,bpz_redshift=bpz_redshift,index=index)
    rr = read_isedfit(isedfit_paramfile,params=pp,isedfit_post=post,index=index)
    srt = reverse(sort(cat.zb1))
    cat = cat[srt]
    rr = rr[srt]
    
    ngal = n_elements(rr)

    yrange = [0,0.2]
    xrange = [0.5,11.5]
;   yrange = [0,max(post.pofz)*1.05]
    
    psfile = isedfit_dir+'pofz_'+suffix+'.ps'
    im_plotconfig, 20, pos, psfile=psfile, charsize=1.3, $
      height=1.6*[1,1,1], xspace=0.02*[1,1,1,1,1], yspace=0.02*[1,1]
    for ii = 0, ngal-1 do begin

       if ii gt 10 then delvarx, xtickname else xtickname = replicate(' ',10)
       
       pofz_ised = post[ii].pofz
       pofz_bpz = cat[ii].pofz*pp.zbin/bpz_dz
;      pofz_bpz = pofz_bpz/total(pofz_bpz)

       pofz_ised = post[ii].pofz/im_integral(pp.redshift,post[ii].pofz)
       pofz_bpz = cat[ii].pofz/im_integral(bpz_redshift,cat[ii].pofz)

;      yrange = [0.0,max(post[ii].pofz)]
       yrange = [0.0,(max(pofz_ised)>max(pofz_bpz))*1.05]
       djs_plot, [0], [0], /nodata, position=pos[*,ii], noerase=ii gt 0, $
         xsty=5, ysty=5, xrange=xrange, yrange=yrange
       polyfill, [pp.redshift,reverse(pp.redshift)], [pofz_ised,pofz_ised*0], $
         /fill, color=cgcolor('grey'), noclip=0
       polyfill, [bpz_redshift,reverse(bpz_redshift)], [pofz_bpz,pofz_bpz*0], $
         /fill, color=cgcolor('tomato'), noclip=0
       djs_oplot, pp.redshift, pofz_ised, thick=4
       djs_oplot, bpz_redshift, pofz_bpz, color=cgcolor('firebrick'), thick=4
       djs_plot, [0], [0], /nodata, position=pos[*,ii], /noerase, $
         xsty=1, ysty=1, xrange=xrange, yrange=yrange, ytickname=replicate(' ',10), $
         xtickname=xtickname, yticks=1

;      djs_oplot, rr[ii].z*[1,1], !y.crange, thick=6
;      djs_oplot, cat[ii].zb1*[1,1], !y.crange, thick=6, line=5
;      djs_oplot, cat[ii].bpz*[1,1], !y.crange, thick=6, line=5, color='orange'
       if cat[ii].group eq 2 then ss = '^{*}' else ss = ''
       
       im_legend, [strtrim(cat[ii].id,2)+ss,$
;        'z_{iSEDfit}='+strtrim(string(rr[ii].z,format='(F12.1)'),2),$
;        'z_{BPZ}='+strtrim(string(cat[ii].zb1,format='(F12.1)'),2)], /left, /top, $
         'z='+strtrim(string(cat[ii].zb1,format='(F12.1)'),2)], /left, /top, $
         box=0, margin=-0.2, charsize=1.1; textcolor=cgcolor(['black','grey']), $
         
    endfor

    xyouts, pos[0,0]-0.03, djs_mean(pos[[1,3],6]), $
      'Redshift Probability', align=0.5, orientation=90, $
      charsize=1.8, /norm
    xyouts, pos[2,14], pos[1,14]-0.1, 'Redshift', align=0.5, $
      charsize=1.8, /norm
;   xyouts, pos[2,13], pos[1,13]-0.1, 'Redshift', align=0.5, $
;     charsize=1.8, /norm
;   xyouts, pos[2,15], pos[1,15]-0.1, 'Redshift', align=0.5, $
;     charsize=1.8, /norm
    
    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

stop    
    
    
    
; --------------------------------------------------
; QAplot: comparison of BPZ & iSEDfit P(z)'s
    isedfit_paramfile = isedfit_dir+'a2744_photoz_paramfile.par'

    cat = read_a2744(/photoz,bpz_dz=bpz_dz,bpz_redshift=bpz_redshift)
    rr = read_isedfit(isedfit_paramfile,params=pp,isedfit_post=post)
    ngal = n_elements(rr)

; make the plot    
    psfile = isedfit_dir+'a2744_photoz_'+suffix+'.ps'

    xrange = [0.0,12.0]
;   yrange = [0,max(post.pofz)*1.05]
    
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, xmargin=[1.3,0.4], width=6.8
    for ii = 0, ngal-1 do begin
       pofz_ised = post[ii].pofz
       pofz_bpz = cat[ii].pofz*pp.zbin/bpz_dz

;      yrange = [0.0,max(post[ii].pofz)]
       yrange = [0.0,(max(pofz_ised)>max(pofz_bpz))*1.05]
       title = cat[ii].galaxy+', z_{iSEDfit}='+strtrim(string(rr[ii].z,format='(F12.3)'),2)+' '+$
         'z_{BPZ}='+strtrim(string(cat[ii].zb1,format='(F12.3)'),2)
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='Redshift', ytitle='Posterior Probability', $
         xrange=xrange, yrange=yrange, title=title
       djs_oplot, pp.redshift, pofz_ised, line=0, thick=8, psym=10, color=cgcolor('firebrick')
       djs_oplot, bpz_redshift, pofz_bpz, line=0, color=cgcolor('dodger blue'), thick=8, psym=10
       djs_oplot, rr[ii].z*[1,1], !y.crange, thick=6
       djs_oplot, cat[ii].zb1*[1,1], !y.crange, thick=6, line=5
       djs_oplot, cat[ii].bpz*[1,1], !y.crange, thick=6, line=5, color='orange'
;      im_legend, [galaxy[ii],'z_{iSEDfit}='+strtrim(string(rr[ii].z,format='(F12.3)'),2),$
;        'z_{BPZ}='+strtrim(string(zz[ii],format='(F12.3)'),2)], /left, /top, $
;        box=0, margin=0, textcolor=cgcolor(['black','firebrick','dodger blue'])
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
; --------------------------------------------------
; write out a table for the paper adopting the BPZ redshifts
    prefix = 'hff'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    rr = read_isedfit(isedfit_paramfile,params=pp,isedfit_post=post)
    ngal = n_elements(rr)

    cat = rsex(isedfit_dir+'flx_iso.jan03')
    cat.bpz = abs(cat.bpz)

    cat = cat[where(cat.bpz ge 7.0 and cat.bpz lt 10.01)]
    cat = cat[reverse(sort(cat.bpz))]
    if total(abs(rr.z-cat.bpz) gt 1E-5) ne 0.0 then message, 'Mismatch!'
;   niceprint, cat.id, cat.bpz, rr.z

    index = where(cat.id eq 383 or cat.id eq 400 or $
      cat.id eq 407 or cat.id eq 429,nobj)
;   index = where(cat.irac_ch1_flux gt 0.0 or cat.irac_ch2_flux gt 0.0,nobj)
    cat = cat[index]
    rr = rr[index]
    post = post[index]

    mag = rsex(isedfit_dir+'model.txt')
    match, cat.id, mag.id, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    mu = mag[m2].nfw ; magnification
    
;   niceprint, cat.id, cat.bpz, rr.z, mu
;   logmu = alog10(mu)
    logmu = alog10(djs_mean(mu))

    rr.mstar_50 = rr.mstar_50-logmu
    rr.sfr_50 = 10^rr.sfr_50-logmu
    rr.sfr_err = rr.sfr_err*rr.sfr_50*alog(10)
    sfrage95 = fltarr(nobj)
    for ii = 0, nobj-1 do sfrage95[ii] = weighted_quantile(post[ii].sfrage,quant=95)

    zform = getredshift(getage(cat.bpz)-sfrage95)
;   splog, 'ID, mass (demag), masserr, SFR (demag), SFRerr, age, ageerr, age95, t(z)'
;   niceprint, cat.id, rr.mstar_50, rr.mstar_err, rr.sfr_50, rr.sfr_err, $
;     rr.sfrage_50*1E3, rr.sfrage_err*1E3, sfrage95*1E3, getage(cat.bpz)*1E3, $
;     zform
;   print

;   use95 = lonarr(nobj)
;   use95[where(cat.id eq 383 or cat.id eq 429 or cat.id eq 292 or $
;     cat.id eq 527 or cat.id eq 120 or cat.id eq 2338)] = 1
;   openw, lun, isedfit_dir+'a2744_isedfit_dec17.txt', /get_lun
;   printf, lun, '# 1 ID'
;   printf, lun, '# 2 Mstar [log Msun]'
;   printf, lun, '# 3 SFR [Msun/yr]'
;   printf, lun, '# 4 Age [Myr]'
;   printf, lun, '# 5 Age95 [Myr]'
;   printf, lun, '# 6 Useage95 [1=yes]'
;   niceprintf, lun, cat.id, rr.mstar_50, rr.sfr_50, $
;     rr.sfrage_50*1E3, sfrage95*1E3, use95
;   free_lun, lun
    
    splog, 'ID, mass, SFR, sSFR, age, age95, t(z), zform'
    niceprint, cat.id, rr.mstar_50, rr.sfr_50, 1D9*rr.sfr_50/10.0^rr.mstar, $
      rr.sfrage_50*1E3, sfrage95*1E3, getage(cat.bpz)*1E3, zform
    print

    splog, 'Mean, minmax mu', djs_median(mu), minmax(mu)
    splog, 'Mean, minmax mass, mean err ', djs_median(10^rr.mstar_50), $
      minmax(10^rr.mstar_50), djs_median(rr.mstar_err*10^rr.mstar*alog(10))
    splog, 'Mean, minmax SFR, mean err ', djs_median(rr.sfr_50), $
      minmax(rr.sfr_50), djs_median(rr.sfr_err)
    splog, 'Mean, minmax zform', djs_median(zform), minmax(zform)
    splog, 'Mean, minmax sfrage95', djs_median(sfrage95), minmax(sfrage95)
    splog, 'Doubling time ', 2*djs_median(10^rr.mstar_50)/djs_median(rr.sfr_50)/1D6
    splog, djs_median(rr.z), getage(djs_median(rr.z))

    
stop

return
end

