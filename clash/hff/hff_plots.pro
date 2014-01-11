pro hff_plots

    isedfit_dir = getenv('CLASH_PROJECTS')+'/hff/'
    montegrids_dir = isedfit_dir+'montegrids/'

    suffix = 'jan03'
    
; --------------------------------------------------
; generate a P(z) QAplot comparing my photoz's with BPZ    
    prefix = 'hff_photoz'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    rr = read_isedfit(isedfit_paramfile,params=pp,isedfit_post=post)
    ngal = n_elements(rr)

; read the catalog    
    cat = rsex(isedfit_dir+'flx_iso.'+suffix)
;   cat.bpz = abs(cat.bpz)
;   cat = cat[where(cat.bpz ge 7.0 and cat.bpz lt 10.01)]
;   cat = cat[reverse(sort(cat.bpz))]

;   cat = cat[where(abs(cat.bpz) ge 7.0 and abs(cat.bpz) lt 10.01)]
;   cat = cat[reverse(sort(abs(cat.bpz)))]

    galaxy = 'ID '+string(cat.id,format='(I0)')
    bpz = rsex(isedfit_dir+'bpz/jan03.bpz')
    pbpz = read_bpz_probs(isedfit_dir+'bpz/jan03.probs',$
      redshift=bpz_redshift,dz=dz)
    niceprint, cat.id, bpz.id, bpz.z_b, rr.z
;   match2, cat.id, bpz.id, m1, m2
;   srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
;   bpz = bpz[m2]

; make the plot    
    psfile = isedfit_dir+'a2744_photoz_'+suffix+'.ps'

    xrange = [0.0,12.0]
;   yrange = [0,max(post.pofz)*1.05]
    
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, xmargin=[1.3,0.4], width=6.8
    for ii = 0, ngal-1 do begin
       pofz_ised = post[ii].pofz
       pofz_bpz = pbpz[ii].pofz*pp.zbin/dz

;      yrange = [0.0,max(post[ii].pofz)]
       yrange = [0.0,(max(pofz_ised)>max(pofz_bpz))*1.05]
       title = galaxy[ii]+', z_{iSEDfit}='+strtrim(string(rr[ii].z,format='(F12.3)'),2)+' '+$
         'z_{BPZ}='+strtrim(string(bpz[ii].z_b,format='(F12.3)'),2)
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='Redshift', ytitle='Posterior Probability', $
         xrange=xrange, yrange=yrange, title=title
       djs_oplot, pp.redshift, pofz_ised, line=0, thick=8, psym=10, color=cgcolor('firebrick')
       djs_oplot, bpz_redshift, pofz_bpz, line=0, color=cgcolor('dodger blue'), thick=8, psym=10
       djs_oplot, rr[ii].z*[1,1], !y.crange, thick=6
       djs_oplot, bpz[ii].z_b*[1,1], !y.crange, thick=6, line=5
       djs_oplot, cat[ii].bpz*[1,1], !y.crange, thick=6, line=5, color='orange'
;      im_legend, [galaxy[ii],'z_{iSEDfit}='+strtrim(string(rr[ii].z,format='(F12.3)'),2),$
;        'z_{BPZ}='+strtrim(string(zz[ii],format='(F12.3)'),2)], /left, /top, $
;        box=0, margin=0, textcolor=cgcolor(['black','firebrick','dodger blue'])
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
; --------------------------------------------------
; make a multi-panel plot for the paper adopting the BPZ redshifts
    prefix = 'hff'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    filt = hff_filterlist(pivotwave=weff,width=hwhm,/useirac)
    weff = weff/1D4
    hwhm = hwhm/1D4
    
    rr = read_isedfit(isedfit_paramfile)
    cat = rsex(isedfit_dir+'flx_iso.jan03')
    cat.bpz = abs(cat.bpz)

    cat = cat[where(cat.bpz ge 7.0 and cat.bpz lt 10.01)]
    cat = cat[reverse(sort(cat.bpz))]

    if total(abs(rr.z-cat.bpz) gt 1E-5) ne 0.0 then message, 'Mismatch!'
;   niceprint, cat.id, cat.bpz, rr.z

    these = where(cat.id eq 383 or cat.id eq 400 or $
      cat.id eq 407 or cat.id eq 429,nobj)
;   these = where(cat.id eq 400 or cat.id eq 270 or $
;     cat.id eq 292 or cat.id eq 2466,nobj)
;   these = where(cat.id eq 407 or cat.id eq 339 or $
;     cat.id eq 120 or cat.id eq 2338,nobj)
;   these = reverse(these)
    cat = cat[these]
    rr = read_isedfit(isedfit_paramfile,index=these,/getmodels)
    
    col = 'dodger blue'
    psfile = isedfit_dir+'a2744_seds_'+suffix+'.eps'
    xrange = [0.3,6.0]
    im_plotconfig, 5, pos, psfile=psfile, xmargin=[1.1,0.2], $
      height=[2.6,2.6], xspace=0.05, yspace=0.05
    for ii = 0, nobj-1 do begin
;      if strmatch(filt,'*acs*') then col = 'orange'
;      if strmatch(filt,'*wfc3*') then col = 'dodger blue'
;      if strmatch(filt,'*irac*') then col = 'brown'
       
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
         strtrim(string(cat[ii].bpz,format='(F12.2)'),2)]
    endfor 
    
    xyouts, pos[0,0]-0.08, pos[1,0], 'AB Magnitude', align=0.5, orientation=90, /normal
    xyouts, pos[0,3], pos[1,3]-0.12, textoidl('Observed-Frame Wavelength (\mu'+'m)'), $
      align=0.5, /normal

    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

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
    
