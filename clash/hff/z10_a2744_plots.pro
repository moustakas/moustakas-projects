pro z10_a2744_plots

    isedfit_dir = getenv('CLASH_PROJECTS')+'/hff/z10_a2744/'
    montegrids_dir = isedfit_dir+'montegrids/'

    filt = hff_filterlist(pivotwave=weff,width=hwhm,/useirac)
    weff = weff/1D4
    hwhm = hwhm/1D4

; --------------------------------------------------
; paper plot: SED + P(z)
    isedfit_paramfile = isedfit_dir+'z10_a2744_paramfile.par'
    isedfit_paramfile_photoz = isedfit_dir+'z10_a2744_photoz_paramfile.par'

    cat = read_z10_a2744()
    these = [0,1,2,n_elements(cat)-1]
    cat = cat[these]
    nobj = n_elements(cat)

    common com_z10_a2744, rr, rr_photoz, post, pp, post1
    if n_elements(rr) eq 0L then begin
       rr = read_isedfit(isedfit_paramfile,index=these,/getmodels,$
         isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir)
       rr_photoz = read_isedfit(isedfit_paramfile_photoz,index=these,$
         isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir,$
         params=pp,isedfit_post=post)
       post1 = post
    endif
    col = 'orange'

; print out the physical properties
    
    
    
; pre-process the photoz's
    bpz = rsex(isedfit_dir+'bpz.out')
    pbpz = read_bpz_probs(isedfit_dir+'bpz.probs',redshift=bpz_zz)
    pbpz = pbpz[these]
    bpz = bpz[these]

    for ii = 0, nobj-1 do begin
       zmin = isedfit_find_zmin(bpz_zz,-2.0*alog(pbpz[ii].pofz>1D-20))
       struct_print, zmin & print & help, bpz[ii], /str
    endfor

; make the plot    
    xrange = [0.3,8.0]
    yrange = [32,23]

    psfile = isedfit_dir+'z10_a2744_isedfit.eps'
    im_plotconfig, 19, pos, psfile=psfile, xmargin=[1.1,0.3], $
      width=[4.5,2.5], height=2.1*[1,1,1,1], xspace=0.1, yspace=[0.05,0.05,0.05]
    for ii = 0, nobj-1 do begin
;      if odd(ii) then ytickname = replicate(' ',10) else delvarx, ytickname
       if ii le 2 then xtickname = replicate(' ',10) else delvarx, xtickname
       ticks = loglevels(xrange)
       if ii eq 3 then xtitle = textoidl('Observed-frame Wavelength (\mu'+'m)')
       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
         xsty=1, ysty=1, /xlog, position=pos[*,2*ii], noerase=ii gt 0, $
         ytickname=ytickname, xtickname=xtickname, xtickv=ticks, $
         xticks=n_elements(ticks)-1, ytickinterval=3, xtitle=xtitle
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
       print, rr[ii].maggies*sqrt(rr[ii].ivarmaggies)
       nsigma = rr[ii].maggies*0+2.0
;      if ii eq 1 then nsigma[1] = 2.49
;      if ii eq 3 then nsigma[2] = 2.1
       mab = maggies2mag(rr[ii].maggies,nsigma=nsigma,$
         ivar=rr[ii].ivarmaggies,magerr=maberr,$
         lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
       used = where(mab gt -90.0,nused)
       upper = where(mab lt -90.0 and mabupper gt -90,nupper)
       if (nused ne 0L) then begin
          oploterror, weff[used], mab[used], hwhm[used], $
            mabhierr[used], psym=symcat(16), $
            symsize=1.7, color=cgcolor(col), /hibar, $
            errcolor=cgcolor(col), errthick=3
          oploterror, weff[used], mab[used], hwhm[used], $
            mabloerr[used], psym=3, color=cgcolor(col), /lobar, $
            errcolor=cgcolor(col), errthick=3
       endif
       if (nupper ne 0) then begin
;         th = 6
;         for uu = 0, nupper-1 do begin
;            cgarrow, weff[upper[uu]], mabupper[upper[uu]], weff[upper[uu]], $
;              mabupper[upper[uu]]+2, /data, /solid, hsize=!d.x_size/50, $
;              thick=th, color=cgcolor('forest green')
;            djs_oplot, weff[upper[uu]]+*[-1,1], thick=th, $
;              mabupper[upper[uu]]*[1,1], color=cgcolor('forest green')
;         endfor
          djs_oplot, weff[upper], mabupper[upper], $
            psym=symcat(11,thick=6), symsize=2, color=im_color('forest green')
       endif
       im_legend, /left, /top, box=0, spacing=1.5, charsize=1.5, $ ; margin=0, $
         strtrim(cat[ii].galaxy,2)

; P(z) inset
       ipofz = post1[ii].pofz/im_integral(pp.redshift,post1[ii].pofz)
       bpofz = pbpz[ii].pofz/im_integral(bpz_zz,pbpz[ii].pofz)
       yr = [0,max(ipofz)>max(bpofz)]*1.02

;      pos2 = [pos[0,ii]+0.05,pos[1,ii]+0.06,pos[2,ii]-0.02,pos[1,ii]+0.18]
;      pos2 = [pos[2,ii]-0.18,pos[1,ii]+0.06,pos[2,ii]-0.02,pos[1,ii]+0.18]
       djs_plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, /norm, $
         position=pos[*,2*ii+1], xrange=[0,12], yrange=yr

       polyfill, [pp.redshift,reverse(pp.redshift)], [ipofz,ipofz*0+!y.crange[0]], $
         /fill, color=cgcolor('tomato'), noclip=0
       polyfill, [bpz_zz,reverse(bpz_zz)], [bpofz,bpofz*0+!y.crange[0]], $
         /fill, color=cgcolor('dodger blue'), noclip=0

       djs_oplot, pp.redshift, ipofz, line=0, color=cgcolor('firebrick')
       djs_oplot, bpz_zz, bpofz, line=0, color=cgcolor('navy')

       if ii eq 3 then xtitle = 'Redshift z'
       djs_plot, [0], [0], /nodata, /noerase, xsty=9, ysty=5, /norm, $
         position=pos[*,2*ii+1], xtitle=xtitle, ytitle='', $ ; charsize=0.8, $
         xrange=[0,12], yrange=yr, yticks=1, $
         ytickname=replicate(' ',10), xtickname=xtickname
;      xyouts, pos2[0,0]-0.01, djs_mean([pos2[1,0],pos2[3,0]]), 'Probability', $
;        align=0.5, orientation=90, charsize=0.8, charthick=3.0, /norm

       im_legend, ['z_{iSEDfit}='+$
         strtrim(string(rr_photoz[ii].z,format='(F12.2)'),2)+$
         '^{+'+strtrim(string(rr_photoz[ii].z_95[0],format='(F12.2)'),2)+'}'+$
         '_{-'+strtrim(string(rr_photoz[ii].z_95[1],format='(F12.2)'),2)+'}',$
         'z_{BPZ}='+strtrim(string(bpz[ii].z_b,format='(F12.2)'),2)+$
         '^{+'+strtrim(string(bpz[ii].z_b_max-bpz[ii].z_b,format='(F12.2)'),2)+'}'+$
         '_{-'+strtrim(string(bpz[ii].z_b-bpz[ii].z_b_min,format='(F12.2)'),2)+'}'], $
         /left, /top, textcolor=cgcolor(['tomato','dodger blue']), $
         color=cgcolor(['tomato','dodger blue']), charsize=1.3, $
         spacing=2.2, box=0, position=[pos[0,2*ii+1]-0.01,pos[3,2*ii+1]-0.02], /norm
;      im_legend, ['z_{iSEDfit}='+strtrim(string(rr_photoz[ii].z,format='(F12.2)'),2),$
;        'z_{BPZ}='+strtrim(string(bpz[ii].z_b,format='(F12.2)'),2)], $
;        /left, /top, textcolor=cgcolor(['tomato','dodger blue']), $
;        color=cgcolor(['tomato','dodger blue']), charsize=1.5, $
;        spacing=2.0, box=0, margin=0
       
    endfor 
    
    xyouts, pos[0,0]-0.08, pos[1,0], 'AB Magnitude', align=0.5, orientation=90, /normal
    xyouts, pos[0,4]-0.08, pos[1,4], 'AB Magnitude', align=0.5, orientation=90, /normal
;   xyouts, pos[0,3], pos[1,3]-0.11, textoidl('Observed-frame Wavelength (\mu'+'m)'), $
;     align=0.5, /normal

    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
; --------------------------------------------------
; paper plot: P(z)'s
    isedfit_paramfile = isedfit_dir+'z10_a2744_photoz_paramfile.par'

    these = [0,1,10]
    rr = read_isedfit(isedfit_paramfile,params=pp,$
      isedfit_post=post,index=these)
    cat = read_z10_a2744(photoz=photoz)
    cat = cat[these]
    ngal = n_elements(cat)

    color = ['orange','tomato','powder blue']
    lcolor = ['red','firebrick','black']
    lthick = [6,6,8]
;   color = ['orange','tomato','powder blue','orchid']
;   lcolor = ['red','firebrick','blue','black']
;   lthick = [6,6,6,8]
;   color = ['orange','tomato','dodger blue','orchid','tan']
;   lcolor = ['red','firebrick','navy','purple','brown']
    
    for ii = 0, ngal-1 do post[ii].pofz = alog(post[ii].pofz/total(post[ii].pofz))
;   for ii = 0, ngal-1 do post[ii].pofz = post[ii].pofz/total(post[ii].pofz)
;   for ii = 0, ngal-1 do post[ii].pofz = post[ii].pofz/$
;     im_integral(pp.redshift,post[ii].pofz)
    pofz_final = post[0].pofz*0+1
    for ii = 0, 1 do pofz_final = pofz_final*post[ii].pofz
    pofz_final = alog(pofz_final/total(pofz_final))
;   pofz_final = pofz_final/total(pofz_final)
;   pofz_final = pofz_final/im_integral(pp.redshift,pofz_final)

    zmin = isedfit_find_zmin(pp.redshift,-2.0*alog(pofz_final>1D-20),nmin=nmin)
    struct_print, zmin
    
    xrange = [0,12]
;   yrange = [0.001,max(pofz_final)*1.01]
    yrange = [-15,0]
;   yrange = [0,max(pofz_final)*1.05]
;   yrange = [0,max(post.pofz)*1.05]
    
    psfile = isedfit_dir+'z10_a2744_pofz.eps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, height=4.5, $
      xmargin=[1.4,0.4], width=6.7
    djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5, $
      xrange=xrange, yrange=yrange;, /ylog
    
;   for ii = ngal-1, 0, -1 do begin
;   for ii = 0, ngal-1 do begin
    for ii = 0, ngal-1 do begin
       polyfill, [pp.redshift,reverse(pp.redshift)], $
         [post[ii].pofz,post[ii].pofz*0+!y.crange[0]], $
;        [post[ii].pofz,post[ii].pofz*0], $
         /fill, color=cgcolor(color[ii]), noclip=0
    endfor

;; redraw JDC
;    ii = 2
;    polyfill, [pp.redshift,reverse(pp.redshift)], $
;      [post[ii].pofz,post[ii].pofz*0], $
;      /fill, color=cgcolor(color[ii]), noclip=0

    for ii = 0, ngal-1 do djs_oplot, pp.redshift, post[ii].pofz, $
      thick=lthick[ii], line=0, color=cgcolor(lcolor[ii]);, psym=10
;   djs_oplot, pp.redshift, pofz_final, line=0, thick=8, $
;     color=cgcolor('black')
    
    djs_plot, [0], [0], /nodata, position=pos, /noerase, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $ ; /ylog, $
      ytitle='log (Likelihood)', xtitle='Redshift';, $
;     ytickname=replicate(' ',10)

    im_legend, cat.galaxy, /left, /top, box=0, charsize=1.5, $
      margin=0, color=color, line=0, thick=8, pspacing=1.7
;     margin=0, color=[color,'black'], line=0, thick=8, pspacing=1.7
;   im_legend, ['JD1A','JD1B','JD1A x JD1B'], /left, /top, box=0, $
;     margin=0, color=[color,'black'], line=0, thick=8, pspacing=1.7
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

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
         'z_{BPZ}='+strtrim(string(cat[ii].z_b_1,format='(F12.3)'),2)
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='Redshift', ytitle='Posterior Probability', $
         xrange=xrange, yrange=yrange, title=title
       djs_oplot, pp.redshift, pofz_ised, line=0, thick=8, psym=10, color=cgcolor('firebrick')
       djs_oplot, bpz_redshift, pofz_bpz, line=0, color=cgcolor('dodger blue'), thick=8, psym=10
       djs_oplot, rr[ii].z*[1,1], !y.crange, thick=6
       djs_oplot, cat[ii].z_b_1*[1,1], !y.crange, thick=6, line=5
       djs_oplot, cat[ii].bpz*[1,1], !y.crange, thick=6, line=5, color='orange'
;      im_legend, [galaxy[ii],'z_{iSEDfit}='+strtrim(string(rr[ii].z,format='(F12.3)'),2),$
;        'z_{BPZ}='+strtrim(string(zz[ii],format='(F12.3)'),2)], /left, /top, $
;        box=0, margin=0, textcolor=cgcolor(['black','firebrick','dodger blue'])
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
; --------------------------------------------------
; paper plot: SED + P(z)
    isedfit_paramfile = isedfit_dir+'z10_a2744_paramfile.par'

    cat = read_z10_a2744()
;   cat.mu = 1
    
    rr = read_isedfit(isedfit_paramfile,index=these,/getmodels,$
      isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir)

;   ul_color = ['powder blue','tomato','tan']
    mag_color = ['navy','firebrick','brown']
    sed_color = ['navy','firebrick','brown']
    ul_color = ['navy','firebrick','brown']
    mag_psym = [16,15,14]
    
    psfile = isedfit_dir+'z10_a2744_isedfit.ps'
    xrange = [0.3,7.0]
    im_plotconfig, 0, pos, psfile=psfile, height=5.0

    ticks = loglevels(xrange)
    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=[33,23], $
      xsty=1, ysty=1, /xlog, position=pos, $
      ytickname=ytickname, xtickname=xtickname, xtickv=ticks, $
      xticks=n_elements(ticks)-1, ytickinterval=2, $
      xtitle='Observed-frame Wavelength (\mu'+'m)', $
      ytitle='AB Magnitude (demagnified)'
    for ii = 0, 2 do djs_oplot, rr[ii].wave/1D4, rr[ii].flux+2.5*alog10(cat[ii].mu), $
      psym=10, color=cgcolor(sed_color[ii])

    im_legend, cat[0:2].galaxy, /left, /top, box=0, charsize=1.5, $
      margin=0, color=mag_color, pspacing=1.7, $
      psym=mag_psym
    
;; A
;    notzero = where(rr[0].bestmaggies gt 0.0)
;    bestmab = -2.5*alog10(rr[0].bestmaggies[notzero])
;    djs_oplot, weff[notzero], bestmab, psym=symcat(6,thick=6), symsize=2.5, $
;      color=cgcolor('grey')
;
;; B    
;    notzero = where(rr[1].bestmaggies gt 0.0)
;    bestmab = -2.5*alog10(rr[1].bestmaggies[notzero])
;    djs_oplot, weff[notzero], bestmab, psym=symcat(6,thick=6), symsize=2.5, $
;      color=cgcolor('orange')

; overplot the data
    for ii = 0, 2 do begin
       print, rr[ii].maggies*sqrt(rr[ii].ivarmaggies)
       nsigma = rr[ii].maggies*0+2.0
;      if ii eq 1 then nsigma[1] = 2.49
;      if ii eq 3 then nsigma[2] = 2.1
       mab = maggies2mag(rr[ii].maggies,nsigma=nsigma,$
         ivar=rr[ii].ivarmaggies,magerr=maberr,$
         lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
       used = where(mab gt -90.0,nused)
       upper = where(mab lt -90.0 and mabupper gt -90,nupper)
       if (nused ne 0L) then begin
          oploterror, weff[used], mab[used]+2.5*alog10(cat[ii].mu), hwhm[used]*0, $
            mabhierr[used], psym=symcat(mag_psym[ii]), $
            symsize=1.7, color=cgcolor(mag_color[ii]), /hibar, $
            errcolor=cgcolor(mag_color[ii]), errthick=!p.thick
          oploterror, weff[used], mab[used]+2.5*alog10(cat[ii].mu), hwhm[used]*0, $
            mabloerr[used], psym=3, color=cgcolor(mag_color[ii]), /lobar, $
            errcolor=cgcolor(mag_color[ii]), errthick=!p.thick
       endif
       if (nupper ne 0) then begin
          djs_oplot, weff[upper], mabupper[upper], $
            psym=symcat(11,thick=6), symsize=3.0, color=im_color(ul_color[ii])
       endif
    endfor

    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    

return
end

;; --------------------------------------------------
;; paper plot: SED + P(z)
;    isedfit_paramfile = isedfit_dir+'z10_a2744_paramfile.par'
;    isedfit_paramfile_photoz = isedfit_dir+'z10_a2744_photoz_paramfile.par'
;
;    cat = read_z10_a2744()
;    these = [0,1,2,n_elements(cat)-1]
;    cat = cat[these]
;    nobj = n_elements(cat)
;
;    common com_z10_a2744, rr, rr_photoz, post, pp, post1
;    if n_elements(rr) eq 0L then begin
;       rr = read_isedfit(isedfit_paramfile,index=these,/getmodels,$
;         isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir)
;       rr_photoz = read_isedfit(isedfit_paramfile_photoz,index=these,$
;         isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir,$
;         params=pp,isedfit_post=post)
;       post1 = post
;    endif
;    col = 'orange'
;
;; print out the physical properties
;    
;    
;    
;; pre-process the photoz's
;    bpz = rsex(isedfit_dir+'bpz.out')
;    pbpz = read_bpz_probs(isedfit_dir+'bpz.probs',redshift=bpz_zz)
;    pbpz = pbpz[these]
;    bpz = bpz[these]
;
;    for ii = 0, nobj-1 do begin
;       zmin = isedfit_find_zmin(bpz_zz,-2.0*alog(pbpz[ii].pofz>1D-20))
;       struct_print, zmin & print & help, bpz[ii], /str
;    endfor
;
;; make the plot    
;    xrange = [0.3,8.0]
;    yrange = [36,23]
;
;    psfile = isedfit_dir+'z10_a2744_isedfit.eps'
;    im_plotconfig, 5, pos, psfile=psfile, xmargin=[1.1,0.2], $
;      height=[2.6,2.6], xspace=0.05, yspace=0.05
;    for ii = 0, nobj-1 do begin
;       if odd(ii) then ytickname = replicate(' ',10) else delvarx, ytickname
;       if ii le 1 then xtickname = replicate(' ',10) else delvarx, xtickname
;       ticks = loglevels(xrange)
;       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
;         xsty=1, ysty=1, /xlog, position=pos[*,ii], noerase=ii gt 0, $
;         ytickname=ytickname, xtickname=xtickname, xtickv=ticks, $
;         xticks=n_elements(ticks)-1, ytickinterval=4
;;      djs_oplot, rr[ii].wave/1D4, rr[ii].flux
;       mask = (rr[ii].wave gt 1214.5*(1+rr[ii].z) and $
;         rr[ii].wave lt 1215*1.0025*(1+rr[ii].z)) eq 1
;       djs_oplot, rr[ii].wave/1D4, djs_maskinterp(rr[ii].flux,mask,rr[ii].wave/1D4)
;; model photometry
;       notzero = where(rr[ii].bestmaggies gt 0.0)
;       bestmab = -2.5*alog10(rr[ii].bestmaggies[notzero])
;       djs_oplot, weff[notzero], bestmab, psym=symcat(6,thick=6), symsize=2.5, $
;         color=cgcolor('grey')
;; overplot the data; distinguish between three different cases, based
;; on the input photometry
;       print, rr[ii].maggies*sqrt(rr[ii].ivarmaggies)
;       nsigma = rr[ii].maggies*0+2.0
;;      if ii eq 1 then nsigma[1] = 2.49
;;      if ii eq 3 then nsigma[2] = 2.1
;       mab = maggies2mag(rr[ii].maggies,nsigma=nsigma,$
;         ivar=rr[ii].ivarmaggies,magerr=maberr,$
;         lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
;       used = where(mab gt -90.0,nused)
;       upper = where(mab lt -90.0 and mabupper gt -90,nupper)
;       if (nused ne 0L) then begin
;          oploterror, weff[used], mab[used], hwhm[used], $
;            mabhierr[used], psym=symcat(16), $
;            symsize=1.7, color=cgcolor(col), /hibar, $
;            errcolor=cgcolor(col), errthick=3
;          oploterror, weff[used], mab[used], hwhm[used], $
;            mabloerr[used], psym=3, color=cgcolor(col), /lobar, $
;            errcolor=cgcolor(col), errthick=3
;       endif
;       if (nupper ne 0) then begin
;;         th = 6
;;         for uu = 0, nupper-1 do begin
;;            cgarrow, weff[upper[uu]], mabupper[upper[uu]], weff[upper[uu]], $
;;              mabupper[upper[uu]]+2, /data, /solid, hsize=!d.x_size/50, $
;;              thick=th, color=cgcolor('forest green')
;;            djs_oplot, weff[upper[uu]]+*[-1,1], thick=th, $
;;              mabupper[upper[uu]]*[1,1], color=cgcolor('forest green')
;;         endfor
;          djs_oplot, weff[upper], mabupper[upper], $
;            psym=symcat(11,thick=6), symsize=2, color=im_color('forest green')
;       endif
;       im_legend, /left, /top, box=0, spacing=1.5, charsize=1.5, margin=0, $
;         strtrim(cat[ii].galaxy,2)
;
;; P(z) inset
;       ipofz = post1[ii].pofz/im_integral(pp.redshift,post1[ii].pofz)
;       bpofz = pbpz[ii].pofz/im_integral(bpz_zz,pbpz[ii].pofz)
;       yr = [0,max(ipofz)>max(bpofz)]*1.02
;
;;      pos2 = [pos[0,ii]+0.05,pos[1,ii]+0.06,pos[2,ii]-0.02,pos[1,ii]+0.18]
;       pos2 = [pos[2,ii]-0.18,pos[1,ii]+0.06,pos[2,ii]-0.02,pos[1,ii]+0.18]
;       djs_plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, /norm, $
;         position=pos2, xrange=[0,12], yrange=yr
;
;       polyfill, [pp.redshift,reverse(pp.redshift)], [ipofz,ipofz*0+!y.crange[0]], $
;         /fill, color=cgcolor('tomato'), noclip=0
;       polyfill, [bpz_zz,reverse(bpz_zz)], [bpofz,bpofz*0+!y.crange[0]], $
;         /fill, color=cgcolor('dodger blue'), noclip=0
;
;       djs_oplot, pp.redshift, ipofz, line=0, color=cgcolor('firebrick')
;       djs_oplot, bpz_zz, bpofz, line=0, color=cgcolor('navy')
;       
;       djs_plot, [0], [0], /nodata, /noerase, xsty=9, ysty=9, /norm, $
;         position=pos2, xtitle='Redshift', ytitle='', charsize=0.8, $
;         xrange=[0,12], yrange=yr, xthick=3, ythick=3, yticks=1, $
;         ytickname=replicate(' ',10), charthick=3.0
;       xyouts, pos2[0,0]-0.01, djs_mean([pos2[1,0],pos2[3,0]]), 'Probability', $
;         align=0.5, orientation=90, charsize=0.8, charthick=3.0, /norm
;
;       im_legend, ['z_{iSEDfit}='+strtrim(string(rr_photoz[ii].z,format='(F12.2)'),2),$
;         'z_{BPZ}='+strtrim(string(bpz[ii].z_b,format='(F12.2)'),2)], $
;         /left, /top, textcolor=cgcolor(['tomato','dodger blue']), $
;         color=cgcolor(['tomato','dodger blue']), charsize=1.0, $
;         pspacing=1.0, box=0, charthick=3.0, margin=-1
;       
;    endfor 
;    
;    xyouts, pos[0,0]-0.08, pos[1,0], 'AB Magnitude', align=0.5, orientation=90, /normal
;    xyouts, pos[0,3], pos[1,3]-0.11, textoidl('Observed-frame Wavelength (\mu'+'m)'), $
;      align=0.5, /normal
;
;    im_plotconfig, psfile=psfile, /psclose, /pdf
;
;stop    
;    
