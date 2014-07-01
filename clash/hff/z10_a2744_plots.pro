pro z10_a2744_plots

    isedfit_dir = getenv('CLASH_PROJECTS')+'/hff/z10_a2744/'
    montegrids_dir = isedfit_dir+'montegrids/'

    filt = hff_filterlist(pivotwave=weff,width=hwhm,/useirac)
    weff = weff/1D4
    hwhm = hwhm/1D4

; --------------------------------------------------
; paper plot: P(z)'s
    isedfit_paramfile = isedfit_dir+'z10_a2744_photoz_paramfile.par'

    rr = read_isedfit(isedfit_paramfile,params=pp,isedfit_post=post,index=index)
    cat = read_z10_a2744(photoz=photoz)
    ngal = n_elements(cat)
    color = ['orange','tomato','powder blue','orchid']
    lcolor = ['red','firebrick','blue','black']
    lthick = [6,6,6,8]
;   color = ['orange','tomato','dodger blue','orchid','tan']
;   lcolor = ['red','firebrick','navy','purple','brown']
    
    for ii = 0, ngal-1 do post[ii].pofz = post[ii].pofz/total(post[ii].pofz)
;   for ii = 0, ngal-1 do post[ii].pofz = post[ii].pofz/$
;     im_integral(pp.redshift,post[ii].pofz)
    pofz_final = post[0].pofz
    for ii = 1, ngal-2 do pofz_final = pofz_final*post[ii].pofz
    pofz_final = pofz_final/total(pofz_final)
;   pofz_final = pofz_final/im_integral(pp.redshift,pofz_final)

    zmin = isedfit_find_zmin(pp.redshift,-2.0*alog(pofz_final>1D-20),nmin=nmin)
    struct_print, zmin
    
    xrange = [0,12]
    yrange = [0.01,max(pofz_final)*1.01]
;   yrange = [0,max(pofz_final)*1.05]
;   yrange = [0,max(post.pofz)*1.05]
    
    psfile = isedfit_dir+'z10_a2744_pofz.eps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, height=4.5, $
      xmargin=[1.3,0.4], width=6.8
    djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5, $
      xrange=xrange, yrange=yrange;, /ylog
    
;   for ii = ngal-1, 0, -1 do begin
;   for ii = 0, ngal-1 do begin
    for ii = 0, ngal-1 do begin
       polyfill, [pp.redshift,reverse(pp.redshift)], $
         [post[ii].pofz,post[ii].pofz*0], $
         /fill, color=cgcolor(color[ii]), noclip=0
    endfor
; redraw JDC
    ii = 2
    polyfill, [pp.redshift,reverse(pp.redshift)], $
      [post[ii].pofz,post[ii].pofz*0], $
      /fill, color=cgcolor(color[ii]), noclip=0

    for ii = 0, ngal-1 do djs_oplot, pp.redshift, post[ii].pofz, $
      thick=lthick[ii], line=0, color=cgcolor(lcolor[ii]);, psym=10
;   djs_oplot, pp.redshift, pofz_final, line=0, thick=8, $
;     color=cgcolor('black')
    
    djs_plot, [0], [0], /nodata, position=pos, /noerase, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $ ; /ylog, $
      ytitle='Redshift Probability', xtitle='Redshift', $
      ytickname=replicate(' ',10)

    im_legend, cat.galaxy, /left, /top, box=0, charsize=1.5, $
      margin=0, color=color, line=0, thick=8, pspacing=1.7
;     margin=0, color=[color,'black'], line=0, thick=8, pspacing=1.7
;   im_legend, ['JD1A','JD1B','JD1A x JD1B'], /left, /top, box=0, $
;     margin=0, color=[color,'black'], line=0, thick=8, pspacing=1.7
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
; --------------------------------------------------
; paper SED + P(z) plot
    isedfit_paramfile = isedfit_dir+'z10_a2744_paramfile.par'

    cat = read_z10_a2744()
    rr = read_isedfit(isedfit_paramfile,index=these,/getmodels,$
      isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir)

;   ul_color = ['powder blue','tomato','tan']
    mag_color = ['navy','firebrick','brown']
    sed_color = ['navy','firebrick','brown']
    ul_color = ['navy','firebrick','brown']
    mag_psym = [16,15,14]
    
    psfile = isedfit_dir+'z10_a2744_isedfit.ps'
    xrange = [0.3,6.0]
    im_plotconfig, 0, pos, psfile=psfile, height=5.0

    ticks = loglevels(xrange)
    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=[33,25], $
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
    

return
end

