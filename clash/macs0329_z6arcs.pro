pro macs0329_z6arcs, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber
; jm11nov08ucsd - build plots for the paper

    common com_post, model, ised, mstar, post, age, Z, tau, sfr0, b100, av, sfrage
    
    isedpath = clash_path(/ised)
    datapath = clash_path(/macs0329_z6arcs)

    isedfit_sfhgrid_dir = clash_path(/monte)
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/clash_sfhgrid.par'

    filters = clash_filterlist(nice=nice_filters,/useirac,weff=weff)
    ndraw = isedfit_ndraw()

    zobj = 6.18
    prefix = 'macs0329_z6arcs'
    paramfile = isedpath+prefix+'_supergrid03_isedfit.par'
;   paramfile = isedpath+prefix+'_supergrid02_isedfit.par'
    
; gather the photometry and restore the results
    cat = read_z6arcs_photometry(adi=adi)
    galaxy = repstr(adi.id,'arc_','Arc ')
    nobj = 4

    if (n_elements(model) eq 0) then $
      model = isedfit_restore(paramfile,ised,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)

    if (n_elements(mstar) eq 0L) then begin
       mstar = isedfit_reconstruct_posterior(paramfile,post=post,$
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath,$
         age=age,Z=Z,tau=tau,sfr0=sfr0,b100=b100,av=av,sfrage=sfrage)
    endif
    ssfr = sfr0-mstar+9 
    
;   niceprint, galaxy, ised.mass_50-alog10(adi.mu), ised.mass_err, $
;     1E3*ised.age_50, 1E3*ised.age_err, ised.z_50/0.019, ised.z_err/0.019, ised.av_50, $
;     ised.av_err, 10^(ised.sfr_50-alog10(adi.mu)), ised.sfr_err*10^(ised.sfr_50-alog10(adi.mu))*alog(10), $
;     1E3*ised.tau_50, 1E3*ised.tau_err

; write out a Table for Adi
    print, 'Arc, mass (log Msun, demagnified), <age> (Myr), metallicity '+$
      '(with respect to solar=0.019), AV (mag), SFR (Msun/yr, demagnified), sSFR (Gyr^-1), formation redshift'
    for ii = 0, nobj-1 do print, galaxy[ii], ised[ii].mass_50-alog10(adi[ii].mu), 1E3*median(sfrage[*,ii]), $
      ised[ii].z_50/0.019, ised[ii].av_50, 10^(ised[ii].sfr_50-alog10(adi[ii].mu)), 10^median(ssfr[*,ii]), $
      getredshift(getage(zobj)-median(sfrage[*,ii]));, format='(A8,2x,F6.3,2x,I0,2x,F7.3,2x,F8.3)'

; ---------------------------------------------------------------------------
; make a plot for the paper    
    this = 1

    readcol, datapath+'lephare_combined_p_of_z.txt', lpzz, lppdf, /silent
;   bpz = rsex(datapath+'Pz.dat')
;   zz = bpz.z & pz = bpz.pcombined
;   readcol, datapath+'z6p18.dat', zz, pz, /silent
    readcol, datapath+'z6p18_v2.dat', zz, pz, /silent
    pz = pz/max(pz)
    lppdf = lppdf/max(lppdf)
    
    psfile = datapath+'z6arc_sed.eps'
    im_plotconfig, 0, pos, psfile=psfile, ymargin=[0.8,5.2], height=4.0

    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
    ytitle1 = textoidl('Apparent AB Magnitude')
;   ytitle1 = textoidl('m_{AB}')

    yrange = [28.0,23.0]
    xrange1 = [2400,60000]
;   xrange1 = [1500,18000]
    xrange2 = xrange1/(1.0+ised[this].zobj)

;   ticks1 = loglevels(xrange1)
    ticks1 = [4000,12000,30000]
    ticks2 = [500,1000,2000,5000]

    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $
      position=pos, xtickformat='(I0)', xticks=n_elements(ticks1)-1, xtickv=ticks1
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2
    im_legend, ['Arc 1.2','z_{phot} = 6.18^{+0.05}_{-0.07}'], /left, /top, box=0, margin=0, charsize=1.7
    oplot, model[this].wave, model[this].flux, line=0
    
; overplot the observed and model photometry
    used = where(cat[this].limit eq 0 and finite(cat[this].mag) and strmatch(filters,'*f606w*',/fold) eq 0,nused)
    upper = where(cat[this].limit eq 1 and finite(cat[this].mag),nupper)
    special = where(cat[this].limit eq 0 and finite(cat[this].mag) and strmatch(filters,'*f606w*',/fold) eq 1,nspecial)

    djs_oplot, cat[this].weff, -2.5*alog10(ised[this].bestmaggies), $
      psym=symcat(6,thick=6), symsize=2.5
;   djs_oplot, cat[this].weff[used], -2.5*alog10(ised[this].bestmaggies[used]), $
;     psym=symcat(6,thick=6), symsize=2.5

    if (nused ne 0L) then begin
       oploterror, cat[this].weff[used], cat[this].mag[used], cat[this].magerr_lo[used], $
         psym=3, color=im_color('firebrick'), /lobar, $
         errcolor=im_color('firebrick'), errthick=!p.thick
       oploterror, cat[this].weff[used], cat[this].mag[used], cat[this].magerr_up[used], $
         psym=symcat(16), symsize=2.0, color=im_color('firebrick'), $
         errcolor=im_color('firebrick'), errthick=!p.thick, /hibar
    endif

    if (nspecial ne 0L) then begin
       oploterror, cat[this].weff[special], cat[this].mag[special], cat[this].magerr_lo[special], $
         psym=symcat(9,thick=6), symsize=2.0, color=im_color('firebrick'), /lobar, $
         errcolor=im_color('firebrick'), errthick=6
       oploterror, cat[this].weff[special], cat[this].mag[special], cat[this].magerr_up[special], $
         psym=symcat(9,thick=6), symsize=2.0, color=im_color('firebrick'), $
         errcolor=im_color('firebrick'), errthick=6, /hibar
    endif

    if (nupper ne 0) then begin
       djs_oplot, cat[this].weff[upper], cat[this].mag[upper], $ ; cat[this].magerr_up[upper], $
         psym=symcat(11,thick=8), symsize=3.0, color=im_color('steel blue')
    endif

; P(z) inset
    pp = [0.6,0.6,0.9,0.77]
    plot, [0], [0], xsty=1, ysty=1, /nodata, /noerase, $
      yrange=[0,1.05], xrange=[5.6,6.601], position=pp, $
      ytitle='Probability', xtitle='Redshift', xtickinterval=0.3, $
      charsize=1.2, ytickname=replicate(' ',10)
;   plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, $
;     position=pp, yrange=[0,1.05], $
;     xrange=[zobj-1,7.0]
;   polyfill, [zz,reverse(zz)], [pz*0,reverse(pz)], /fill, $
;     noclip=0, /data, color=im_color('grey80')
    legend, ['BPZ','LPZ'], /left, /top, box=0, margin=0, $
      line=[0,5], thick=4, pspacing=1.7, charsize=1.2
    djs_oplot, zz, pz, thick=6, line=0;, color=im_color('midnight blue')
    djs_oplot, lpzz, lppdf, thick=6, line=5;, color=im_color('forest green')
;   plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.05], $
;     xrange=[zobj-1,7.0], position=pp, ytitle='Probability', $
;     xtitle='Redshift', xtickinterval=0.5, ytickname=replicate(' ',10), $
;     charsize=1.3
    
;   im_plotconfig, /psclose

; posterior distributions
;   psfile = 'z6arc_post.eps'

    im_plotconfig, 17, pos, xspace=[0.3,0.3], yspace=1.0, ymargin=[8.2,1.0], $
      xmargin=[1.2,0.4], height=2.1*[1,1], width=2.7*[1,1,1]

    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, position=pos[*,0], yrange=[0,1.05], $
      xrange=[7.0,10.3], ytitle='', charsize=1.4, xtitle=''
    im_plothist, mstar[*,this]-alog10(adi[this].mu), bin=0.15, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, (ised[this].mass_50-alog10(adi[this].mu))*[1,1], !y.crange, $
      line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, position=pos[*,0], yrange=[0,1.05], $
      xrange=[7.0,10.3], ytitle='Probability', charsize=1.4, $
      xtitle=textoidl('log (M_{*}/M'+sunsymbol()+')'), $
      xtickinterval=1, ytickname=replicate(' ',10)
    
    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.05], $
      xrange=[-0.04,0.92], position=pos[*,1], ytitle='', charsize=1.4, $
      xtitle=''
    im_plothist, sfrage[*,this], bin=0.05, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, median(sfrage[*,this])*[1,1], !y.crange, line=5, thick=8
;   im_plothist, age[*,this], bin=0.05, /peak, /overplot, $
;     /fill, fcolor=im_color('grey80')
;   djs_oplot, ised[this].age_50*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.05], $
      xrange=[-0.04,0.92], position=pos[*,1], ytitle='', charsize=1.4, $
      xtitle='<Age> (Gyr)', ytickname=replicate(' ',10), $
      xtickinterval=0.2

    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[-0.04,1.6], position=pos[*,2], ytitle='', $
      xtitle='', charsize=1.4
    im_plothist, Z[*,this]/0.02, bin=0.1, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, ised[this].Z_50/0.02*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[-0.04,1.6], position=pos[*,2], ytitle='', $
      xtitle='Z/Z'+sunsymbol(), ytickname=replicate(' ',10), charsize=1.4, $
      xtickinterval=0.5

    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[-0.03,0.62], position=pos[*,3], ytitle='', $
      xtitle=textoidl('A_{V} (mag)'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=0.2
    im_plothist, av[*,this], bin=0.04, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, ised[this].av_50*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[-0.03,0.62], position=pos[*,3], ytitle='Probability', $
      xtitle=textoidl('A_{V} (mag)'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=0.2

    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[1,11], position=pos[*,4], ytitle='', $
      xtitle=textoidl('SFR (M'+sunsymbol()+' yr^{-1})'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=2
    im_plothist, 10^(sfr0[*,this]-alog10(adi[this].mu)), bin=0.5, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, 10^(ised[this].sfr_50-alog10(adi[this].mu))*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[1,11], position=pos[*,4], ytitle='', $
      xtitle=textoidl('SFR (M'+sunsymbol()+' yr^{-1})'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=2

    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[0.0,10.0], position=pos[*,5], ytitle='', $
      ytickname=replicate(' ',10), charsize=1.4
    im_plothist, 10^ssfr[*,this], bin=0.6, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, median(10^ssfr[*,this])*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[0.0,10.0], position=pos[*,5], ytitle='', $
      xtitle=textoidl('sSFR (Gyr^{-1})'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=2

;   plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.1], $
;     xrange=[-0.05,5.1], position=pos[*,5], ytitle='', $
;     xtitle=textoidl('\tau (Gyr)'), $
;     ytickname=replicate(' ',10), charsize=1.4, xtickinterval=0.5
;   im_plothist, tau[*,this], bin=0.3, /peak, /overplot, $
;     /fill, fcolor=im_color('grey80')
;   djs_oplot, ised[this].tau_50*[1,1], !y.crange, line=5, thick=8
;   plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
;     xrange=[-0.05,5.1], position=pos[*,5], ytitle='', $
;     xtitle=textoidl('\tau (Gyr)'), $
;     ytickname=replicate(' ',10), charsize=1.4, xtickinterval=1

    im_plotconfig, /psclose, psfile=psfile, /pskeep, /pdf

stop    
    
    
; ---------------------------------------------------------------------------
; plot the NxN distributions of parameters as an upper triangle
    for ii = 0, n_elements(cat)-1 do begin
       psfile = datapath+strlowcase(adi[ii].id)+'_manyd.ps'
       splog, 'Writing '+psfile
       manyd = transpose([[mstar[*,ii]],[sfrage[*,ii]],[Z[*,ii]/0.019],[av[*,ii]],[sfr0[*,ii]]])
       label = textoidl(['log (M/M'+sunsymbol()+')','Age_{w} (Gyr)','Z/Z'+sunsymbol(),'A_{V} (mag)','log (\psi)'])

       im_manyd_scatterplot, fltarr(ndraw)+1, manyd, psfile, label=label, $
         axis_char_scale=1.4, /internal, outliers=1, $
         /nogrey, levels=errorf((dindgen(2)+1)/sqrt(2)), /upper, $
         title=galaxy[ii]
       spawn, 'ps2pdf '+psfile, /sh
    endfor       

; merge everything into a single postscript file and then convert to PDF
    outfile = datapath+'z6arc_manyd.pdf'
    allfiles = [datapath+strlowcase(adi.id)+'_manyd.ps']
    splog, 'Writing '+outfile
    spawn, 'gs -q -dNOPAUSE -sDEVICE=pdfwrite '+$
      '-sOutputFile='+outfile+' -dBATCH '+strjoin(allfiles,' '), /sh

stop    
    
; ---------------------------------------------------------------------------
; QAplot comparing the demagnified photometry of all the arcs
    psfile = datapath+'z6arc_demagnified.ps'
    im_plotconfig, 0, pos, psfile=psfile, ymargin=[1.0,1.1], height=5.0

    cat = read_z6arcs_photometry(adi=adi)
    
    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
    ytitle1 = textoidl('Demagnified AB Magnitude')

    ticks1 = [4000,12000,30000]
    ticks2 = [500,1000,2000,5000]
    
    yrange = [30.5,25.5]
;   yrange = [32,20.0]
    xrange1 = [2000,70000]
    xlog = 1

    xrange2 = xrange1/(1.0+zobj)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=xlog, $
      position=pos, xtickformat='(I0)', $
      xticks=n_elements(ticks1)-1, xtickv=ticks1
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=xlog, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2

    color = ['dodger blue','firebrick','forest green','black']
    for ii = 0, nobj-1 do begin
       used = where(cat[ii].limit eq 0 and finite(cat[ii].mag),nused)
       upper = where(cat[ii].limit eq 1 and finite(cat[ii].mag),nupper)

       if (nused ne 0) then begin
          oploterror, cat[ii].weff[used], cat[ii].mag[used]+2.5*alog10(adi[ii].mu), cat[ii].magerr[used], $
            psym=symcat(16,thick=4), symsize=2.0, color=im_color(color[ii]), $
            errcolor=im_color(color[ii]), errthick=!p.thick
       endif
       if (nupper ne 0) then begin
          djs_oplot, cat[ii].weff[upper], cat[ii].mag[upper]+2.5*alog10(adi[ii].mu), $
            psym=symcat(11,thick=4), symsize=2.0, color=im_color(color[ii])
       endif
    endfor
    im_legend, textoidl(galaxy), /right, /bottom, box=0, charsize=1.7, $
      margin=0, spacing=2.2, textcolor=color
    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

stop    
    
; ---------------------------------------------------------------------------
; QAplot comparing the SE and Postman photometry
    cat1 = read_macs0329_z6arcs(adi=adi,/usemag) ; sextractor
    cat2 = read_macs0329_z6arcs(adi=adi) ; marc

    psfile = datapath+'qaplot_z6arc_photometry.ps'
    im_plotconfig, 0, pos, psfile=psfile, ymargin=[1.0,1.1], height=5.0

    clash_to_maggies, cat1, maggies1, ivarmaggies1, /usemag, /useirac
    clash_to_maggies, cat2, maggies2, ivarmaggies2, /useirac
    
    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
    ytitle1 = textoidl('m_{AB}')

    yrange = [32,20.0]
    xrange1 = [2000,70000]
    xlog = 1
    for ii = 0, nobj-1 do begin
       xrange2 = xrange1/(1.0+zobj)
       plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
         xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=xlog, $
         xtickinterval=3000, position=pos, xtickformat='(I0)'
       axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=xlog

; SExtractor       
       used = where((maggies1[*,ii] gt 0.0) and $
         (ivarmaggies1[*,ii] gt 0.0),nused)
       upper = where((maggies1[*,ii] le 0.0) and $
         (ivarmaggies1[*,ii] gt 0.0),nupper)

       if (nused ne 0) then begin
          mab = maggies2mag(maggies1[used,ii],ivar=ivarmaggies1[used,ii],magerr=mab_err)
          oploterror, weff[used], mab, mab_err, psym=symcat(16,thick=4), $
            symsize=2.0, color=im_color('forest green'), $
            errcolor=im_color('forest green'), errthick=!p.thick
       endif
       if (nupper ne 0) then begin
          mab = maggies2mag(1.0/sqrt(ivarmaggies1[upper,ii]))
          djs_oplot, weff[upper], mab, psym=symcat(11,thick=4), $
            symsize=2.0, color=im_color('forest green')
       endif
    
; Marc   
       used = where((maggies2[*,ii] gt 0.0) and $
         (ivarmaggies2[*,ii] gt 0.0),nused)
       upper = where((maggies2[*,ii] le 0.0) and $
         (ivarmaggies2[*,ii] gt 0.0),nupper)

       if (nused ne 0) then begin
          mab = maggies2mag(maggies2[used,ii],ivar=ivarmaggies2[used,ii],magerr=mab_err)
          oploterror, weff[used], mab, mab_err, psym=symcat(16,thick=4), $
            symsize=2.0, color=im_color('midnight blue'), $
            errcolor=im_color('midnight blue'), errthick=!p.thick
       endif
       if (nupper ne 0) then begin
          mab = maggies2mag(1.0/sqrt(ivarmaggies2[upper,ii]))
          djs_oplot, weff[upper], mab, psym=symcat(11,thick=4), $
            symsize=2.0, color=im_color('midnight blue')
       endif

       legend, textoidl(galaxy[ii]), /left, /top, box=0, charsize=1.7, $
         margin=0, spacing=2.2
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

; ---------------------------------------------------------------------------
; make a QAplot showing the best-fitting SEDs and posterior distributions
    isedfitpsfile = datapath+prefix+'.ps'
;   isedfitpsfile = datapath+prefix+'_isedfit.ps'
    im_plotconfig, 0, pos, psfile=isedfitpsfile, ymargin=[1.0,1.1], $
      height=5.0

    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
    ytitle1 = textoidl('m_{AB}')

    yrange = [31,21.5]
    xrange1 = [2000,70000]
    xlog = 1
;   yrange = [30.5,24.0]
;   xrange1 = [2000,17000]
    for ii = 0, nobj-1 do begin
       xrange2 = xrange1/(1.0+ised[ii].zobj)
       plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
         xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=xlog, $
         xtickinterval=3000, position=pos, xtickformat='(I0)'
       axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=xlog

       oplot, model[ii].wave, model[ii].flux, line=0
       
; overplot the observed and model photometry
       used = where((ised[ii].maggies gt 0.0) and $ ; used in the fitting
         (ised[ii].ivarmaggies gt 0.0),nused)
       notused = where((ised[ii].maggies gt 0.0) and $ ; not used in the fitting
         (ised[ii].ivarmaggies eq 0.0),nnotused)
       nodata = where((ised[ii].maggies eq 0.0) and $ ; no measurement
         (ised[ii].ivarmaggies eq 0.0),nnodata)

       djs_oplot, weff[used], -2.5*alog10(ised[ii].bestmaggies[used]), $
         psym=symcat(6,thick=6), symsize=2.5

       used = where((ised[ii].maggies gt 0.0) and $ ; used in the fitting
         (ised[ii].ivarmaggies gt 0.0),nused)
       upper = where((ised[ii].maggies le 0.0) and $ ; upper limit
         (ised[ii].ivarmaggies gt 0.0),nupper)

       if (nused ne 0L) then begin
          mab = maggies2mag(ised[ii].maggies[used],$
            ivar=ised[ii].ivarmaggies[used],magerr=mab_err)
          oploterror, weff[used], mab, mab_err, psym=symcat(16), $
            symsize=2.0, color=im_color('firebrick'), $
            errcolor=im_color('firebrick'), errthick=!p.thick
       endif

       if (nupper ne 0) then begin
          mab = maggies2mag(1.0/sqrt(ised[ii].ivarmaggies[upper]))
          oploterror, weff[upper], mab, mab*0.0, psym=symcat(11,thick=6), $
            symsize=3.0, color=im_color('steel blue'), $
            errcolor=im_color('steel blue'), errthick=!p.thick
       endif

;; overplot uvis, acs, and ir separately 
;       uvis = where(strmatch(nice_filters,'*uvis*',/fold))
;       mab = maggies2mag(ised[ii].maggies[uvis],$
;         ivar=ised[ii].ivarmaggies[uvis],magerr=mab_err)
;       oploterror, filtinfo[uvis].weff, mab, mab_err, psym=symcat(16), $
;         symsize=2.0, color=fsc_color('dodger blue',101), $
;         errcolor=fsc_color('dodger blue',101), errthick=!p.thick
;       
;       acs = where(strmatch(nice_filters,'*acs*',/fold))
;       mab = maggies2mag(ised[ii].maggies[acs],$
;         ivar=ised[ii].ivarmaggies[acs],magerr=mab_err)
;       oploterror, filtinfo[acs].weff, mab, mab_err, psym=symcat(16), $
;         symsize=2.0, color=fsc_color('tan',101), $
;         errcolor=fsc_color('tan',101), errthick=!p.thick
;       
;       ir = where(strmatch(nice_filters,'*ir*',/fold))
;       mab = maggies2mag(ised[ii].maggies[ir],$
;         ivar=ised[ii].ivarmaggies[ir],magerr=mab_err)
;       oploterror, filtinfo[ir].weff, mab, mab_err, psym=symcat(16), $
;         symsize=2.0, color=fsc_color('firebrick',101), $
;         errcolor=fsc_color('firebrick',101), errthick=!p.thick
    
       label = galaxy[ii]
;      label = [galaxy[ii],'z = '+string(ised[ii].zobj,format='(F6.4)')]
       legend, textoidl(label), /left, /top, box=0, charsize=1.7, $
         margin=0, spacing=2.2

       label = [$
         'log (M/M'+sunsymbol()+') = '+strtrim(string(ised[ii].mass_50-alog10(adi[ii].mu),format='(F12.2)'),2),$
         'Age = '+strtrim(string(ised[ii].age_50*1E3,format='(I0)'),2)+' Myr',$
         'Z/Z'+sunsymbol()+' = '+strtrim(string(ised[ii].Z_50/0.019,format='(F12.2)'),2),$
         'A_{V} = '+strtrim(string(ised[ii].av_50,format='(F12.2)'),2)+' mag',$
         '\tau = '+strtrim(string(ised[ii].tau_50,format='(F12.2)'),2)+' Gyr',$
         '\psi = '+strtrim(string(10^(ised[ii].sfr_50-alog10(adi[ii].mu)),$
         format='(F12.2)'),2)+' M'+sunsymbol()+' yr^{-1}']
;        'log (b_{100}) = '+strtrim(string(ised[ii].b100_50,format='(F12.3)'),2)]
       legend, textoidl(label), /left, /top, box=0, spacing=1.9, charsize=1.4, $
         position=[15000.0,yrange[0]-(yrange[0]-yrange[1])*0.44], /data
    endfor
    im_plotconfig, psfile=isedfitpsfile, /psclose, /pdf, /pskeep

; ---------------------------------------------------------------------------
; P(M) QAplot
    col = ['black','forest green','dodger blue','firebrick']
    line = [0,3,4,5]
    psfile = datapath+'qa_pofm.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.5
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=[6.8,11.2], yrange=[0,200], xtitle=mzplot_masstitle(), $
      ytitle='Probability', ytickname=replicate(' ',10)
    for ii = 0, 3 do im_plothist, mstar[*,ii]-alog10(adi[ii].mu), $
      line=line[ii], /overplot, thick=8, color=im_color(col[ii]), bin=0.1
    im_legend, galaxy, /right, /top, box=0, line=line, $
      color=col, charsize=1.5
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
stop    
    
return
end


;; ---------------------------------------------------------------------------
;; make a plot for the paper    
;    this = 1
;    
;    psfile = 'z6arc_sed.eps'
;    im_plotconfig, 0, pos, psfile=psfile, ymargin=[1.0,1.1], height=5.0
;
;    xtitle1 = textoidl('Observed Wavelength (\AA)')
;    xtitle2 = textoidl('Rest Wavelength (\AA)')
;    ytitle1 = textoidl('m_{AB}')
;
;    yrange = [29.0,23.6]
;    xrange1 = [1500,21000]
;    xrange2 = xrange1/(1.0+ised[this].zobj)
;    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
;      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=0, $
;      xtickinterval=5000, position=pos, xtickformat='(I0)'
;    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=0
;
;    oplot, model[this].wave, model[this].flux, line=0
;    
;; overplot the observed and model photometry
;    used = where((ised[this].maggies gt 0.0) and $ ; used in the fitting
;      (ised[this].ivarmaggies gt 0.0),nused)
;    notused = where((ised[this].maggies gt 0.0) and $ ; not used in the fitting
;      (ised[this].ivarmaggies eq 0.0),nnotused)
;    nodata = where((ised[this].maggies eq 0.0) and $ ; no measurement
;      (ised[this].ivarmaggies eq 0.0),nnodata)
;
;    djs_oplot, weff[used], -2.5*alog10(ised[this].bestmaggies[used]), $
;      psym=symcat(6,thick=6), symsize=2.5
;
;    used = where((ised[this].maggies gt 0.0) and $ ; used in the fitting
;      (ised[this].ivarmaggies gt 0.0),nused)
;    upper = where((ised[this].maggies le 0.0) and $ ; upper limit
;      (ised[this].ivarmaggies gt 0.0),nupper)
;
;    if (nused ne 0L) then begin
;       mab = maggies2mag(ised[this].maggies[used],$
;         ivar=ised[this].ivarmaggies[used],magerr=mab_err)
;       oploterror, weff[used], mab, hwhm[used], mab_err, psym=symcat(16), $
;         symsize=2.0, color=im_color('firebrick'), $
;         errcolor=im_color('firebrick'), errthick=!p.thick
;;      oploterror, weff[used], mab, hwhm[used], mab_err, psym=symcat(16), $
;;        symsize=2.0, color=im_color('firebrick'), $
;;        errcolor=im_color('firebrick'), errthick=!p.thick
;    endif
;
;    if (nupper ne 0) then begin
;       mab = maggies2mag(1.0/sqrt(ised[this].ivarmaggies[upper]))
;       djs_oplot, weff[upper], mab, psym=symcat(11,thick=8), $
;         symsize=3.0, color=im_color('steel blue')
;;      oploterror, weff[upper], mab, hwhm[upper], mab*0.0, psym=symcat(11,thick=6), $
;;        symsize=3.0, color=im_color('steel blue'), $
;;        errcolor=im_color('steel blue'), errthick=!p.thick
;    endif
;
;; inset with posterior distributions
;    x1 = 0.50
;    y1 = 0.50
;    y2 = 0.24
;    hgt = 0.16
;    wid = 0.20
;    fact = 1.1
;    H0 = '70'
;    
;    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.05], $
;      xrange=[7.5,9.9], position=[x1,y1,x1+wid,y1+hgt], ytitle='', $
;      xtitle=textoidl('log (M_{*}/M'+sunsymbol()+')'), $
;      xtickinterval=1, ytickname=replicate(' ',10), $
;      charsize=1.1
;    im_plothist, mstar[*,this]-alog10(adi[this].mu), bin=0.15, /peak, /overplot, $
;      /fill, fcolor=im_color('grey80')
;    djs_oplot, (ised[this].mass_50-alog10(adi[this].mu))*[1,1], !y.crange, $
;      line=5, thick=6
;    
;    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.05], $
;      xrange=[-0.02,0.9], position=[x1+wid*fact,y1,x1+wid*fact+wid,y1+hgt], ytitle='', $
;      xtitle='Age (Gyr)', ytickname=replicate(' ',10), charsize=1.1, $
;      xtickinterval=0.3
;    im_plothist, age[*,this], bin=0.05, /peak, /overplot, $
;      /fill, fcolor=im_color('grey80')
;    djs_oplot, ised[this].age_50*[1,1], !y.crange, line=5, thick=6
;
;    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
;      xrange=[-0.02,1.6], position=[x1,y2,x1+wid,y2+hgt], ytitle='', $
;      xtitle='Z/Z'+sunsymbol(), ytickname=replicate(' ',10), charsize=1.1, $
;      xtickinterval=0.5
;    im_plothist, Z[*,this]/0.02, bin=0.1, /peak, /overplot, $
;      /fill, fcolor=im_color('grey80')
;    djs_oplot, ised[this].Z_50/0.02*[1,1], !y.crange, line=5, thick=6
;
;    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
;      xrange=[-0.02,1.0], position=[x1+wid*fact,y2,x1+wid*fact+wid,y2+hgt], ytitle='', $
;      xtitle=textoidl('A_{V} (mag)'), $
;      ytickname=replicate(' ',10), charsize=1.1, xtickinterval=0.5
;    im_plothist, av[*,this], bin=0.1, /peak, /overplot, $
;      /fill, fcolor=im_color('grey80')
;    djs_oplot, ised[this].av_50*[1,1], !y.crange, line=5, thick=6
;
;;   plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
;;     xrange=[0,200], position=[x1+wid*fact,y2,x1+wid*fact+wid,y2+hgt], ytitle='', $
;;     xtitle=textoidl('\psi (M'+sunsymbol()+' yr^{-1})'), $
;;     ytickname=replicate(' ',10), charsize=1.1, xtickinterval=100
;;   im_plothist, 10^(sfr0[*,this]-alog10(adi[this].mu)), bin=10, /peak, /overplot, $
;;     /fill, fcolor=im_color('grey80')
;;   djs_oplot, 10^(ised[this].sfr_50-alog10(adi[this].mu))*[1,1], !y.crange, line=5, thick=6
;
;    im_plotconfig, /psclose
;
