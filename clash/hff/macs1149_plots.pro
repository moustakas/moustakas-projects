;;+
;; Calls ja_plot of ja_oplot but designed to do histograms
;; properly (not like psym=10. Especially useful with log axes
;; N.B. will loose top bin - set max to max(xval)+binsize to  avoid this
;; INPUTS:
;;  x - locations from histogram (i.e. the lower limit of the bin)
;;  y - histogram values
;;-
;pro ja_histplot,x,y,over=over,_extra=_extra
;
;  xr=reform(transpose(rebin(x,n_elements(x),2)),n_elements(x)*2)
;  yr=shift(reform(transpose(rebin(y,n_elements(y),2)),n_elements(y)*2),+1)
;  yr[0]=0
;  yr[n_elements(yr)-1]=0
;
;  if keyword_set(over) then ja_oplot,xr,yr,_extra=_extra $
;  else ja_plot,xr,yr,_extra=_extra
;end
;
;IDL> h=histogram(alog10(test),locations=l,binsize=0.5,max=4)
;IDL> ja_histplot,10^l,h,/xlog,yrange=[0,2]

function age2zform, axis, index, value
; convert from age to formation redshift
    zz = getredshift(getage(9.44)-value/1E3)
;   print, value, zz
    return, string(zz, format='(F0.1)')
end

pro render_postplot, xx, pos, xrange=xrange, yrange=yrange, $
  xtitle=xtitle, ytitle=ytitle, binsize=binsize, noerase=noerase, $
  xtickinterval=xtickinterval, nomedian=nomedian, monte=monte, $
  charsize=charsize, nonorm=nonorm, color_fill=color_fill, $
  color_outline=color_outline, color_monte=color_monte, $
  fill_monte=fill_monte, xlog=xlog, logbins=logbins, _extra=extra, $
  xstyle=xstyle

    yrange = [0,1.1]
;   yrange = [0,1.05]

    if (n_elements(xrange) eq 0) then xrange = minmax(xx)*[0.9,1.1]
    if (n_elements(binsize) eq 0) then begin
       if keyword_set(logbins) then $
         binsize = (alog10(xrange[1])-alog10(xrange[0]))/ceil(0.3*sqrt(n_elements(xx))) else $
           binsize = (xrange[1]-xrange[0])/ceil(0.3*sqrt(n_elements(xx)))
    endif

    if n_elements(color_fill) eq 0 then color_fill = 'tan' ; 'pale turquoise'
    if n_elements(color_outline) eq 0 then color_outline = 'tan' ; 'pale turquoise'
    if n_elements(color_monte) eq 0 then color_monte = 'grey80'

    plot, [0], [0], xsty=5, ysty=5, /nodata, position=pos, $
      yrange=yrange, xrange=xrange, noerase=noerase, charsize=charsize, $
      xlog=xlog
    im_plothist, xx, bin=binsize, /peak, /overplot, /fill, $
      fcolor=im_color(color_fill), xhist, yhist, charsize=charsize, $
      color=im_color(color_outline,255), logbins=logbins, xlog=xlog
;   if keyword_set(nomedian) eq 0 then $
;     djs_oplot, median(xx)*[1,1], !y.crange, line=5, thick=8
    if n_elements(monte) ne 0 then begin
       im_plothist, monte, bin=binsize*2, mxhist, myhist, /noplot, $
         _extra=extra, logbins=logbins, xlog=xlog
       if keyword_set(nonorm) then normfactor = 1D else $
         normfactor = max(myhist)/(yrange[1]*0.9)
       if keyword_set(fill_monte) then begin
          im_plothist, monte, bin=binsize*2, /overplot, /fill, $
            normfactor=normfactor, fcolor=im_color(color_monte), $
            color=im_color(color_monte), logbins=logbins, xlog=xlog
       endif else begin
          im_plothist, monte, bin=binsize*2, /overplot, line=1, $
            normfactor=normfactor, color=im_color(color_monte), $
            logbins=logbins, xlog=xlog
;           normfactor=max(myhist)/max(yhist)/1.2
       endelse
    endif
    if n_elements(xstyle) eq 0 then xstyle = 1
    plot, [0], [0], xsty=xstyle, ysty=1, /noerase, /nodata, position=pos, $
      yrange=yrange, xrange=xrange, ytitle=ytitle, $
      xtitle=textoidl(xtitle), xtickinterval=xtickinterval, $
      ytickname=replicate(' ',10), charsize=charsize, xlog=xlog, $
      _extra=extra
return
end

pro macs1149_plots
; jm16mar15siena - plots for the MACS1149 paper

; I plan to use your figure for M1149-JD (source 663), as the new data allow
; improved fitting. Most other candidates are faint, and the results of SED
; fitting may not significant. I intend to plot the one with fixed z=9.44, to
; avoid potential confusion as the difference is at ~ 2 sigma level. Could you
; please write some text about the parameters in the lower panel: The difference
; between blue and brown histograms.  Why the Age(SFR) is smaller than "Age"
; while the SFR rate is lower? I thought it is the other way around. The
; formation epoch would be 0.52 - 0.2=0.32 Gyr, or z~14. Correct?  The meaning
; of b,Tau, and SFR(100Myr) The EW sub-panel seems truncated by the label.
    
    common com_macs1149, allpost, allmodel

    isedfit_dir = getenv('CLASH_PROJECTS')+'/hff/macs1149/'
    montegrids_dir = isedfit_dir+'montegrids/'

    isedfit_paramfile = isedfit_dir+'macs1149_paramfile.par'
    suffix = 'jan28'

    filt = hff_filterlist(pivotwave=weff,width=hwhm,/useirac,/forwei)
    weff = weff/1D4
    hwhm = hwhm/1D4/2

    cat = rsex(isedfit_dir+'flx_iso.'+suffix)
    cat = struct_addtags(cat,replicate({logmu: 1.0},n_elements(cat)))
    
    this = where(cat.id eq 663)
    cat = cat[this]
    cat.logmu = alog10(10.0) ; hack!

    nobj = n_elements(cat)
    
; restore the models
    sfhgrid = [1]
    nsfhgrid = n_elements(sfhgrid)
    if n_elements(allpost) eq 0L then begin
       for ii = 0, nsfhgrid-1 do begin
          model1 = read_isedfit(isedfit_paramfile,index=this,/getmodels,$
            thissfhgrid=ii+1,isedfit_dir=isedfit_dir,isedfit_post=post1,$
            montegrids_dir=montegrids_dir,fnu=fnu)
          if ii eq 0 then allpost = post1 else allpost = [allpost,post1]
          if ii eq 0 then allmodel = model1 else allmodel = [allmodel,model1]
       endfor
    endif

    mstar_50 = allmodel.mstar_50-cat.logmu
    mstar_err = allmodel.mstar_err
    sfr_50 = 10^(allmodel.sfr_50-cat.logmu)
    sfr_err = allmodel.sfr_err*sfr_50*alog(10)
    sfrage_50 = allmodel.sfrage_50
    sfrage95 = fltarr(nobj)
    for ii = 0, nobj-1 do sfrage95[ii] = weighted_quantile(allpost[ii].sfrage*1E3,quant=95)

    zform = getredshift(getage(cat.bpz)-sfrage_50)
    
    splog, 'ID, mass, SFR, sSFR, age, age95, t(z), zform'
    niceprint, cat.id, mstar_50, sfr_50, 1D9*sfr_50/10.0^mstar_50, $
      sfrage_50*1E3, sfrage95*1E3, getage(cat.bpz)*1E3, zform

    splog, 'Mean, minmax logmu', djs_median(cat.logmu), minmax(cat.logmu)
    splog, 'Mean, minmax mass, mean err ', djs_median(10^mstar_50), $
      minmax(10^mstar_50), djs_median(mstar_err*10^mstar_50*alog(10))
    splog, 'Mean, minmax SFR, mean err ', djs_median(sfr_50), $
      minmax(sfr_50), djs_median(sfr_err)
    splog, 'Mean, minmax zform', djs_median(zform), minmax(zform)
    splog, 'Mean, minmax sfrage95', djs_median(sfrage95), minmax(sfrage95)
    splog, 'Doubling time ', 2*djs_median(10^mstar_50)/djs_median(sfr_50)/1D6
    splog, djs_median(allmodel.z), getage(djs_median(allmodel.z))

; --------------------------------------------------
; paper plot: posteriors on stellar mass, SFR, and SFR-weighted age
    psfile = isedfit_dir+'macs1149_id663_post.ps'
    im_plotconfig, 13, pos, psfile=psfile, xspace=0.2*[1,1], yspace=0.0, $
      xmargin=[0.8,0.4], height=2.0, width=2.3*[1,1,1], charsize=1.2, $
      color_outline='black'

    render_postplot, 10^(allpost.mstar-cat.logmu), pos[*,0], xrange=[5D7,4D9], $
      ytitle='', /logbins, binsize=0.15, /xlog, $
      xtitle='Stellar Mass [(\mu/10)^{-1} M'+sunsymbol()+']', $
      color_outline='black'
    im_legend, repstr(strtrim(string(10^mstar_50,format='(E12.2)'),2),'E+0','\times10^{')+$
      '} M'+sunsymbol(), /left, /top, box=0, margin=0, charsize=1.1
    djs_oplot, 10^mstar_50*[1,1], !y.crange, line=5, thick=6
    
    render_postplot, 10^(allpost.sfr-cat.logmu), pos[*,1], /noerase, xrange=[0.2,2.3], $
      xtitle='SFR [(\mu/10)^{-1} M'+sunsymbol()+' yr^{-1}]', bin=0.15, $
      color_outline='black' ;, monte=allmodel.bigsfr+10
    im_legend, strtrim(string(sfr_50,format='(F12.2)'),2)+' M'+sunsymbol()+' yr^{-1}', $
      /left, /top, box=0, margin=0, charsize=1.1
    djs_oplot, sfr_50*[1,1], !y.crange, line=5, thick=6

    xr = [0,350]
    agelabel = [0.0, 100, 200, 300] ; [Myr]

    render_postplot, allpost.sfrage*1E3, pos[*,2], xrange=xr, /noerase, $
      xtitle='<t>_{SFR} (Myr)', xtickinterval=100, color_outline='black', bin=40, $
      xstyle=9
    im_legend, strtrim(string(1E3*sfrage_50,format='(I0)'),2)+' Myr', $
      /left, /top, box=0, margin=0, charsize=1.1
    djs_oplot, sfrage_50*1E3*[1,1], !y.crange, line=5, thick=6, /save
    
     xtickv = getredshift(getage(cat.bpz)-agelabel/1E3)
     xticks = n_elements(xtickv)-1
     zxrange = getredshift(getage(cat.bpz)-!x.crange/1E3)
    
    axis, /xaxis, xsty=1, xtitle='Formation Redshift', xtickformat='age2zform', $
      xtickv=agelabel, xticks=xticks;, xrange=zxrange
;     xtickv=[11.1, 13.6, 18.1], xticks=2
;   axis, /xaxis, xsty=1, xrange=zxrange, xtitle='Formation Redshift', $
;     xtickv=xtickv, xticks=xticks

    xyouts, pos[0,0]-0.02, (pos[3,0]-pos[1,0])/2.0+pos[1,0], 'Probability', $
      align=0.5, orientation=90, /norm
    
    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

stop    

; --------------------------------------------------
; paper plot: SED plot adopting the BPZ redshift
    psfile = isedfit_dir+'macs1149_id663_sed.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=4.5

    col = 'firebrick'
    xrange = [0.3,7.0]
    ticks = loglevels(xrange)
    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=[30,22], $
      xsty=1, ysty=1, /xlog, position=pos, xtickv=ticks, $
      xticks=n_elements(ticks)-1, ytickinterval=2, $
      ytitle='AB Magnitude', xtitle='Observed-Frame Wavelength (\mu'+'m)'
    mask = (allmodel.wave gt 1214.5*(1+allmodel.z) and $
      allmodel.wave lt 1215*1.0025*(1+allmodel.z)) eq 1
    djs_oplot, allmodel.wave/1D4, djs_maskinterp(allmodel.flux,mask,allmodel.wave/1D4)
; model photometry
    notzero = where(allmodel.bestmaggies gt 0.0)
    bestmab = -2.5*alog10(allmodel.bestmaggies[notzero])
    djs_oplot, weff[notzero], bestmab, psym=symcat(6,thick=6), symsize=2.5, $
      color=cgcolor('dodger blue')
; overplot the data; distinguish between three different cases, based
; on the input photometry
    print, allmodel.maggies*sqrt(allmodel.ivarmaggies)
    nsigma = allmodel.maggies*0+2.0
;   nsigma[2] = 3.11 ; hack!
    print, nsigma
;   if ii eq 3 then nsigma[2] = 2.1
    mab = maggies2mag(allmodel.maggies,nsigma=nsigma,$
      ivar=allmodel.ivarmaggies,magerr=maberr,$
      lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
    used = where(mab gt -90.0,nused)
    upper = where(mab lt -90.0 and mabupper gt -90,nupper)
    if (nused ne 0L) then begin
       oploterror, weff[used], mab[used], hwhm[used], $
         mabhierr[used], psym=symcat(16), $
         symsize=1.7, color=im_color(col), /hibar, $
         errcolor=im_color(col), errthick=!p.thick
       oploterror, weff[used], mab[used], hwhm[used], $
         mabloerr[used], psym=3, color=im_color(col), /lobar, $
         errcolor=im_color(col), errthick=!p.thick
    endif
    if (nupper ne 0) then begin
       djs_oplot, weff[upper], mabupper[upper], $
         psym=symcat(11,thick=6), symsize=3.0, color=im_color('orange')
    endif
    im_legend, /left, /top, box=0, spacing=1.8, charsize=1.5, $
      ['M1149-JD','z = '+strtrim(string(cat.bpz,format='(F12.2)'),2)]

    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
return
end

