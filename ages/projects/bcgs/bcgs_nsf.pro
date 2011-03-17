pro nsf_makeplot, model, isedfit, filters=filters, psfile=psfile, $
  xrange=xrange1, yrange=yrange1, galaxy=galaxy
; support routine to make the actual plot
    
    filtinfo = im_filterspecs(filterlist=filters)
    nfilt = n_elements(filters)
    
    xtitle1 = textoidl('Observed Wavelength (\mu'+'m)')
    xtitle2 = textoidl('Rest Wavelength (\mu'+'m)')
    ytitle1 = textoidl('Apparent AB Magnitude')
;   ytitle1 = textoidl('m_{AB}')

    im_plotconfig, 8, pos, psfile=psfile, ymargin=[0.8,1.1]

; make the plot       
;   yrange = [28.5,18.5]
;   xrange1 = [3000.0,6E4]
    xrange2 = xrange1/(1.0+isedfit.zobj)
    plot, [0], [0], /nodata, xrange=xrange1/1D4, yrange=yrange1, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, /xlog, $
      xtickinterval=1000, position=pos
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2/1D4, /xlog
    djs_oplot, model.wave/1D4, model.flux, line=0, color='grey'

; overplot the filter names
    nice_filters = strarr(nfilt)
    for ii = 0, nfilt-1 do nice_filters[ii] = repstr(strmid(filters[ii],$
      strpos(filters[ii],'_',/reverse_search)+1),'.par','')
    nice_filters = repstr(nice_filters,'0','')
    nice_filters = repstr(repstr(nice_filters,'ch1','[3.6]'),'ch2','[4.5]')
    nice_filters = repstr(repstr(nice_filters,'ch3','[5.8]'),'ch4','[8.0]')
    for ii = 0, nfilt-1 do xyouts, filtinfo[ii].weff/1D4, $
      !y.crange[0]-0.3, nice_filters[ii], /data, align=0.5, charsize=1.6
    
; overplot the observed and model photometry
    used = where((isedfit.maggies gt 0.0) and $ ; used in the fitting
      (isedfit.ivarmaggies gt 0.0),nused)
    notused = where((isedfit.maggies gt 0.0) and $ ; not used in the fitting
      (isedfit.ivarmaggies eq 0.0),nnotused)
    nodata = where((isedfit.maggies eq 0.0) and $ ; no measurement
      (isedfit.ivarmaggies eq 0.0),nnodata)

    djs_oplot, filtinfo.weff/1D4, -2.5*alog10(isedfit.bestmaggies), $
      psym=symcat(6,thick=6), symsize=2.5
;   djs_oplot, filtinfo[used].weff, -2.5*alog10(isedfit.bestmaggies[used]), $
;     psym=symcat(6,thick=6), symsize=2.5

; overplot the various photometry
    ndwfs = where(strmatch(filters,'*ndwfs*',/fold))
    if (ndwfs[0] ne -1) then begin
       mab = maggies2mag(isedfit.maggies[ndwfs],$
         ivar=isedfit.ivarmaggies[ndwfs],magerr=mab_err)
       oploterror, filtinfo[ndwfs].weff/1D4, mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=fsc_color('dodger blue',101), $
         errcolor=fsc_color('dodger blue',101), errthick=!p.thick
    endif
    sdss = where(strmatch(filters,'*sdss*',/fold))
    if (sdss[0] ne -1) then begin
       mab = maggies2mag(isedfit.maggies[sdss],$
         ivar=isedfit.ivarmaggies[sdss],magerr=mab_err)
       oploterror, filtinfo[sdss].weff/1D4, mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=fsc_color('dodger blue',101), $
         errcolor=fsc_color('dodger blue',101), errthick=!p.thick
    endif
    newfirm = where(strmatch(filters,'*newfirm*',/fold))
    if (newfirm[0] ne -1) then begin
       mab = maggies2mag(isedfit.maggies[newfirm],$
         ivar=isedfit.ivarmaggies[newfirm],magerr=mab_err)
       oploterror, filtinfo[newfirm].weff/1D4, mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=fsc_color('black',101), $
         errcolor=fsc_color('black',101), errthick=!p.thick
    endif
    twomass = where(strmatch(filters,'*twomass*',/fold))
    if (twomass[0] ne -1) then begin
       mab = maggies2mag(isedfit.maggies[twomass],$
         ivar=isedfit.ivarmaggies[twomass],magerr=mab_err)
       oploterror, filtinfo[twomass].weff/1D4, mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=fsc_color('black',101), $
         errcolor=fsc_color('black',101), errthick=!p.thick
    endif
       
    irac = where(strmatch(filters,'*irac*',/fold))
    if (irac[0] ne -1) then begin
       mab = maggies2mag(isedfit.maggies[irac],$
         ivar=isedfit.ivarmaggies[irac],magerr=mab_err)
       oploterror, filtinfo[irac].weff/1D4, mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=fsc_color('firebrick',101), $
         errcolor=fsc_color('firebrick',101), errthick=!p.thick
    endif
       
    label = [$
      galaxy,$
      'z = '+string(isedfit.zobj,format='(F6.4)'),$
      'log (M/M'+sunsymbol()+') = '+strtrim(string(isedfit.mass_avg,format='(F12.2)'),2),$
;     '\pm'+strtrim(string(isedfit.mass_err,format='(F12.2)'),2),$
      'Age = '+strtrim(string(isedfit.age_avg,format='(F12.1)'),2)+' Gyr'];+$
;     '\pm'+strtrim(string(isedfit.age_err,format='(F12.2)'),2)+' Gyr']
     legend, textoidl(label), /left, /top, box=0, charsize=1.5, $
       margin=0, spacing=2.2

    im_plotconfig, psfile=psfile, /psclose

return
end

pro bcgs_nsf
; jm10oct12ucsd - pretty SED plots for the NSF proposal

    bcgspath = ages_path(/projects)+'bcgs/'

; lo-z example    
    index = 46
    lozparamfile = bcgspath+'bcgs_sdss_isedfit.par'
    lozparams = read_isedfit_paramfile(lozparamfile,sfhgrid=2)
    lozmodel = isedfit_restore(lozparamfile,lozisedfit,$
      params=lozparams,iopath=bcgspath,index=index)

    ss = mrdfits(bcgspath+'sdss_sample_v1.fits.gz',1)
    galaxy = 'SDSS J'+strmid(strcompress(im_dec2hms(ss[index].ra/15D),/rem),0,4)+$
      strmid(strcompress(im_dec2hms(ss[index].dec),/rem),0,5)
    
    nsf_makeplot, lozmodel, lozisedfit, $
      filters=strtrim(lozparams.filterlist,2), $
      psfile=bcgspath+'qaplot_nsf_loz.eps', $
      xrange=[3000,2.6E4], yrange=[19.5,12.0], galaxy=galaxy

; hi-z example
    index = 44
    hizparamfile = bcgspath+'bcgs_isedfit.par'
    hizparams = read_isedfit_paramfile(hizparamfile,sfhgrid=2)
    hizmodel = isedfit_restore(hizparamfile,hizisedfit,$
      params=hizparams,iopath=bcgspath,index=index)

    ss = rsex(bcgspath+'bcgs_sample_v3.sex')
    galaxy = 'NDWFS J'+strmid(strcompress(im_dec2hms(ss[index].ra/15D),/rem),0,4)+'+'+$
      strmid(strcompress(im_dec2hms(ss[index].dec),/rem),0,4)
    
    nsf_makeplot, hizmodel, hizisedfit, $
      filters=strtrim(hizparams.filterlist,2), $
      psfile=bcgspath+'qaplot_nsf_hiz.eps', $
      xrange=[3000,1E5], yrange=[25,18.5], galaxy=galaxy

return
end
    
