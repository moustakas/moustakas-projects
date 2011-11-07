pro talk_11oct_clash, keynote=keynote
; jm11oct14ucsd - miscellaneous plots for the 2011 Oct CLASH team
; meeting 

    datapath = getenv('IM_RESEARCH_DIR')+'/talks/2011/11oct_clash/'
    if keyword_set(keynote) then talkpath = datapath+'keynote/' else $
      talkpath = datapath

    isedpath = clash_path(/isedfit)
    isedfit_sfhgrid_dir = clash_path(/montegrids)

    filters = strtrim(clash_filterlist(nice=nice_filters),2)
    filtinfo = im_filterspecs(filterlist=filters)

    ndraw = isedfit_ndraw()

; --------------------------------------------------
; BCG SED-fitting example
    paramfile = isedpath+'bcgs_supergrid10_isedfit.par'
    model = isedfit_restore(paramfile,isedfit,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)
    
    mstar = isedfit_reconstruct_posterior(paramfile,post=post,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath,$
      age=age,Z=Z,tau=tau,sfr100=sfr100,b100=b100,av=av)
    
    xtitle1 = textoidl('Observed Wavelength (\AA)')
    xtitle2 = textoidl('Rest Wavelength (\AA)')
    ytitle1 = textoidl('m_{AB}')

    yrange = [25.5,13.5]
    xrange1 = [2500,2.0E4]
    xrange2 = xrange1/(1.0+isedfit[0].zobj)
    
    psfile = talkpath+'abell2261_bcg_isedfit.ps'
    im_plotconfig, 8, pos, psfile=psfile, ymargin=[1.0,1.1], keynote=keynote

    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'

    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, /xlog, $
      xtickinterval=1000, position=pos, color=im_color(keycolor)
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, /xlog, color=im_color(keycolor)

    col = ['goldenrod','dark green','midnight blue']
    for ii = 0, 2 do oplot, model[ii].wave, $
      model[ii].flux, line=0, color=im_color(col[ii])

    for ii = 0, 2 do begin
; overplot the observed and model photometry
       used = where((isedfit[ii].maggies gt 0.0) and $ ; used in the fitting
         (isedfit[ii].ivarmaggies gt 0.0),nused)
       notused = where((isedfit[ii].maggies gt 0.0) and $ ; not used in the fitting
         (isedfit[ii].ivarmaggies eq 0.0),nnotused)
       nodata = where((isedfit[ii].maggies eq 0.0) and $ ; no measurement
         (isedfit[ii].ivarmaggies eq 0.0),nnodata)

       djs_oplot, filtinfo[used].weff, -2.5*alog10(isedfit[ii].bestmaggies[used]), $
         psym=symcat(6,thick=6), symsize=2.5

; overplot galex, cfhtls, and irac
       uvis = where(strmatch(nice_filters,'*uvis*',/fold))
       mab = maggies2mag(isedfit[ii].maggies[uvis],$
         ivar=isedfit[ii].ivarmaggies[uvis],magerr=mab_err)
       oploterror, filtinfo[uvis].weff, mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=fsc_color('dodger blue',101), $
         errcolor=fsc_color('dodger blue',101), errthick=!p.thick
       
       acs = where(strmatch(nice_filters,'*acs*',/fold))
       mab = maggies2mag(isedfit[ii].maggies[acs],$
         ivar=isedfit[ii].ivarmaggies[acs],magerr=mab_err)
       oploterror, filtinfo[acs].weff, mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=fsc_color('tan',101), $
         errcolor=fsc_color('tan',101), errthick=!p.thick
       
       ir = where(strmatch(nice_filters,'*ir*',/fold))
       mab = maggies2mag(isedfit[ii].maggies[ir],$
         ivar=isedfit[ii].ivarmaggies[ir],magerr=mab_err)
       oploterror, filtinfo[ir].weff, mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=fsc_color('firebrick',101), $
         errcolor=fsc_color('firebrick',101), errthick=!p.thick
    endfor
       
    label = ['Abell 2261 BCG','z = '+string(isedfit[0].zobj,format='(F6.4)')]
    legend, textoidl(label), /left, /top, box=0, charsize=1.7, $
      margin=0, textcolor=im_color(keycolor), spacing=2.2

; inset with P(M)
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.05], $
      xrange=[11.0,12.5], position=[0.52,0.25,0.9,0.53], ytitle='P(M)', $
      xtitle='Stellar Mass', xtickinterval=0.5, ytickname=replicate(' ',10), $
      charsize=1.5, color=im_color(keycolor)
    im_plothist, mstar[*,0], bin=0.03, /peak, /overplot, /fill, fcolor=im_color('goldenrod')
    im_plothist, mstar[*,1], bin=0.03, /peak, /overplot, /fill, fcolor=im_color('dark green')
    im_plothist, mstar[*,2], bin=0.03, /peak, /overplot, /fill, fcolor=im_color('midnight blue')

; maximum likelihood values
    oplot, isedfit[0].mass*[1,1], !y.crange, line=5, thick=5, color=im_color(keycolor)
    oplot, isedfit[1].mass*[1,1], !y.crange, line=5, thick=5, color=im_color(keycolor)
    oplot, isedfit[2].mass*[1,1], !y.crange, line=5, thick=5, color=im_color(keycolor)
    
    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

; plot the NxN distributions of parameters for the BCG and the knot
    psfile = talkpath+'abell2261_bcg_manyd.ps'
    indx = 2 ; =total flux
    manyd = transpose([[mstar[*,indx]],[age[*,indx]],[Z[*,indx]/0.02],[b100[*,indx]]])
    label = textoidl(['log (M/M'+sunsymbol()+')','Age (Gyr)','Z/Z'+sunsymbol(),'log (b_{100})'])
;   manyd = transpose([[mstar[*,indx]],[age[*,indx]],[Z[*,indx]/0.02],[tau[*,indx]],[b100[*,indx]]])
;   label = textoidl(['log (M/M'+sunsymbol()+')','Age (Gyr)','Z/Z'+sunsymbol(),$
;     '\tau (Gyr)','b_{100} (\times10^{-3})'])
;   range = [[12.0,12.5],[7.5,11.5],[0.75,1.75],[-0.1,1.3],[-0.5,5.0]]

    im_manyd_scatterplot, fltarr(ndraw)+1, manyd, psfile, label=label, $
      axis_char_scale=1.5, /internal, nsig=4D, /outliers, $
      /nogrey, keynote=keynote, levels=errorf((dindgen(2)+1)/sqrt(2)), $
      /upper, in_nticks=3

; plot the SED-fit for the knot and for the NxN distribution of
; parameters
    psfile = talkpath+'abell2261_knot_isedfit.ps'
    im_plotconfig, 8, pos, psfile=psfile, ymargin=[1.0,1.1], keynote=keynote

    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'

    yrange = [26.5,19.0]
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, /xlog, $
      xtickinterval=1000, position=pos, color=im_color(keycolor)
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, /xlog, color=im_color(keycolor)

;   oplot, model[2].wave, model[2].flux, line=1
    oplot, model[3].wave, model[3].flux, line=0, color=im_color(keycolor)

; overplot the observed and model photometry
    used = where((isedfit[3].maggies gt 0.0) and $ ; used in the fitting
      (isedfit[3].ivarmaggies gt 0.0),nused)
    notused = where((isedfit[3].maggies gt 0.0) and $ ; not used in the fitting
      (isedfit[3].ivarmaggies eq 0.0),nnotused)
    nodata = where((isedfit[3].maggies eq 0.0) and $ ; no measurement
      (isedfit[3].ivarmaggies eq 0.0),nnodata)

    djs_oplot, filtinfo[used].weff, -2.5*alog10(isedfit[3].bestmaggies[used]), $
      psym=symcat(6,thick=6), symsize=2.5

; overplot galex, cfhtls, and irac
    uvis = where(strmatch(nice_filters,'*uvis*',/fold))
    mab = maggies2mag(isedfit[3].maggies[uvis],$
      ivar=isedfit[3].ivarmaggies[uvis],magerr=mab_err)
    oploterror, filtinfo[uvis].weff, mab, mab_err, psym=symcat(16), $
      symsize=2.0, color=fsc_color('dodger blue',101), $
      errcolor=fsc_color('dodger blue',101), errthick=!p.thick
    
    acs = where(strmatch(nice_filters,'*acs*',/fold))
    mab = maggies2mag(isedfit[3].maggies[acs],$
      ivar=isedfit[3].ivarmaggies[acs],magerr=mab_err)
    oploterror, filtinfo[acs].weff, mab, mab_err, psym=symcat(16), $
      symsize=2.0, color=fsc_color('tan',101), $
      errcolor=fsc_color('tan',101), errthick=!p.thick
    
    ir = where(strmatch(nice_filters,'*ir*',/fold))
    mab = maggies2mag(isedfit[3].maggies[ir],$
      ivar=isedfit[3].ivarmaggies[ir],magerr=mab_err)
    oploterror, filtinfo[ir].weff, mab, mab_err, psym=symcat(16), $
      symsize=2.0, color=fsc_color('firebrick',101), $
      errcolor=fsc_color('firebrick',101), errthick=!p.thick

    label = ['Abell 2261 BCG knot','z = '+string(isedfit[3].zobj,format='(F6.4)')]
    legend, textoidl(label), /left, /top, box=0, charsize=1.7, $
      margin=0, textcolor=im_color(keycolor), spacing=2.2

; inset with P(M)
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.05], $
      xrange=[9.7,10.15], position=[0.52,0.25,0.9,0.53], ytitle='P(M)', $
      xtitle='Stellar Mass', xtickinterval=0.1, ytickname=replicate(' ',10), $
      charsize=1.5, color=im_color(keycolor)
    im_plothist, mstar[*,3], bin=0.03, /peak, /overplot, /fill, fcolor=im_color('goldenrod')
    oplot, isedfit[3].mass*[1,1], !y.crange, line=5, thick=5, color=im_color(keycolor)
    
    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

    psfile = talkpath+'abell2261_knot_manyd.ps'
    indx = 3 ; =total flux
    manyd = transpose([[mstar[*,indx]],[age[*,indx]],[Z[*,indx]/0.02],[b100[*,indx]]])
    label = textoidl(['log (M/M'+sunsymbol()+')','Age (Gyr)','Z/Z'+sunsymbol(),'log (b_{100})'])

    im_manyd_scatterplot, fltarr(ndraw)+1, manyd, psfile, label=label, $
      axis_char_scale=1.5, /internal, nsig=4D, /outliers, $
      /nogrey, keynote=keynote, levels=errorf((dindgen(2)+1)/sqrt(2)), $
      /upper
    
stop
    
; --------------------------------------------------
; arc SED-fitting examples
    galaxy = ['MACS1206-01274','MACS1206-01794']
    index = [9,10]
    
    paramfile = isedpath+'arcs_supergrid01_isedfit.par'
    model = isedfit_restore(paramfile,isedfit,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,index=index)
    
    xtitle1 = textoidl('Observed Wavelength (\AA)')
    xtitle2 = textoidl('Rest Wavelength (\AA)')
    ytitle1 = textoidl('m_{AB}')

    yrangeall = [[24.0,19.0],[28.5,21.0]]
    xrangeall = [[2500,2.0E4],[2500,2.0E4]]
    for ii = 0, 1 do begin
       yrange = yrangeall[*,ii]
       xrange1 = xrangeall[*,ii]
       xrange2 = xrange1/(1.0+isedfit[ii].zobj)

       psfile = talkpath+strlowcase(galaxy[ii])+'_isedfit.ps'
       im_plotconfig, 8, pos, psfile=psfile, ymargin=[1.0,1.1], keynote=keynote

       if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'

       plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
         xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, /xlog, $
         xtickinterval=1000, position=pos, color=im_color(keycolor)
       axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, /xlog, $
         color=im_color(keycolor)

       oplot, model[ii].wave, model[ii].flux, line=0, color=im_color(keycolor)
       
; overplot the observed and model photometry
       used = where((isedfit[ii].maggies gt 0.0) and $ ; used in the fitting
         (isedfit[ii].ivarmaggies gt 0.0),nused)
       notused = where((isedfit[ii].maggies gt 0.0) and $ ; not used in the fitting
         (isedfit[ii].ivarmaggies eq 0.0),nnotused)
       nodata = where((isedfit[ii].maggies eq 0.0) and $ ; no measurement
         (isedfit[ii].ivarmaggies eq 0.0),nnodata)

       djs_oplot, filtinfo[used].weff, -2.5*alog10(isedfit[ii].bestmaggies[used]), $
         psym=symcat(6,thick=6), symsize=2.5

; overplot galex, cfhtls, and irac
       uvis = where(strmatch(nice_filters,'*uvis*',/fold))
       mab = maggies2mag(isedfit[ii].maggies[uvis],$
         ivar=isedfit[ii].ivarmaggies[uvis],magerr=mab_err)
       oploterror, filtinfo[uvis].weff, mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=fsc_color('dodger blue',101), $
         errcolor=fsc_color('dodger blue',101), errthick=!p.thick
       
       acs = where(strmatch(nice_filters,'*acs*',/fold))
       mab = maggies2mag(isedfit[ii].maggies[acs],$
         ivar=isedfit[ii].ivarmaggies[acs],magerr=mab_err)
       oploterror, filtinfo[acs].weff, mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=fsc_color('tan',101), $
         errcolor=fsc_color('tan',101), errthick=!p.thick
       
       ir = where(strmatch(nice_filters,'*ir*',/fold))
       mab = maggies2mag(isedfit[ii].maggies[ir],$
         ivar=isedfit[ii].ivarmaggies[ir],magerr=mab_err)
       oploterror, filtinfo[ir].weff, mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=fsc_color('firebrick',101), $
         errcolor=fsc_color('firebrick',101), errthick=!p.thick
    
       label = [galaxy[ii],'z = '+string(isedfit[ii].zobj,format='(F6.4)')]
       legend, textoidl(label), /left, /top, box=0, charsize=1.7, $
         margin=0, textcolor=im_color(keycolor), spacing=2.2

       label = [$
         'log (M/M'+sunsymbol()+') = '+strtrim(string(isedfit[ii].mass,format='(F12.2)'),2),$
         'Age = '+strtrim(string(isedfit[ii].age,format='(F12.2)'),2)+' Gyr',$
         'Z/Z'+sunsymbol()+' = '+strtrim(string(isedfit[ii].Z/0.02,format='(F12.2)'),2),$
         'A_{V} = '+strtrim(string(isedfit[ii].av*isedfit[ii].mu,format='(F12.3)'),2)+' mag',$
;        '\tau = '+strtrim(string(isedfit[ii].tau,format='(F12.1)'),2)+' Gyr',$
;        'log \psi_{100} = '+strtrim(string(isedfit[ii].sfr100,format='(F12.3)'),2)+' M'+sunsymbol()+' yr^{-1}']
         'log (b_{100}) = '+strtrim(string(isedfit[ii].b100,format='(F12.3)'),2)]
       legend, textoidl(label), /left, /top, box=0, spacing=1.9, charsize=1.5, $
         position=[10000.0,yrange[0]-(yrange[0]-yrange[1])*0.3], /data, $
         textcolor=im_color(keycolor)
       im_plotconfig, psfile=psfile, /psclose, keynote=keynote
    endfor

; plot the NxN distributions of parameters
    common com_arc, arc_mstar, arc_post, arc_age, arc_Z, arc_tau, arc_sfr100, arc_b100, arc_av
    if (n_elements(arc_mstar) eq 0L) then begin
       arc_mstar = isedfit_reconstruct_posterior(paramfile,index=index,post=arc_post,$
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath,$
         age=arc_age,Z=arc_Z,tau=arc_tau,sfr100=arc_sfr100,b100=arc_b100,av=arc_av)
    endif

; upper triangle
    for ii = 0, 1 do begin
       psfile = talkpath+strlowcase(galaxy[ii])+'_manyd.ps'
       splog, 'Writing '+psfile
       manyd = transpose([[arc_mstar[*,ii]],[arc_age[*,ii]],[arc_Z[*,ii]/0.02],[arc_av[*,ii]],[arc_b100[*,ii]]])
       label = textoidl(['log (M/M'+sunsymbol()+')','Age (Gyr)','Z/Z'+sunsymbol(),'A_{V} (mag)','log (b_{100})'])
;      manyd = transpose([[arc_mstar[*,ii]],[arc_Z[*,ii]/0.02],[arc_tau[*,ii]],[arc_av[*,ii]],[arc_sfr100[*,ii]]])
;      label = textoidl(['log (M/M'+sunsymbol()+')','Z/Z'+sunsymbol(),$
;        '\tau (Gyr)','A_{V} (mag)','\psi_{100} (M'+sunsymbol()+' yr^{-1})'])

       im_manyd_scatterplot, fltarr(ndraw)+1, manyd, psfile, label=label, $
         axis_char_scale=1.5, /internal, nsig=4D, outliers=0, $
         /nogrey, keynote=keynote, levels=errorf((dindgen(2)+1)/sqrt(2)), $
         /upper
    endfor       

return
end
