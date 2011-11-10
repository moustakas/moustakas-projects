pro macs0329_z6arcs, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber
; jm11nov08ucsd - build plots for the paper

    common com_post, mstar, post, age, Z, tau, sfr0, b100, av
    
    isedpath = clash_path(/ised)
    datapath = clash_path(/macs0329_z6arcs)

    isedfit_sfhgrid_dir = clash_path(/monte)
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/clash_sfhgrid.par'

    filters = strtrim(clash_filterlist(nice=nice_filters),2)
    filtinfo = im_filterspecs(filterlist=filters)

    ndraw = isedfit_ndraw()

    prefix = 'macs0329_z6arcs'
    paramfile = isedpath+prefix+'_supergrid02_isedfit.par'
    
; gather the photometry
    alladi = rsex(datapath+prefix+'.cat')
    allcat = read_clash_catalog('macs0329',/arcs)
    spherematch, allcat.ra, allcat.dec, 15D*hms2dec(alladi.ra), hms2dec(alladi.dec), 1D/3600, m1, m2
    cat = allcat[m1]
    adi = alladi[m2]
    struct_print, adi

; ---------------------------------------------------------------------------
; make a QAplot showing the best-fitting SEDs and posterior distributions
    galaxy = repstr(adi.id,'arc_','Arc ')
    model = isedfit_restore(paramfile,isedfit,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)

    isedfitpsfile = datapath+prefix+'_isedfit.ps'
    im_plotconfig, 0, pos, psfile=isedfitpsfile, ymargin=[1.0,1.1], $
      height=5.0

    xtitle1 = textoidl('Observed Wavelength (\AA)')
    xtitle2 = textoidl('Rest Wavelength (\AA)')
    ytitle1 = textoidl('m_{AB}')

    yrange = [30.5,24.0]
    xrange1 = [2000,17000]
    for ii = 0, n_elements(cat)-1 do begin
       xrange2 = xrange1/(1.0+isedfit[ii].zobj)
       plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
         xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=0, $
         xtickinterval=3000, position=pos, xtickformat='(I0)'
       axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=0

       oplot, model[ii].wave, model[ii].flux, line=0
       
; overplot the observed and model photometry
       used = where((isedfit[ii].maggies gt 0.0) and $ ; used in the fitting
         (isedfit[ii].ivarmaggies gt 0.0),nused)
       notused = where((isedfit[ii].maggies gt 0.0) and $ ; not used in the fitting
         (isedfit[ii].ivarmaggies eq 0.0),nnotused)
       nodata = where((isedfit[ii].maggies eq 0.0) and $ ; no measurement
         (isedfit[ii].ivarmaggies eq 0.0),nnodata)

       djs_oplot, filtinfo[used].weff, -2.5*alog10(isedfit[ii].bestmaggies[used]), $
         psym=symcat(6,thick=6), symsize=2.5

; overplot uvis, acs, and ir separately 
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
         margin=0, spacing=2.2

       label = [$
         'log (M/M'+sunsymbol()+') = '+strtrim(string(isedfit[ii].mass,format='(F12.2)'),2),$
         'Age = '+strtrim(string(isedfit[ii].age*1E3,format='(F12.1)'),2)+' Myr',$
         'Z/Z'+sunsymbol()+' = '+strtrim(string(isedfit[ii].Z/0.02,format='(F12.2)'),2),$
         'A_{V} = '+strtrim(string(isedfit[ii].av,format='(F12.3)'),2)+' mag',$
         '\tau = '+strtrim(string(isedfit[ii].tau,format='(F12.1)'),2)+' Gyr',$
         'log \psi_{100} = '+strtrim(string(isedfit[ii].sfr100,format='(F12.3)'),2)+' M'+sunsymbol()+' yr^{-1}',$
         'log (b_{100}) = '+strtrim(string(isedfit[ii].b100,format='(F12.3)'),2)]
       legend, textoidl(label), /left, /top, box=0, spacing=1.9, charsize=1.4, $
         position=[10000.0,yrange[0]-(yrange[0]-yrange[1])*0.44], /data
    endfor
    im_plotconfig, psfile=isedfitpsfile, /psclose

; plot the NxN distributions of parameters
    if (n_elements(mstar) eq 0L) then begin
       mstar = isedfit_reconstruct_posterior(paramfile,post=post,$
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath,$
         age=age,Z=Z,tau=tau,sfr0=sfr0,b100=b100,av=av)
    endif

; upper triangle
    for ii = 0, n_elements(cat)-1 do begin
       psfile = datapath+strlowcase(adi[ii].id)+'_manyd.ps'
       splog, 'Writing '+psfile
       manyd = transpose([[mstar[*,ii]],[age[*,ii]],[Z[*,ii]/0.02],[av[*,ii]],[sfr0[*,ii]]])
       label = textoidl(['log (M/M'+sunsymbol()+')','Age (Gyr)','Z/Z'+sunsymbol(),'A_{V} (mag)','log (\psi)'])

       im_manyd_scatterplot, fltarr(ndraw)+1, manyd, psfile, label=label, $
         axis_char_scale=1.5, /internal, outliers=0, $
         /nogrey, keynote=keynote, levels=errorf((dindgen(2)+1)/sqrt(2)), $
         /upper
    endfor       

; merge everything into a single postscript file and then convert to
; PDF
    outfile = datapath+prefix+'.pdf'
    allfiles = [isedfitpsfile,datapath+strlowcase(adi.id)+'_manyd.ps']
    spawn, 'gs -q -dNOPAUSE -sDEVICE=pdfwrite '+$
      '-sOutputFile='+outfile+' -dBATCH '+strjoin(allfiles,' '), /sh
    
stop    
    
return
end
