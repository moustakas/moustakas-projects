function clash_restore, supergrid, index=index, isedfit=isedfit

    isedpath = clash_path(/ised)
    sfhgrid_basedir = clash_path(/monte)
    prefix = 'a383'

    paramfile = isedpath+prefix+'_supergrid'+string(supergrid,$
      format='(i2.2)')+'_isedfit.par'
    model = isedfit_restore(paramfile,isedfit,params=params,$
      iopath=isedpath,index=index,sfhgrid_basedir=sfhgrid_basedir)
    
return, model
end

pro qaplot_clash_isedfit
; jm11may05ucsd - build some QAplots comparing the various supergrids 

    clashpath = clash_path(/ised)
    sample = mrdfits(clash_path(/cat)+'a383.fits.gz',1)
    sample = [sample,sample] ; HACK!
    galaxy = sample.galaxy

; restore the filters and best-fitting models
    filterlist = clash_filterlist()
    filtinfo = im_filterspecs(filterlist=filterlist)
    weff = filtinfo.weff
    fwhm = filtinfo.fwhm
    hwhm = filtinfo.fwhm/2.0
    nfilt = n_elements(filterlist)

    model_bc03 = clash_restore(1,isedfit=bc03)
    model_fsps = clash_restore(2,isedfit=fsps)
    model_mara = clash_restore(3,isedfit=mara)
    ngal = n_elements(bc03)
    
; generate the QA-plot
    xtitle1 = 'Observed Wavelength (\AA)'
    xtitle2 = 'Rest Wavelength (\AA)'
    ytitle1 = 'm_{AB}'

    bc03_color = 'black'
    ss_color = 'dodger blue'
    ae_color = 'firebrick'
    
    psfile = clashpath+'qaplot_a383_isedfit.ps'
    im_plotconfig, 8, pos, psfile=psfile, ymargin=[1.0,1.1]
    for igal = 0, ngal-1 do begin
       if ((igal mod 10) eq 0) then print, igal+1L, ngal, string(13b), $
         format='("Building QAplot for galaxy ",I0,"/",I0,A10,$)'

; redshift and best-fit model
       zobj = bc03[igal].zobj

       xrange1 = [min(weff-1.3*fwhm),max(weff+2*fwhm)]
       xrange1[0] = xrange1[0]>(300.0*(1+zobj))
       xrange1[1] = xrange1[1]<max(model_bc03[igal].wave)

       get_element, model_bc03[igal].wave, xrange1, xx
       xrange2 = xrange1/(1.0+zobj)
       yrange = [max(model_bc03[igal].flux[xx[0]:xx[1]]),min(model_bc03[igal].flux[xx[0]:xx[1]])]
       pad = (yrange[0]-yrange[1])*0.1
       yrange = yrange+pad*[1,-1]
       
       djs_plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
         xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, /xlog, $
         xtickinterval=1000, position=pos ;ymargin=[4,3], 
       axis, /xaxis, xsty=1, xtitle=textoidl(xtitle2), xrange=xrange2, /xlog

       label = textoidl([strtrim(repstr(galaxy[igal],'_',' '),2),$
         'z = '+string(zobj,format='(F6.4)')])
       legend, label, /left, /top, box=0, spacing=1.5, charsize=1.4;, margin=0

; overplot the best-fitting models       
       djs_oplot, model_bc03[igal].wave, model_bc03[igal].flux, line=0, color=fsc_color(bc03_color,101)
       djs_oplot, model_fsps[igal].wave, model_fsps[igal].flux, line=0, color=fsc_color(ss_color,101)
       djs_oplot, model_mara[igal].wave, model_mara[igal].flux, line=0, color=fsc_color(ae_color,101)

;      im_legend, ['BC03','BaSTI-ss','BaSTI-ae'], /right, /top, box=0, $
;        charsize=1.4, color=[bc03_color,ss_color,ae_color], margin=0, $
;        line=0, thick=4, pspacing=1.6
       
;      bestmab = -2.5*alog10(isedfit[igal].bestmaggies)
;      djs_oplot, weff, bestmab, psym=symcat(6,thick=6), symsize=2.5

; overplot the photometry
       used = where((bc03[igal].maggies gt 0.0) and $ ; used in the fitting
         (bc03[igal].ivarmaggies gt 0.0),nused)
       notused = where((bc03[igal].maggies gt 0.0) and $ ; not used in the fitting
         (bc03[igal].ivarmaggies eq 0.0),nnotused)
       nodata = where((bc03[igal].maggies eq 0.0) and $ ; no measurement
         (bc03[igal].ivarmaggies eq 0.0),nnodata)
       upper = where((bc03[igal].maggies le 0.0) and $ ; upper limit
         (bc03[igal].ivarmaggies gt 0.0),nupper)

       if (nused ne 0L) then begin
          mab = maggies2mag(bc03[igal].maggies[used],$
            ivar=bc03[igal].ivarmaggies[used],magerr=mab_err)
          oploterror, weff[used], mab, hwhm[used], mab_err, psym=symcat(6,thick=6), $
            symsize=3.0, color=djs_icolor('black'), $
            errcolor=djs_icolor('black'), errthick=!p.thick
       endif

       if (nnotused ne 0L) then begin
          mab = maggies2mag(bc03[igal].maggies[notused])
          oploterror, weff[notused], mab, hwhm[notused], mab*0.0, $
            psym=symcat(4,thick=6.0), symsize=3.0, color=djs_icolor('tan'), $
            errcolor=djs_icolor('tan'), errthick=!p.thick
       endif

       if (nupper ne 0) then begin
          mab = maggies2mag(1.0/sqrt(bc03[igal].ivarmaggies[upper]))
          oploterror, weff[upper], mab, hwhm[upper], mab*0.0, psym=symcat(18), $
            symsize=3.0, color=djs_icolor('blue'), $
            errcolor=djs_icolor('blue'), errthick=!p.thick
       endif

; legend
       label = [$
         'BC03: log(M/M'+sunsymbol()+')='+strtrim(string(bc03[igal].mass,format='(F12.2)'),2)+$
         ', Z/Z_{'+sunsymbol()+'}='+strtrim(string(bc03[igal].Z/0.02,format='(F12.3)'),2)+$
;        ', Age='+strtrim(string(bc03[igal].age,format='(F12.2)'),2)+' Gyr'+$
         ', b_{100}='+strtrim(string(bc03[igal].b100,format='(F12.3)'),2)+$
;        ', \psi_{100}='+strtrim(string(bc03[igal].sfr100,format='(F12.2)'),2)+' M_{'+sunsymbol()+'} yr^{-1}'+$
         ', \chi^2='+strtrim(string(bc03[igal].chi2,format='(F12.2)'),2),$
         'FSPS: log(M/M'+sunsymbol()+')='+strtrim(string(fsps[igal].mass,format='(F12.2)'),2)+$
         ', Z/Z_{'+sunsymbol()+'}='+strtrim(string(fsps[igal].Z/0.02,format='(F12.3)'),2)+$
;        ', Age='+strtrim(string(fsps[igal].age,format='(F12.2)'),2)+' Gyr'+$
         ', b_{100}='+strtrim(string(fsps[igal].b100,format='(F12.3)'),2)+$
;        ', \psi_{100}='+strtrim(string(fsps[igal].sfr100,format='(F12.2)'),2)+' M_{'+sunsymbol()+'} yr^{-1}'+$
         ', \chi^2='+strtrim(string(fsps[igal].chi2,format='(F12.2)'),2),$
         'Maraston+05: log(M/M'+sunsymbol()+')='+strtrim(string(mara[igal].mass,format='(F12.2)'),2)+$
         ', Z/Z_{'+sunsymbol()+'}='+strtrim(string(mara[igal].Z/0.02,format='(F12.3)'),2)+$
;        ', Age='+strtrim(string(mara[igal].age,format='(F12.2)'),2)+' Gyr'+$
         ', b_{100}='+strtrim(string(mara[igal].b100,format='(F12.3)'),2)+$
;        ', \psi_{100}='+strtrim(string(mara[igal].sfr100,format='(F12.2)'),2)+' M_{'+sunsymbol()+'} yr^{-1}'+$
         ', \chi^2='+strtrim(string(mara[igal].chi2,format='(F12.2)'),2)]
       im_legend, textoidl(label), /right, /bottom, box=0, spacing=1.7, charsize=1.2, margin=0, $
         textcolor=[bc03_color,ss_color,ae_color]

    endfor 
    im_plotconfig, psfile=psfile, /psclose, /gzip

stop    

return
end
