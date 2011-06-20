function photoztemplates_restore, info, vname=vname
; internal function to rebuild the best-fitting SED in the observed 
; frame  

    ngal = n_elements(info)
    light = 2.99792458D18       ; speed of light [A/s]

    k_load_vmatrix, vmatrix, lambda, vfile=vfile, $
      lfile=lfile, vpath=vpath, vname=vname
    wave = k_lambda_to_centers(lambda)
    model = {wave: wave, flux: wave*0.0}
    model = replicate(temporary(model),ngal)
    model.flux = vmatrix#info.k_coeffs
    for igal = 0L, ngal-1L do begin
       model[igal].wave = wave*(1.0+info[igal].k_zobj) ; [A]
       model[igal].flux = model[igal].flux/(1.0+info[igal].k_zobj)
    endfor

    model.flux = model.flux*model.wave^2/light ; observed frame
    model.flux = -2.5*alog10(model.flux>1D-50)-48.6 ; [AB mag]
    
return, model
end

pro photoztemplates_qaplot, emlines=emlines

    path = sings_path(/proj)+'photoztemplates/'
    if keyword_set(emlines) then begin
       vname = 'default'
       suffix = '_emlines'
    endif else begin
       vname = 'default.nolines'
       suffix = ''
    endelse
    info = mrdfits(path+'kcorr'+suffix+'.fits.gz',1)
    spec = mrdfits(path+'optspec.fits.gz',1)
    
    ngal = n_elements(info)
    if (ngal eq 0L) then begin
       doc_library, 'photoztemplates_qaplot'
       return
    endif

; restore the filters and best-fitting models
    in_filterlist = photoztemplates_filterlist()
    in_filtinfo = im_filterspecs(filterlist=in_filterlist)
    in_nfilt = n_elements(in_filterlist)
    model = photoztemplates_restore(info,vname=vname)

; generate the QA-plot
    xtitle1 = 'Observed Wavelength (\AA)'
    xtitle2 = 'Rest Wavelength (\AA)'
    ytitle1 = 'm_{AB}'

    psfile = 'qaplot_photoztemplates'+suffix+'.ps'
    im_plotconfig, 8, pos, height=5.5, psfile=psfile, $
      ymargin=[0.9,1.1]

    for igal = 0L, ngal-1L do begin
       if ((igal mod 10) eq 0) then print, igal+1L, ngal, string(13b), $
         format='("Building QAplot for galaxy ",I0,"/",I0,A10,$)'
       
; redshift and best-fit model
       zobj = info[igal].k_zobj
       wave = model[igal].wave    ; [A]
       flux = model[igal].flux    ; [AB]

; get the x-range
;      xrange1 = [min(in_filtinfo.weff-1.3*in_filtinfo.fwhm),$
;        max(in_filtinfo.weff+2.1*in_filtinfo.fwhm)]
;      xrange1[0] = xrange1[0]>(500.0*(1+zobj))
;      xrange1[1] = (xrange1[1]<max(wave));<4E5
       xrange1 = [850,2E6]
       get_element, wave, xrange1, xx

; input photometry: distinguish between three different cases
       used = where((info[igal].k_maggies gt 0.0) and $ ; used in the fitting
         (info[igal].k_ivarmaggies gt 0.0),nused)
       notused = where((info[igal].k_maggies gt 0.0) and $ ; not used in the fitting
         (info[igal].k_ivarmaggies eq 0.0),nnotused)
       nodata = where((info[igal].k_maggies eq 0.0) and $ ; no measurement
         (info[igal].k_ivarmaggies eq 0.0),nnodata)

       if (nused ne 0) then mab_used = maggies2mag(info[igal].k_maggies[used],$
         ivar=info[igal].k_ivarmaggies[used],magerr=mab_used_err)
       if (nnotused ne 0L) then mab_notused = maggies2mag(info[igal].k_maggies[notused])
       
; best-fit photometry       
       bestmab = -2.5*alog10(info[igal].k_bestmaggies)
;      yrange = [max(bestmab)*1.05,min(bestmab)*0.95]
       yrange = [(max(mab_used)>max(mab_notused))*1.1,(min(mab_used)<min(mab_notused))*0.9]

; make the plot          
       djs_plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
         xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, /xlog, $
         position=pos
       axis, /xaxis, xsty=1, xtitle=textoidl(xtitle2), xrange=xrange2, /xlog
       djs_oplot, wave, flux, line=0, color='grey'
       isdata = where(spec.isdata[*,igal],nisdata)
       if (nisdata ne 0L) then begin
          notdata_blue = where(spec.isdata[*,igal] eq 0 and spec.wave lt min(spec.wave[isdata]))
          notdata_red = where(spec.isdata[*,igal] eq 0 and spec.wave gt max(spec.wave[isdata]))
          djs_oplot, spec.wave[isdata], spec.flux[isdata,igal], color='orange', thick=4
          djs_oplot, spec.wave[notdata_blue], spec.flux[notdata_blue,igal], color='blue', thick=4
          djs_oplot, spec.wave[notdata_red], spec.flux[notdata_red,igal], color='red', thick=4
       endif
       djs_oplot, in_filtinfo.weff, bestmab, psym=symcat(6,thick=6), symsize=2.0

       if (nused ne 0) then begin
          oploterror, in_filtinfo[used].weff, mab_used, in_filtinfo[used].fwhm/2.0, $
            mab_used_err, psym=symcat(16), symsize=1.5, color=djs_icolor('dark green'), $
            errcolor=djs_icolor('dark green'), errthick=!p.thick
       endif
       if (nnotused ne 0) then begin
          oploterror, in_filtinfo[notused].weff, mab_notused, in_filtinfo[notused].fwhm/2.0, $
            mab_notused*0.0, psym=symcat(4,thick=6.0), symsize=2.0, color=djs_icolor('blue'), $
            errcolor=djs_icolor('blue'), errthick=!p.thick
       endif

; legend
       legend, textoidl([strtrim(spec.galaxy[igal],2)]), /left, /top, $
         box=0, spacing=1.5, charsize=1.4
    endfor       
    im_plotconfig, psfile=psfile, /psclose, /gzip
    
return
end
