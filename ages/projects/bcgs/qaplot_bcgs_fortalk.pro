pro qaplot_bcgs_fortalk, all=all
; jm10jul23ucsd - make some pretty plots for talks

    bcgspath = ages_path(/projects)+'bcgs/'
    paramfile = bcgspath+'bcgs_isedfit.par'

    params = read_isedfit_paramfile(paramfile)
    params.sfhgrid = 2 ; 2 ; solar!

; *first* restore the best-fitting models and *then* build the
; filenames 
    model = isedfit_restore(paramfile,isedfit,params=params,$
      iopath=bcgspath,index=index,silent=silent)
    ngal = n_elements(isedfit)
    fp = isedfit_filepaths(params,iopath=bcgspath,ngalaxy=ngal)

;; read all the models
;    nchunk = n_elements(fp.isedfit_models_chunkfiles)
;    for ichunk = 0, nchunk-1 do begin
;       chunkfile = fp.modelspath+fp.isedfit_models_chunkfiles[ichunk]+'.gz'
;       if (keyword_set(silent) eq 0) then splog, 'Reading '+chunkfile
;       modelgrid1 = mrdfits(chunkfile,1,/silent)
;       modelgrid1 = struct_trimtags(temporary(modelgrid1),$
;         except=['MODELMAGGIES'])
;       if (n_elements(modelgrid) eq 0) then $
;         modelgrid = temporary(modelgrid1) else $
;           modelgrid = [temporary(modelgrid),temporary(modelgrid1)]
;    endfor
;    age = modelgrid.modelage
;    nmodel = n_elements(modelgrid)
;    nage = n_elements(modelgrid[0].modelage)
;    
;; read the chi2 grid (assumes just one gchunk file!)
;    chi2gridfile = fp.modelspath+fp.chi2grid_gchunkfiles[0]+'.gz'
;    splog, 'Reading '+chi2gridfile
;    chi2grid = mrdfits(chi2gridfile,1)
;    chi2grid = reform(chi2grid,nage,nmodel,ngal)
    
; restore some filter info
    filterlist = strtrim(params.filterlist,2)
    filtinfo = im_filterspecs(filterlist=filterlist)
    nfilt = n_elements(filterlist)

; generate one plot for each galaxy
    xtitle1 = 'Observed Wavelength (\AA)'
    xtitle2 = 'Rest Wavelength (\AA)'
    ytitle1 = 'm_{AB}'

    if keyword_set(all) then begin
       psfile = bcgspath+'fortalk/bcgs_allseds.ps'
       im_plotconfig, 8, pos, psfile=psfile, bits_per_pixel=24
    endif
;   for ii = 0, 4 do begin

;   for igal = 2, 2 do begin
    for igal = 0L, ngal-1L do begin
       if (keyword_set(all) eq 0) then begin
          psfile = bcgspath+'fortalk/bcg'+string(igal,format='(I2.2)')+'_sed.ps'
          im_plotconfig, 8, psfile=psfile, /keynote
       endif

; redshift and best-fit model
; ToDo: add IGM absorption!!       
       zobj = isedfit[igal].zobj
       wave = model[igal].wave    ; [A]
       flux = model[igal].flux    ; [AB]

; make the plot       
       xrange1 = [min(filtinfo.weff-1.3*filtinfo.fwhm),$
         max(filtinfo.weff+2*filtinfo.fwhm)]
       xrange2 = xrange1/(1.0+zobj)
       if (isedfit[igal].chi2 ge 1E6) then begin
          djs_plot, [0], [0], /nodata, xsty=9, ysty=1, xrange=xrange1, $
            yrange=yrange, /xlog, xtitle=xtitle1, ytitle=ytitle1, $
            ytickname=replicate(' ',10), ymargin=[4,3], xtickinterval=1000
          axis, /xaxis, xsty=1, xtitle=textoidl(xtitle2), xrange=xrange2, /xlog
          legend, ['No mass estimate available'], /left, /top, $
            box=0, spacing=1.5, charsize=1.6
          
          label = [strtrim(galaxy[igal],2),'z = '+string(zobj,format='(F6.4)')]
          legend, label, /right, /bottom, box=0, spacing=1.5, charsize=1.6
       endif else begin
          xrange1[0] = xrange1[0]>(900.0*(1+zobj))
          xrange1[1] = xrange1[1]<max(wave)
          get_element, wave, xrange1, xx

          weff = filtinfo.weff
          hwhm = filtinfo.fwhm/2.0
          bestmab = -2.5*alog10(isedfit[igal].bestmaggies)

          yrange = fltarr(2)
          yrange[0] = ((max(bestmab)>max(flux[xx[0]:xx[1]]))*1.05)<30.0
          yrange[1] = (min(bestmab)<min(flux[xx[0]:xx[1]]))*0.93

          djs_plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
            xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, /xlog, $
            ymargin=[4,3], xtickinterval=1000
          axis, /xaxis, xsty=1, xtitle=textoidl(xtitle2), xrange=xrange2, /xlog

          djs_oplot, wave, flux, line=0, color='grey'
          djs_oplot, weff, bestmab, psym=symcat(6,thick=6), symsize=2.5

; overplot the data; distinguish between three different cases, based
; on the input photometry
          used = where((isedfit[igal].maggies gt 0.0) and $ ; used in the fitting
            (isedfit[igal].ivarmaggies gt 0.0),nused)
          notused = where((isedfit[igal].maggies gt 0.0) and $ ; not used in the fitting
            (isedfit[igal].ivarmaggies eq 0.0),nnotused)
          nodata = where((isedfit[igal].maggies eq 0.0) and $ ; no measurement
            (isedfit[igal].ivarmaggies eq 0.0),nnodata)

          if (nused ne 0L) then begin
             mab = maggies2mag(isedfit[igal].maggies[used],$
               ivar=isedfit[igal].ivarmaggies[used],magerr=mab_err)
             oploterror, weff[used], mab, hwhm[used], mab_err, psym=symcat(16), $
               symsize=2.0, color=djs_icolor('dark green'), $
               errcolor=djs_icolor('dark green'), errthick=!p.thick
          endif

          if (nnotused ne 0L) then begin
             mab = maggies2mag(isedfit[igal].maggies[notused])
             oploterror, weff[notused], mab, hwhm[notused], mab*0.0, $
               psym=symcat(4,thick=6.0), symsize=3.0, color=djs_icolor('red'), $
               errcolor=djs_icolor('red'), errthick=!p.thick
          endif

;; make the inset
;          chi2 = chi2grid[*,*,igal].chi2
;          mass = chi2grid[*,*,igal].mass
;          like = exp(-(chi2-isedfit[igal].chi2))
;          contour, like, age, alog10(mass>1E8)
          
          label1 = 'z = '+string(zobj,format='(F6.4)')
          label2 = 'log (M/M'+sunsymbol()+') = '+strtrim(string(isedfit[igal].mass_avg,format='(F12.2)'),2)
          im_legend, label2, /left, /top, box=0, spacing=1.7, charsize=1.6;, margin=0
          im_legend, label1, /right, /bottom, box=0, spacing=1.5, charsize=1.6;, margin=0
       endelse
       if (keyword_set(all) eq 0) then $
         im_plotconfig, psfile=psfile, /psclose, /pdf, /keynote
    endfor 
    if keyword_set(all) then im_plotconfig, psfile=psfile, /gzip, /psclose

return
end
