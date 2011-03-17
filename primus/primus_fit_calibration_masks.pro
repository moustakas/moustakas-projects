pro primus_fit_calibration_masks, dofit=dofit, doplot=doplot
; jm08sep18nyu - fit all the PRIMUS calibration masks

; echo "primus_fit_calibration_masks" | idl > & log & 

    rerun = '0216'
    rootname = ['d23a0043','d23a0045','vvds0184']

    for jj = 0L, n_elements(rootname)-1L do begin

;      index = 8
       index = [3,4,8,9,10,11]
       if keyword_set(dofit) then begin
          primus_fit_redshift, rootname[jj], rerun, /nostar, $
            /noagn, /nopowerlaw, index=index, /clobber
       endif

       if keyword_set(doplot) then begin

          primus_read_1dinputs, rootname[jj], rerun, extract=extract, $
            slits=slits, oned=oned, photoinfo=photo
          if (n_elements(index) eq 0L) then begin
             index = where((slits.currz gt 0.0) and (oned.zbest_gal gt 0.0))
          endif
;         splog, slits[index].currz
;         niceprint, oned[index].zmin_gal, oned[index].zmin_photoz_gal

          im_plotfaves, /post
          dfpsplot, '~/d23a0043_photoz_test.ps', /color, /landscape
          for kk = 0L, n_elements(index)-1L do begin
             chi2min = oned[index[kk]].chi2best_gal
             djs_plot, oned[index[kk]].zgrid_gal, oned[index[kk]].chi2_photoz_gal-chi2min, $
               xr=[0,1.2], xsty=3, ysty=3, ps=10, color='dark red', $
               xtitle='Redshift', ytitle='\chi^{2}-\chi^{2}_{min}', $
               title='Object '+string(index[kk],format='(I4.4)'), $
               yrange=[0.0,(max(oned[index[kk]].chi2_photoz_gal)>max(oned[index[kk]].chi2_gal))-chi2min]
             djs_oplot, oned[index[kk]].zgrid_gal, oned[index[kk]].chi2_gal-chi2min, ps=10, color='blue'
             djs_oplot, slits[index[kk]].currz*[1,1], !y.crange, color='dark green', thick=10.0
             djs_oplot, oned[index[kk]].zbest_photoz_gal*[1,1], !y.crange, color='dark red', thick=10.0, line=5
             djs_oplot, oned[index[kk]].zbest_gal*[1,1], !y.crange, color='blue', thick=10.0, line=5
          endfor
          dfpsclose & im_plotfaves

stop          
          
          plot, slits[index].currz, oned[index].zbest_gal, ps=4, $
            xr=[0,1.1], yr=[0,1.1]
          plot, slits[index].currz, oned[index].zbest_photoz_gal, ps=4, $
            xr=[0,1.1], yr=[0,1.1]

       endif

stop       
       
    endfor
       
stop    
    
return
end
    
