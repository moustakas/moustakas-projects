pro plot_brown_postburst
; jm08nov18nyu - plot up Michael's post-starbursts
    
    path = ages_path(/projects)+'postburst/brown/'
    ispec = read_ages(/ancillary)

; plot the rejected objects
;   readcol, path+'mbrown_reject.list', rejpass, rejaper, $
;     format='I,I', comment='#', /silent
;   help, rejpass
;   rej = speclinefit_locate(ispec,'ages_'+string(rejpass,format='(I3.3)')+$
;     '/'+string(rejaper,format='(I3.3)'))
;   ages_display_spectrum, ispec[rej], /plotobswave, $
;     /postscript, labeltype=5, psname=path+'psb_rejected_08nov18.ps', /pdf

; plot the sample of good objects

    readcol, path+'psb_auto_2008nov19.list', agesid, ra, dec, finalpass, $
      finalaper, format='I,D,D,I,I', comment='#', /silent
    final = speclinefit_locate(ispec,'ages_'+string(finalpass,format='(I3.3)')+$
      '/'+string(finalaper,format='(I3.3)'))
; without the models
    ages_display_spectrum, ispec[final[0]], /plotobswave, plottype=0, nsmooth=0, $
      labeltype=5;, /postscript;, psname=path+'psb_selected_08nov19_nomodels.ps', /pdf
stop
; with the models
    ages_display_spectrum, ispec[final], /plotobswave, plottype=1, nsmooth=0, $
      /postscript, labeltype=5, psname=path+'psb_selected_08nov19.ps', /pdf

stop    
    
return
end
    
