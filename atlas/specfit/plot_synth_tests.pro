pro plot_synth_tests, agebins=agebins, postscript=postscript
; jm03sep18uofa
    
    restore, 'synth_tests.idlsave'
    snratio = result[0].snratio
    
    if n_elements(agebins) eq 0L then agebins = [100.0,2000.0] ; [Myr]

    nagebins = n_elements(agebins)+1
    natlas = n_elements(result)
    nmonte = n_elements(result[0].chi2_fit)

    age = eigeninfo.type/1E6          ; template ages [Myr]
;   age = eigeninfo.template_type/1E6 ; template ages [Myr]

    stats = {$
      mean:     fltarr(nagebins), $
      variance: fltarr(nagebins), $
      sigma:    fltarr(nagebins), $
      truemean: fltarr(nagebins), $
      resid:    fltarr(nagebins)}
    stats = replicate(stats,natlas)

    strage = strarr(nagebins)
    
    for i = 0L, nagebins-1L do begin

       case i of
          0L: begin
             indx = where(age le agebins[i],nindx)
             strage[i] = 'Age < '+string(agebins[i],format='(I0)')+' Myr'
          end
          nagebins-1L: begin
             indx = where(age ge agebins[i-1L],nindx)
             strage[i] = 'Age > '+string(agebins[i-1L],format='(I0)')+' Myr'
          end
          else: begin
             indx = where((age gt agebins[i-1L]) and (age lt agebins[i]),nindx)
             strage[i] = string(agebins[i-1L],format='(I0)')+' > Age < '+$
               string(agebins[i],format='(I0)')+' Myr'
          endelse
       endcase

       stats.truemean[i] = total(result.coeff_true[indx],1)

       bintotal = total(result.coeff[indx,*,*],1)
       
       stats.mean[i]     = total(bintotal,1)/nmonte
       smean = transpose(stats.mean[i] # (fltarr(nmonte)+1))
       
       stats.variance[i] = total((bintotal-smean)^2.0,1)/(nmonte-1)
       stats.sigma[i] = sqrt(stats.variance[i])
          
    endfor

    stats.resid = stats.mean-stats.truemean

; generate plots
    
    yrange = 100*(abs(min(stats.resid))>max(stats.resid))*[-1.1,1.1]
    xrange = 100*minmax(stats.mean)
    
    plot, [0], [0], /nodata, xsty=3, ysty=3, $
      yrange=yrange, xthick=2.0, ythick=2.0, charsize=2.0, $
      charthick=2.0, xtitle='Mean Percent Contribution', $
      ytitle='Percent Residuals', xrange=xrange
    djs_oplot, !x.crange, [0,0], line=0, thick=3.0

    plotsym, 0, 1
    djs_oplot, 100*stats.mean[0], 100*stats.resid[0], ps=8
    djs_oplot, 100*stats.mean[1], 100*stats.resid[1], ps=6
    djs_oplot, 100*stats.mean[2], 100*stats.resid[2], ps=7

    legend, ['Young','Intermediate','Old'], psym=[8,6,7], /right, $
      /bottom, box=0, charsize=2.0, charthick=2.0
;   legend, strage, psym=[8,6,7], /right, /bottom, box=0, charsize=2.0, charthick=2.0

return
end    
