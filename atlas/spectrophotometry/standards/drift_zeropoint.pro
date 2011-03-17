pro drift_zeropoint, postscript=postscript, paper=paper
; jm03arp24uofa
; jm04apr30uofa - new computation of this quantity

    if keyword_set(paper) then begin
       psname = 'qaplot_drift_zeropoint.eps'
       postscript = 1L
       encapsulated = 1L
    endif else begin
       psname = 'qaplot_drift_zeropoint.ps'
       encapsulated = 0L
    endelse

    rootpath = atlas_path(/analysis)+'spectrophotometry/standards/'

    if keyword_set(paper) then $
      outpath = atlas_path(/papers)+'atlas/FIG_ATLAS/' else outpath = rootpath
    logpath = rootpath+'qalog_photo/'

    driftlogs = file_search(logpath+'*drift*.log')
    logs = repstr(driftlogs,'drift','4.5')
    nlogs = n_elements(logs)

    zero = {drift: ptr_new(), drifterr: 0.0, sens: 0.0, $
      senserr: 0.0, date: '', zpt: ptr_new(), zpterr: ptr_new(), $
      seeing: ptr_new(), energy: ptr_new(), result: 0.0, $
      result_err: 0.0, driftnum: 0L, photonum: 0L}
    zero = replicate(zero,nlogs)

    for j = 0L, nlogs-1L do begin
    
       drift = iparse_qalogfile(driftlogs[j])
       ndrift = n_elements(drift)
       
       photo = iparse_qalogfile(logs[j])
       nstds = n_elements(photo)
       
       zero[j].drift = ptr_new(drift.senszero) 
       
       zero[j].date = drift[0].date

       djs_iterstat, photo.senszero, mean=mn, sigma=sig, mask=mask
       zero[j].sens = mn
       zero[j].senserr = sig

       zero[j].zpt = ptr_new(drift.senszero-zero[j].sens)
       zero[j].zpterr = ptr_new(replicate(zero[j].senserr,ndrift))

       zero[j].seeing = ptr_new(replicate(djs_mean(photo[where(mask)].seeing),ndrift))
       zero[j].energy = ptr_new(replicate(djs_mean(photo[where(mask)].energy),ndrift))

       zero[j].result = djs_mean(drift.senszero)-djs_mean(photo.senszero)
       zero[j].result_err = sqrt((djsig(drift.senszero)/sqrt(ndrift))^2+(djsig(photo.senszero)/sqrt(nstds))^2)
       zero[j].driftnum = ndrift
       zero[j].photonum = nstds

    endfor

    srt = sort(fitsdate2jd(zero.date))
    zero = zero[srt]
    
    for j = 0L, nlogs-1L do begin
       if (j eq 0L) then begin
          zpt = *zero[j].zpt 
          zpterr = *zero[j].zpterr
          seeing = *zero[j].seeing
          energy = *zero[j].energy
       endif else begin
          zpt = [zpt,*zero[j].zpt]
          zpterr = [zpterr,*zero[j].zpterr]
          seeing = [seeing,*zero[j].seeing]
          energy = [energy,*zero[j].energy]
       endelse
    endfor

    xtitle = 'Observing Run Number'
    xrange = [0.1,nlogs+0.9]

    meanzpt = djs_mean(zero.result)
    sigzpt = stddev(zero.result)
    x = lindgen(nlogs)

    if keyword_set(postscript) then begin
       dfpsplot, outpath+psname, /square, encapsulated=encapsulated
       postthick = 8.0
    endif else begin
       window, 0, xs=500, ys=500
       postthick = 2.0
    endelse

    plotsym, 8, 1.2, fill=0, thick=postthick
    ploterror, lindgen(nlogs)+1L, zero.result, zero.result_err, ps=8, xsty=3, ysty=3, $
      syms=2.5, xthick=postthick, ythick=postthick, thick=postthick, errthick=postthick, $
      charsize=1.8, charthick=postthick, $
      xtitle=xtitle, ytitle='Sensitivity Difference [20" Scan - 4.5" Slit, mag]', $
      yrange=[0.05,0.20], xrange=xrange
    for i = 0L, nlogs-1L do begin
       xyouts, i+1, zero[i].result-zero[i].result_err-0.005, string(zero[i].driftnum,format='(I0)')+'/'+$
         string(zero[i].photonum,format='(I0)'), /data, charthick=postthick, charsize=1.0, $
         align=0.5
;      xyouts, i+1, zero[i].result-zero[i].result_err-0.01, strmid(zero[i].date,0,7), $
;        /data, charthick=postthick, charsize=1.0, align=0.5
    endfor
      
    oplot, !x.crange, meanzpt*[1.0,1.0], line=0, thick=postthick
    oplot, !x.crange, (meanzpt+sigzpt)*[1.0,1.0], line=2, thick=postthick
    oplot, !x.crange, (meanzpt-sigzpt)*[1.0,1.0], line=2, thick=postthick

    legtxt = string(meanzpt,format='(F4.2)')+' \pm '+string(sigzpt,format='(F4.2)')+' mag'
;   legend, textoidl(legtxt), /right, /top, charsize=1.8, charthick=postthick, box=0
;   legend, textoidl(['\Delta = '+strn(meanzpt,format='(F4.2)')+' mag',$
;     '\sigma_{\Delta} = '+strn(sigzpt,format='(F5.3)')+' mag']), /left, $
;     /top, charsize=2.0, charthick=postthick, box=0

    if keyword_set(postscript) then dfpsclose

    niceprint, zero.date

    ptr_free, zero.drift, zero.zpt, zero.zpterr, zero.seeing, zero.energy

return
end    
