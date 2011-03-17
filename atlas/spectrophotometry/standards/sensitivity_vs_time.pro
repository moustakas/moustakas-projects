pro sensitivity_vs_time, postscript=postscript
; jm03apr23uofa - written
; jm05jun28uofa - updated to iSPEC2d v2.0 reductions

    psname = 'qaplot_sensitivity_vs_time.ps'

    outpath = atlas_path(/analysis)+'spectrophotometry/standards/'
    logpath = outpath+'qalog_photo/'

    loglist = file_search(logpath+'qalog_sens_4.5*.log',count=logcount)
    photo = iparse_qalogfile(loglist)

    udate = photo[uniq(photo.date,sort(photo.date))].date
    ndate = n_elements(udate)

    mjd = fltarr(ndate)
    meansens = mjd*0.0
    maxsens = mjd*0.0
    medsens = mjd*0.0
    
    for idate = 0L, ndate-1L do begin

       match = where(udate[idate] eq photo.date,nmatch)

       mjd[idate] = djs_mean(fitsdate2jd(photo[match].date))-2400000L
       meansens[idate] = djs_mean(photo[match].senszero)
       medsens[idate] = djs_median(photo[match].senszero)
       maxsens[idate] = max(photo[match].senszero)
              
    endfor
    
    if keyword_set(postscript) then begin
       dfpsplot, outpath+psname, /square, /isolatin1
       postthick = 8.0
    endif else begin
       window, 0, xs=500, ys=500
       postthick = 2.0
    endelse

    plot, mjd, maxsens, ps=-4, xsty=3, ysty=3, xthick=postthick, ythick=postthick, $
      charsize=1.3, charthick=postthick, line=2, thick=postthick, xtitle='Modified Julian Date', $
      ytitle='Sensitivity'
    oplot, mjd, meansens, ps=-6, line=0, thick=postthick

    legend, ['Mean Sensitivity','Maximum Sensitivity'], line=[0,2], psym=[-4,-6], $
      /right, /top, box=0, charsize=1.6, charthick=postthick, thick=postthick
    
    if keyword_set(postscript) then dfpsclose

return
end


