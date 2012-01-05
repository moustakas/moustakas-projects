pro qaplot_alpha_arcfit_night, mike, setup, side, linlist=linlist, qafile=qafile
; jm10jan04ucsd - 

    if (n_elements(linlist) eq 0) then linlist = getenv('XIDL_DIR')+$
      '/Spec/Arcs/Lists/hires_thar.lst'

    ordr_str = mike_getfil('ordr_str',setup,side=side)
    arcs = where(mike.type EQ 'ARC' AND mike.flg_anly NE 0 $
      AND mike.setup EQ setup AND mike.side EQ side,narc)
    thesearcs = strtrim(mike[arcs].img_root,2)

    for iarc = 0, narc-1 do begin
; generate the arc name and the IDL save set name from the raw arcfile
       arc_root = 'Raw/'+thesearcs[iarc]
       arc_fil = 'Arcs/Arc_m'+thesearcs[iarc]
       arc_saveset = 'Arcs/Fits/m'+repstr(thesearcs[iarc],'.fits','_fit.idl')
; restore the IDL save set and read the line-list used in X_ARCFIT
       restore, arc_saveset
       x_arclist, linlist, lines
       lines = lines[sort(lines.wave)]
       lineid = lindgen(n_elements(lines)) ; unique ID number (zero-indexed)
; loop on each order
       for ii = 0L, n_elements(ordr_str)-1L do begin
          if rejstr[ii].ngdf EQ 0L then continue
          these = lindgen(rejstr[ii].ngdf) ; lines used
          gdfit = rejstr[ii].gdfpt[these]
          if (rejstr[ii].nrej NE 0L) then $
            rejpt = rejstr[ii].rejpt[0:rejstr[ii].nrej-1] else $
            rejpt = -1L
          fit = 10^x_calcfit(double(rejstr[ii].gdfpx[these]),FITSTR=all_arcfit[ii])
          goodbad = these*0B+1B ; 1=good, 0=bad
          if (rejpt[0] ne -1L) then goodbad[rejpt] = 0B
; pack into a structure
          out1 = {arc: thesearcs[iarc], order: -1, lineid: -1, pixel: 0.0, $
                  lambda_true: 0.0D, lambda_fit: 0.0D, good: 0B}
          out1 = replicate(out1,n_elements(gdfit))
          out1.order = ordr_str[ii].order
          out1.pixel = rejstr[ii].gdfpx[these]
          out1.lambda_true = lines[gdfit].wave
          out1.lambda_fit = fit
          out1.lineid = lineid[gdfit]
          out1.good = goodbad
          if (iarc eq 0) and (ii eq 0) then $
             out = out1 else out = [out,out1]
       endfor ; close order
    endfor ; close arc lamp

; write out
    if (n_elements(qafile) ne 0) then begin
       outfile = repstr(repstr(qafile,'.ps','.fits'),'qaplot_','')
       im_mwrfits, out, outfile, /clobber
    endif
    
; make the plot    
    good = where(out.good,ngood,comp=bad,ncomp=nbad)
    lam = out.lambda_true
    resid = 1000.0*(out.lambda_fit-out.lambda_true)
    
    xrange = minmax(out.lambda_true)
;   xrange = [5600,5700]
;   yrange = max(abs(out.lambda_true-out.lambda_fit))*[-1,1]*1E3 ; [mA]
    yrange = 20*[-1,1]

    im_plotconfig, 0, pos, psfile=qafile, xmargin=[1.2,0.3]
; in the first plot show all the lines used in the fitting (plotting
; repeat lines multiple times), but distinguishing between lines that
; were rejected/retained
    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xtitle='Wavelength (\AA)', ytitle='Residuals (m\AA)'
    djs_oplot, lam[good], resid[good], psym=symcat(16), symsize=0.4
    djs_oplot, lam[bad], resid[bad], psym=symcat(9,thick=3), symsize=0.5, color='red'
    im_legend, ['N_{arc}='+string(narc,format='(I0)'),$
      'N_{good}='+string(ngood,format='(I0)'),$
      'N_{bad}='+string(nbad,format='(I0)')], $
      /right, /top, box=0, charsize=1.4
; in the second plot group the individual lines together; make the
; symbol size proportional to the fraction of the time that the line
; was rejected
    allid = out.lineid
    uindx = uniq(allid,sort(allid))
    id = allid[uindx]
    nid = n_elements(id)
    badfrac = fltarr(nid)
    
    for jj = 0, nid-1 do begin
       these = where(id[jj] eq allid,nthese)
       badfrac[jj] = total(out[these].good eq 0)/float(nthese)
;      symsize = total(out[these].good eq 1)/float(nthese)
;      zero = where(symsize eq 0.0,nzero)
;      if (nzero ne 0) then symsize[zero] = min(symsize)
    endfor

    xrange = minmax(out[uindx].lambda_true)
    yrange = [0.0,max(badfrac)]
    
    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xtitle='Wavelength (\AA)', ytitle='F = Fraction of Time Line Rejected'
    djs_oplot, out[uindx].lambda_true, badfrac, psym=symcat(16), $
      symsize=1.2
    good = where(badfrac lt 0.05,ngood,comp=bad,ncomp=nbad)
    im_legend, [$
      'N_{lines}='+string(n_elements(lines),format='(I0)'),$
      'N(F<0.05)='+string(ngood,format='(I0)'),$
      'N(F>0.05)='+string(nbad,format='(I0)')], $
      /right, /top, box=0, charsize=1.4

    im_plotconfig, psfile=qafile, /psclose, /pdf
    
return
end
    
