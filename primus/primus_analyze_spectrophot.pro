pro primus_analyze_spectrophot, ps=ps
; jm08may27nyu - analyze the output from PRIMUS_CHECK_SPECTROPHOT

    path = primus_path()

    if keyword_set(ps) then begin
       debug = 0L
       postthick1 = 4.0
       postthick2 = 3.0
       postthick3 = 8.0
    endif else begin
       im_window, 0, xratio=0.7, yratio=0.6
       postthick1 = 2.0
       postthick2 = 2.0
       postthick3 = 2.0
    endelse

    allmask = file_search(path+'*.fits')
    psname = path+'qa_spectrophot.ps'

    xpage = 8.5 & ypage = 9.8
    im_openclose, psname, postscript=keyword_set(ps), $
      xsize=xpage, ysize=ypage
    
    pagemaker, nx=1, ny=2, width=6.7, height=3.5*[1,1], xmargin=[1.5,0.3], $
      ymargin=[0.8,1.0], xspace=0.0, yspace=1.0, xpage=xpage, ypage=ypage, $
      position=pos, /normal

    xtitle1 = 'B-R (DEEP2)' & ytitle1 = 'B-R (PRIMUS)'
    xtitle2 = 'R-I (DEEP2)' & ytitle2 = 'R-I (PRIMUS)'
    xrange1 = [-0.8,3.8] & yrange1 = xrange1
    xrange2 = [-0.8,2.3] & yrange2 = xrange2

    for imask = 0L, n_elements(allmask)-1L do begin

       splog, 'Reading '+allmask[imask]
       result = mrdfits(allmask[imask],1)

       good = where(result.chi2 gt 0.0,ngood) ; avoid -999.0
       bright = where(result[good].mag lt 22.0,nbright,comp=faint,ncomp=nfaint)
       
       br = -2.5*alog10(result[good].maggies[0]/result[good].maggies[1])
       br_synth = -2.5*alog10(result[good].model_maggies_synth[0]/result[good].model_maggies_synth[1])
       ri = -2.5*alog10(result[good].maggies[1]/result[good].maggies[2])
       ri_synth = -2.5*alog10(result[good].model_maggies_synth[1]/result[good].model_maggies_synth[2])
;      niceprint, br, br_synth, ri, ri_synth
       
       plot, [0], [0], /nodata, charsize=1.8, charthick=postthick2, $
         xthick=postthick1, ythick=postthick1, xtitle=xtitle1, ytitle=ytitle1, $
         title=strtrim(result[0].field,2)+'/'+strtrim(result[0].mask,2)+'/'+$
         strtrim(result[0].date,2), xsty=1, ysty=1, xrange=xrange1, yrange=yrange1, $
         position=pos[*,0]
       oplot, !x.crange, !y.crange, line=0, thick=postthick1
       plotsym, 0, 1.0, fill=0, thick=postthick1
       if (nfaint ne 0L) then djs_oplot, br[faint], br_synth[faint], $
         ps=8, color='orange'
       plotsym, 0, 1.0, fill=1, thick=postthick1
       if (nbright ne 0L) then djs_oplot, br[bright], br_synth[bright], $
         ps=8, color='dark green'

       plot, [0], [0], /nodata, /noerase, charsize=1.8, charthick=postthick2, $
         xthick=postthick1, ythick=postthick1, xtitle=xtitle2, ytitle=ytitle2, $
         xsty=1, ysty=1, xrange=xrange2, yrange=yrange2, $
         position=pos[*,1]
       oplot, !x.crange, !y.crange, line=0, thick=postthick1

       plotsym, 0, 1.0, fill=0, thick=postthick1
       if (nfaint ne 0L) then djs_oplot, ri[faint], ri_synth[faint], $
         ps=8, color='orange'
       plotsym, 0, 1.0, fill=1, thick=postthick1
       if (nbright ne 0L) then djs_oplot, ri[bright], ri_synth[bright], $
         ps=8, color='dark green'

       if (not keyword_set(ps)) then cc = get_kbrd(1)
       
    endfor
    
    im_openclose, postscript=keyword_set(ps), /close

stop       
       
return
end
    
