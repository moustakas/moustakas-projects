pro qaplot_id_arclines, number
; jm10jan04ucsd - make QAplots
    
    arcfile = 'Arcs/Fits/mr0001_fit.idl'
    splog, 'Restoring '+arcfile
    restore, arcfile

    case number of
       0: begin
          linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/alpha_thar_custom_01.lst'
          psfile = getenv('DEEP2_ALPHA_DIR')+'/qaplots/qaplot_arclines_murphy.ps'
       end
       1: begin
          linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/hires_thar.lst'
          psfile = getenv('DEEP2_ALPHA_DIR')+'/qaplots/qaplot_arclines_hires.ps'
       end
    endcase

    x_arclist, linlist, lines
    lines = lines[sort(lines.wave)]
    
    im_plotconfig, 8, pos, psfile=psfile, ymargin=[2.3,1.2], $
      width=8.0
    
    nordr = n_elements(guess_ordr)
    for ii = 0, nordr-1 do begin
       sz = size(sv_aspec, /dimen)
       wv = x_calcfit(dindgen(sz[0]),fitstr=all_arcfit[ii])
       if (wv[0] lt 10.0) then wv = 10.0^wv
       flux = sv_aspec[*,ii]*1D-4

       xrange = minmax(wv)
       yrange = [min(flux),max(flux)*0.6]

       these = lindgen(rejstr[ii].ngdf) ; lines used
       gdfit = rejstr[ii].gdfpt[these]
       if (rejstr[ii].nrej eq 0) then rejpt = -1 else $
         rejpt = rejstr[ii].rejpt[0:rejstr[ii].nrej-1]
       goodbad = these*0B+1B    ; 1=good, 0=bad
       if (rejpt[0] ne -1) then goodbad[rejpt] = 0B

       extend_color = replicate(djs_icolor('default'),n_elements(these))
       bad = where(goodbad eq 0,nbad)
       if (nbad ne 0) then extend_color[bad] = djs_icolor('red')
       
       mlines = lines[gdfit]

       label = strtrim(string(mlines.wave,format='(F12.3)'),2)
;      label = strtrim(mlines.name,2)
       im_lineid_plot, wv, flux, mlines.wave, label, $
         xsty=3, ysty=3, xtitle=textoidl('Wavelength (\AA)'), $
         ytitle=textoidl('Counts (\times 10^{-4})'), xrange=xrange, $
         yrange=yrange, extend_thick=2.0, extend_color=extend_color, $
         label_color=extend_color, position=pos, lcharthick=1.0, $
         lcharsize=1.2
;      djs_plot, wv, flux, xsty=3, ysty=3, position=pos, $
;        xtitle='Wavelength (\AA)', ytitle='Counts (\times 10^{-4})'
       im_legend, 'Order='+string(guess_ordr[ii],format='(I0)'), $
         /right, /top, box=0, charsize=1.5
    endfor

    im_plotconfig, psfile=psfile, /psclose, /gzip

return
end
    
