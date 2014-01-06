pro plotbcgsfhs_sfhs, pdf=pdf, build_models=build_models, clobber=clobber
; jm13dec30siena - plot the iSEDfit fitting results

    sample = read_bcgsfhs_sample()
    ncl = n_elements(sample)
    
    bcgpath = bcgsfhs_path()
    if keyword_set(pdf) then begin
       paperpath = bcgsfhs_path()
       suffix = '.ps'
    endif else begin
       paperpath = bcgsfhs_path(/paper)
       suffix = '.eps'
    endelse

    spsmodels = 'fsps_v2.4_miles'
    imf = 'salp'

    outfile = bcgpath+'ssps_'+spsmodels+'_'+imf+'.fits'

    
; make the plot    
    ageline = 0
    Zmetalline = 0
    
    agecolor = 'firebrick'
    Zmetalcolor = 'forest green'

    xrange = [0.6,1.4]
    yrange = [2.0,3.0]

    psfile = paperpath+'bcg_ssps'+suffix
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.3

    pos = im_getposition(nx=5,ny=3,/landscape,yspace=0.1,xspace=0.1,$
      xmargin=[1.0,0.4])

    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)

       if ic gt 9 then begin
          delvarx, xtickname
       endif else begin
          xtickname = replicate(' ',10)
       endelse

       if (ic mod 5) eq 0 then begin
          delvarx, ytickname
       endif else begin
          ytickname = replicate(' ',10)
       endelse
       
       djs_plot, [0], [0], /nodata, position=pos[*,ic], noerase=ic gt 0, $
         xsty=1, ysty=1, yrange=yrange, xrange=xrange, $
         ytickname=ytickname, xtickname=xtickname, ytickinterval=0.3, $
         xtickinterval=0.3
       im_legend, strupcase(cluster), /left, /top, box=0, margin=0, charsize=1.0

; overplot the model grids       
       for iZ = 0, nZmetal-1 do djs_oplot, ssp[iZ].rj, ssp[iZ].ur, $
         line=Zmetalline, color=cgcolor(Zmetalcolor), thick=2
       for ia = 0, nage-1 do djs_oplot, ssp.rj[ia], ssp.ur[ia], $
         line=ageline, color=cgcolor(agecolor), thick=2
    endfor

    xyouts, min(pos[0,*])-0.05, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
      textoidl('u - r'), orientation=90, align=0.5, charsize=1.4, /norm
    xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.09, $
      textoidl('r - J'), align=0.5, charsize=1.4, /norm
    
    im_plotconfig, psfile=psfile, /psclose, pdf=pdf

return
end
    
    
