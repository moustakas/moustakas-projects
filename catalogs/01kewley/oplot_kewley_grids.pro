pro oplot_kewley_grids, plotnumber=plotnumber, _extra=extra
; jm02jul17uofa
; can add continuous and starburst99 models later

; this routine assumes:
;    [O III] = [O III] 5007    
;    [N II]  = [N II]  6584
    
    if n_elements(plotnumber) eq 0L then plotnumber = 1L
    
; correction factors    
    
    nratio = 3.054
    oratio = 2.984
    
    ncor = nratio/(1.0+nratio)
    ocor = oratio/(1.0+oratio)

    grids = read_kewley_grids(Z=Z,U=U,_extra=extra)

    nZ = n_elements(Z)
    nU = n_elements(U)
    
    case plotnumber of

       1L: begin
          x = grids.oiii_oii    ; [O III]/[O II]
          y = grids.oii_h_alpha ; [O II]/Ha
       end
       2L: begin
          x = grids.nii_oii     ; [N II]/[O II]
          y = grids.oii_h_alpha ; [O II]/Ha
       end
       3L: begin
          x = grids.nii_sii     ; [N II]/[S II]
          y = grids.oii_h_alpha ; [O II]/Ha
       end
       4L: begin
          x = grids.nii_oii  ; [N II]/[O II]
          y = grids.oiii_oii ; [O III]/[O II]
       end
       5L: begin
          x = grids.nii_h_alpha ; [N II]/Ha
          y = grids.oii_h_alpha ; [O II]/Ha
       end
       else:
    endcase
    
    for i = 0L, nZ-1L do begin
       zidx = i*nZ + lindgen(nZ)
       oplot, x[zidx], y[zidx], thick=1.5, color=djs_icolor('dark green')
       for j=0,nU-1 do begin
          qidx = i + lindgen(nU)*nZ
          oplot, x[qidx], y[qidx], thick=1.5, color=djs_icolor('dark red')
       endfor
    endfor

; label the first and last grid points

    xyouts, x[nU-1L], y[nU-1L], textoidl('Z/Z'+sunsymbol()+'=')+strn(ztext[0],length=3), $
      align=0.0, /data, charsize=1.4, charthick=5.0
    xyouts, x[(nZ-1L)*nU+nU-1L], y[(nZ-1L)*nU+nU-1L], textoidl('Z/Z'+sunsymbol()+'=')+strn(ztext[nZ-1L],length=3), $
      align=0.0, /data, charsize=1.4, charthick=5.0

;   for i = 0L, nZ-2L do xyouts, x[i*nU+nU-1L], y[i*nU+nU-1L], strn(ztext[i],length=3), $
;     align=0.0, /data, charsize=1.4, charthick=2.0

    xyouts, x[0], y[0], 'log U='+strn(qtext[0],format='(F5.2)'), /data, align=1.0, $
      charsize=1.4, charthick=5.0
    xyouts, x[nU-1L], y[nU-1L], 'log U='+strn(qtext[nU-1L],format='(F5.2)')+'; ', /data, align=1.0, $
      charsize=1.4, charthick=5.0
    
;   xyouts, x[0], y[0], 'q='+strn(qtext[0],format='(E7.1)'), /data, align=1.0, $
;     charsize=1.4, charthick=2.0
;   xyouts, x[nU-1L], y[nU-1L], 'q='+strn(qtext[nU-1L],format='(E7.1)')+'; ', /data, align=1.0, $
;     charsize=1.4, charthick=2.0
    
;   for i = 1L, nU-1L do xyouts, x[i], y[i], strn(qtext[i],format='(E7.1)'), /data, align=1.0, $
;     charsize=1.4, charthick=2.0

; label all the grid points
    
;    for i = 0L, nZ-2L do xyouts, x[i*nU+nU-1L], y[i*nU+nU-1L], strn(ztext[i],length=3), $
;      align=0.0, /data, charsize=1.4, charthick=2.0
;    xyouts, x[(nZ-1L)*nU+nU-1L], y[(nZ-1L)*nU+nU-1L], 'Z='+strn(ztext[nZ-1L],length=3), $
;      align=0.0, /data, charsize=1.4, charthick=2.0
;
;    xyouts, x[0], y[0], 'q='+strn(qtext[i],format='(E8.2)'), /data, align=1.0, $
;      charsize=1.4, charthick=2.0
;    for i = 1L, nU-1L do xyouts, x[i], y[i], strn(qtext[i],format='(E8.2)'), /data, align=1.0, $
;      charsize=1.4, charthick=2.0

return
end    



