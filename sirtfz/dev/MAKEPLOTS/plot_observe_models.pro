pro plotfluxtitles

    xyouts, [0.52,0.52], [0.05,0.05], 'log z', align=0.5, charthick=2.0, $
      charsize=1.8, /normal
    xyouts, [0.02,0.02], [0.52,0.52], 'log '+textoidl('f_{\nu}')+' (mJy)', align=0.5, charthick=2.0, $
      charsize=1.8, /normal, orientat=90
    
return
end

pro plotfluxes, z, sed, position, xy, ylim, jstart

    for j = jstart, jstart+3L do begin

       k = j-jstart
       pos = position[*,k]
    
       case k of

          0L : plot, [min(z)/1.1,max(z)*1.3], [ylim[0,k],ylim[1,k]], /nodata, /noerase, xsty=3, ysty=3, $
            /ylog, /xlog, position=pos
          1L : begin
             plot, [min(z)/1.1,max(z)*1.3], [ylim[0,k],ylim[1,k]], /nodata, /noerase, xsty=3, ysty=7, $
               /ylog, /xlog, position=pos, ytickname=replicate(' ',10)
             axis, yaxis=1, ysty=3, /ylog
          endcase
          2L : plot, [min(z)/1.1,max(z)*1.3], [ylim[0,k],ylim[1,k]], /nodata, /noerase, xsty=3, ysty=3, $
            /ylog, /xlog, position=pos, xtickname=replicate(' ',10)
          3L : begin
             plot, [min(z)/1.1,max(z)*1.3], [ylim[0,k],ylim[1,k]], /nodata, /noerase, xsty=3, ysty=7, $
               /ylog, /xlog, position=pos, xtickname=replicate(' ',10)
             axis, yaxis=1, ysty=3, /ylog
          endcase             

       endcase

       oplot, z, sed[j].f3_6mu, color=2, line=0
       oplot, z, sed[j].f4_5mu, color=3, line=1
       oplot, z, sed[j].f5_8mu, color=4, line=2
       oplot, z, sed[j].f8mu, color=5, line=3
       oplot, z, sed[j].f24mu, color=6, line=4
       oplot, z, sed[j].f70mu, color=7, line=5
       oplot, z, sed[j].f160mu, color=16, line=6      

       xyouts, [xy[0,k],xy[0,k]], [xy[1,k],xy[1,k]], sed[j].model_name, align=0.5, $
         charthick=1.5, charsize=1.3, /normal
       xyouts, [xy[0,k],xy[0,k]], [xy[1,k]-0.03,xy[1,k]-0.03], $
         string(sed[j].lum60_sun,format='(F5.2)'), align=0.5, $
         charthick=1.5, charsize=1.3, /normal

    endfor

    plotfluxtitles
    
return
end

;pro plotcolors
;
;return
;end

pro plot_observe_models, ps=ps

; jm00oct30uofa
; plot up the various SIRTF fluxes of the template SEDs as a function
; of redshift.  reads the output from TOY_UNIVERSE.PRO.
    
    colortable2
    plotfaves, pthick=1.8, charsize=1.8, charthick=1.9
    plotsym, 0, 1, /fill

    path = (sirtf_datapath())[0]
    sed = sread(path+'data/observe_models.dat') ; restore the "observations"

    zsample = redshift_array() ; redshift structure
    z = zsample.zarray
    nz = n_elements(z)
    
    position = [ [0.1,0.1,0.52,0.52],  $
                 [0.52,0.1,0.94,0.52], $
                 [0.1,0.52,0.52,0.94], $
                 [0.52,0.52,0.94,0.94] ]                 

    xy = [ [0.43,0.47], $
           [0.85,0.47], $
           [0.43,0.89], $
           [0.85,0.89] ]

; ----------------------------------------------------------------------
    if keyword_set(ps) then begin
       ps_open, path+'plots/fzmodels1', /portrait, /ps_fonts
       device, /inches, /times
    endif else window, 0, xs=550, ys=550
    jstart = 0L
    ylim = [ [3E-7,2E-2], $
             [5E-7,3E-2], $
             [5E-7,4E-2], $
             [5E-7,1E-1] ]
    plotfluxes, z, sed, position, xy, ylim, jstart ; models 1-4
    if keyword_set(ps) then ps_close
; ----------------------------------------------------------------------

; ----------------------------------------------------------------------
    if keyword_set(ps) then begin
       ps_open, path+'plots/fzmodels2', /portrait, /ps_fonts
       device, /inches, /times
    endif else window, 1, xs=550, ys=550
    jstart = 4L
    ylim = [ [3E-7,0.5], $
             [3E-7,1], $
             [3E-7,1], $
             [3E-7,0.5] ]
    plotfluxes, z, sed, position, xy, ylim, jstart ; models 5-8
    if keyword_set(ps) then ps_close
; ----------------------------------------------------------------------

; ----------------------------------------------------------------------
    if keyword_set(ps) then begin
       ps_open, path+'plots/fzmodels3', /portrait, /ps_fonts
       device, /inches, /times
    endif else window, 2, xs=550, ys=550
    jstart = 8L
    ylim = [ [5E-5,100], $
             [5E-5,5000], $
             [1E-4,5000], $
             [5E-5,5000] ]
    plotfluxes, z, sed, position, xy, ylim, jstart ; models 9-12
    if keyword_set(ps) then ps_close
; ----------------------------------------------------------------------

; ----------------------------------------------------------------------
    if keyword_set(ps) then begin
       ps_open, path+'plots/fzmodels4', /portrait, /ps_fonts
       device, /inches, /times
    endif else window, 3, xs=550, ys=550
    jstart = 12L
    ylim = [ [5E-5,5000], $
             [5E-5,5000], $
             [1E-4,5000], $
             [5E-5,5000] ]
    plotfluxes, z, sed, position, xy, ylim, jstart ; models 13-16
    if keyword_set(ps) then ps_close
; ----------------------------------------------------------------------

    plotfaves, /restore
    
return
end
