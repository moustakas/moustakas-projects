;+
; NAME:
;       PLOT_KEWLEY_GRIDS
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;       log12oh   - 
;       siioffset - shift the [S II] 6716,6731 ratio to [S II] 6716
;                   assuming an electron density of ~100/cm^3 (not
;                   generalized: see code below)
;
; OUTPUTS:
;
; PROCEDURES USED:
;       READ_KEWLEY_GRIDS()
;
; COMMENTS:
;       Add continuous and starburst99 models later. 
;
;       This routine assumes that [O III] = [O III] 5007 and [N II]  =
;       [N II] 6584.
;
;       A [S II] 6716 / [S II] 6731 ratio of 1.3 gives an electron
;       density of 135/cm^3.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 July 17, U of A
;       jm04januofa - major updates and generalizations 
;       jm04dec07uofa - more plots
;       jm05mar02uofa - added LOG12OH, ULABELCOLOR, and ZLABELCOLOR
;                       inputs 
;       jm05mar14uofa - added UCHARSIZE and ZCHARSIZE inputs
;       jm05oct23uofa - added SIIOFFSET keyword
;-

pro plot_kewley_grids, plotnumber=plotnumber, overplot=overplot, labeltype=labeltype, $
  Umin=Umin, Umax=Umax, Zmin=Zmin, Zmax=Zmax, noZgrid=noZgrid, noUgrid=noUgrid, $
  Zgridcolor=Zgridcolor, Ugridcolor=Ugridcolor, Ulabelcolor=Ulabelcolor, $
  Zlabelcolor=Zlabelcolor, Ucharsize=Ucharsize, Zcharsize=Zcharsize, $
  postscript=postscript, hiiregions=hiiregions, log12oh=log12oh, $
  _extra=extra, siioffset=siioffset

    light = 2.99792458D10        ; speed of light [cm/s]

    if n_elements(plotnumber) eq 0L then plotnumber = 1L
    if n_elements(Zgridcolor) eq 0L then Zgridcolor = 'dark red' 
    if n_elements(Ugridcolor) eq 0L then Ugridcolor = 'dark green'
    if n_elements(Ulabelcolor) eq 0L then Ulabelcolor = 'dark green'
    if n_elements(Zlabelcolor) eq 0L then Zlabelcolor = 'dark red'
    if n_elements(Ucharsize) eq 0L then Ucharsize = 1.0
    if n_elements(Zcharsize) eq 0L then Zcharsize = 1.0

    if keyword_set(hiiregions) then regions = read_hii_regions(/sample)
    
; correction factors    
    
    nratio = 3.054
    oratio = 2.984
    
    ncor = nratio/(1.0+nratio)
    ocor = oratio/(1.0+oratio)

    Zsun = 8.69

    siiratio = alog10(1.3)
    
    grids = read_kewley_grids(Z=Z,U=U,_extra=extra)

; allow the user to only plot certain ranges of the metallicity and
; ionization parameter
    
    if (n_elements(Umax) ne 0L) then begin
       w = where(U le Umax,nw)
       if (nw ne 0L) then begin
          U = U[w]
          grids = grids[*,w]
       endif
    endif
    
    if (n_elements(Umin) ne 0L) then begin
       w = where(U ge Umin,nw)
       if (nw ne 0L) then begin
          U = U[w]
          grids = grids[*,w]
       endif
    endif

    if (n_elements(Zmax) ne 0L) then begin
       w = where(Z le Zmax,nw)
       if (nw ne 0L) then begin
          Z = Z[w]
          grids = grids[w,*]
       endif
    endif
    
    if (n_elements(Zmin) ne 0L) then begin
       w = where(Z ge Zmin,nw)
       if (nw ne 0L) then begin
          Z = Z[w]
          grids = grids[w,*]
       endif
    endif

    q = light*10^U
    nZ = n_elements(Z)
    nU = n_elements(U)

    if keyword_set(log12oh) then begin
       Z = alog10(Z)+Zsun
       Ztitle = '12 + log (O/H)'
       Zrange = [7.3,9.2]
    endif else begin
       Z = alog10(Z)
       Ztitle = 'log Z/Z'+sunsymbol()
       Zrange = [-1.1,0.5]
    endelse
    
; allowed plots    
    
    case plotnumber of

       1L: begin

          x = grids.nii_6584_h_alpha
          y = grids.oiii_5007_h_beta

          xrange = [-2.7,0.2]
          yrange = [-2.0,1.5]
          xtitle = 'log ([N II] \lambda6548 / H\alpha)'
          ytitle = 'log ([O III] \lambda5007 / H\beta)'

          if keyword_set(hiiregions) then begin
             good = where((regions.nii_6584_h_alpha gt -9000.0) and $
               (regions.oiii_5007_h_beta gt -9000.0))
             xregion = regions[good].nii_6584_h_alpha
             yregion = regions[good].oiii_5007_h_beta
          endif
             
       end

       2L: begin

          x = grids.sii_h_alpha
          y = grids.oiii_5007_h_beta

          xrange = [-2.4,0.2]
          yrange = [-2.0,1.5]
          xtitle = 'log ([S II] \lambda\lambda6716,6731 / H\alpha)'
          ytitle = 'log ([O III] \lambda5007 / H\beta)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.sii_h_alpha gt -9000.0) and $
               (regions.oiii_5007_h_beta gt -9000.0))
             xregion = regions[good].sii_h_alpha
             yregion = regions[good].oiii_5007_h_beta
          endif
             
       end

       3L: begin

          x = grids.oi_6300_h_alpha
          y = grids.oiii_5007_h_beta

          xrange = [-3.3,-0.5]
          yrange = [-2.0,1.5]
          xtitle = 'log ([O I] \lambda6300 / H\alpha)'
          ytitle = 'log ([O III] \lambda5007 / H\beta)'
          
          if keyword_set(hiiregions) then begin
             splog, 'No HII region data on [O I].'
;            good = where((regions.oi_6300_h_alpha gt -9000.0) and $
;              (regions.oiii_5007_h_beta gt -9000.0))
;            xregion = regions[good].oi_6300_h_alpha
;            yregion = regions[good].oiii_5007_h_beta
          endif
             
       end

       4L: begin

          x = grids.oiii_5007_oii
          y = grids.oii_h_alpha

          xrange = [-2.0,1.5]
          yrange = [-1.4,0.7]
          xtitle = 'log ([O III] \lambda5007 / [O II] \lambda3727)'
          ytitle = 'log ([O II] \lambda3727 / H\alpha)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.oiii_5007_oii gt -9000.0) and $
               (regions.oii_h_alpha gt -9000.0))
             xregion = regions[good].oiii_5007_oii
             yregion = regions[good].oii_h_alpha
          endif
             
       end

       5L: begin

          x = grids.nii_6584_sii
          y = grids.oii_h_alpha

          xrange = [-1.0,0.9]
          yrange = [-1.4,0.4]
          xtitle = 'log ([N II] \lambda6548 / [S II] \lambda\lambda6716,6731)'
          ytitle = 'log ([O II] \lambda3727 / H\alpha)'
         
          if keyword_set(hiiregions) then begin
             good = where((regions.nii_6584_sii gt -9000.0) and $
               (regions.oii_h_alpha gt -9000.0))
             xregion = regions[good].nii_6584_sii
             yregion = regions[good].oii_h_alpha
          endif
             
       end

       6L: begin

          x = grids.nii_6584_oii
          y = grids.oiii_5007_oii

          xrange = [-2.0,1.0]
          yrange = [-2.0,1.5]
          xtitle = 'log ([N II] \lambda6548 / [O II] \lambda3727)'
          ytitle = 'log ([O III] \lambda5007 / [O II] \lambda3727)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.nii_6584_oii gt -9000.0) and $
               (regions.oiii_5007_oii gt -9000.0))
             xregion = regions[good].nii_6584_oii
             yregion = regions[good].oiii_5007_oii
          endif
             
       end

       7L: begin

          x = grids.nii_6584_oii
          y = grids.oiii_5007_h_beta

          xrange = [-2.0,1.0]
          yrange = [-2.0,1.5]

          xtitle = 'log ([N II] \lambda6548 / [O II] \lambda3727)'
          ytitle = 'log ([O III] \lambda5007 / H\beta)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.nii_6584_oii gt -9000.0) and $
               (regions.oiii_5007_h_beta gt -9000.0))
             xregion = regions[good].nii_6584_oii
             yregion = regions[good].oiii_5007_h_beta
          endif
             
       end

       8L: begin

          x = grids.R23
          y = grids.nii_6584_oii

          xrange = [-0.4,1.2]
          yrange = [-2.0,1.0]
          xtitle = 'log R_{23}'
          ytitle = 'log ([N II] \lambda6548 / [O II] \lambda3727)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.nii_6584_oii gt -9000.0) and (regions.R23 gt -9000.0))
             xregion = regions[good].R23
             yregion = regions[good].nii_6584_oii
          endif
             
       end

       9L: begin

          x = grids.oii_h_beta
          y = grids.oiii_5007_h_beta

          xrange = [-2.0,1.2]
          yrange = [-1.8,1.5]
          xtitle = 'log ([O II] \lambda3727 / H\beta)'
          ytitle = 'log ([O III] \lambda5007 / H\beta)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.oii_h_beta gt -9000.0) and $
               (regions.oiii_5007_h_beta gt -9000.0))
             xregion = regions[good].oii_h_beta
             yregion = regions[good].oiii_5007_h_beta
          endif
             
       end

       10L: begin

          x = grids.R23
          y = grids.oiii_5007_oii

          xrange = [-1.0,1.2]
          yrange = [-2.0,1.5]
          xtitle = 'log R_{23}'
          ytitle = 'log ([O III] \lambda5007 / [O II] \lambda3727)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.R23 gt -9000.0) and $
               (regions.oiii_5007_oii gt -9000.0))
             xregion = regions[good].R23
             yregion = regions[good].oiii_5007_oii
          endif
             
       end

       11L: begin

          x = grids.nii_6584_h_alpha
          y = grids.oii_h_alpha

          xrange = [-2.0,0.1]
          yrange = [-1.4,0.4]
          xtitle = 'log ([N II] \lambda6548 / H\alpha)'
          ytitle = 'log ([O II] \lambda3727 / H\alpha)'
         
       end

       12L: begin

          x = grids.nii_6584_oiii_5007
          y = grids.oii_h_alpha

          xrange = [-3.0,1.5]
          yrange = [-1.4,0.4]
          xtitle = 'log ([N II] \lambda6548 / [O III] \lambda5007)'
          ytitle = 'log ([O II] \lambda3727 / H\alpha)'
         
       end

       13L: begin

          x = grids.oiii_5007_h_beta
          y = grids.oii_h_alpha

          xrange = [-1.8,1.2]
          yrange = [-1.3,0.8]
          xtitle = 'log ([O III] \lambda5007 / H\beta)'
          ytitle = 'log ([O II] \lambda3727 / H\alpha)'
          
       end

       14L: begin

          x = grids.nii_6584_oii
          y = grids.oii_h_alpha

          xrange = [-1.8,0.6]
          yrange = [-1.4,0.4]
          xtitle = 'log ([N II] \lambda6548 / [O II] \lambda3727)'
          ytitle = 'log ([O II] \lambda3727 / H\alpha)'
         
       end

       15L: begin

;         x = transpose(q # (fltarr(nZ)+1))
;         xrange = [4E6,4E8]
;         xtitle = 'q [cm/s]'
;         xlog = 1

          x = transpose(U # (fltarr(nZ)+1))
          xrange = [-1.5,-4.2]
          xtitle = 'log U'
          xlog = 0
          
          y = grids.oii_h_alpha
          yrange = [-1.4,0.4]
          ytitle = 'log ([O II] \lambda3727 / H\alpha)'
         
       end

       16L: begin

;         x = transpose(q # (fltarr(nZ)+1))
;         xrange = [4E6,4E8]
;         xtitle = 'q [cm/s]'
;         xlog = 1

          x = transpose(U # (fltarr(nZ)+1))
          xrange = [-1.5,-4.2]
          xtitle = 'log U'
          xlog = 0
          
          y = grids.oiii_5007_oii
          yrange = [-2.5,2.0]
          ytitle = 'log ([O III] \lambda5007 / [O II] \lambda3727)'

       end

       17L: begin

          x = Z # (fltarr(nU)+1)
          xrange = Zrange
          xtitle = Ztitle

          y = grids.oii_h_alpha
          yrange = [-1.4,0.4]
          ytitle = 'log ([O II] \lambda3727 / H\alpha)'

       end

       18L: begin

          x = transpose(U # (fltarr(nZ)+1))
          xrange = [-1.5,-4.2]
          xtitle = 'log U'
          xlog = 0
          
          y = grids.oiii_5007_oii
          yrange = [-1.4,1.4]
          ytitle = 'log ([O III] \lambda5007 / [O II] \lambda3727)'
         
       end

       19L: begin

          x = Z # (fltarr(nU)+1)
          xrange = Zrange
          xtitle = Ztitle
          
          y = grids.oiii_5007_oii
          yrange = [-1.4,1.4]
          ytitle = 'log ([O III] \lambda5007 / [O II] \lambda3727)'

       end

       20L: begin

          x = grids.nii_6584_oii
          y = grids.oii_h_beta_oii_sii

          xrange = [-2.0,1.0]
          yrange = [-1.0,2.5]
          xtitle = 'log ([N II] \lambda6548 / [O II] \lambda3727)'
          ytitle = 'log ([O II] / H\beta) + log ([O II] / [S II])'

          if keyword_set(hiiregions) then begin
             good = where((regions.nii_6584_oii gt -9000.0) and $
               (regions.oii_h_beta_oii_sii gt -9000.0))
             xregion = regions[good].nii_6584_oii
             yregion = regions[good].oii_h_beta_oii_sii
          endif
             
       end

       21L: begin

          x = grids.R23
          y = grids.oiii_5007_h_beta

          xrange = [-1.0,1.2]
          yrange = [-2.0,1.5]
          xtitle = 'log R_{23}'
          ytitle = 'log ([O III] \lambda5007 / H\beta)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.R23 gt -9000.0) and $
               (regions.oiii_5007_h_beta gt -9000.0))
             xregion = regions[good].R23
             yregion = regions[good].oiii_5007_h_beta
          endif
             
       end

       22L: begin

          x = grids.R23
          y = grids.oii_h_alpha

          xrange = [-1.0,1.2]
          yrange = [-1.4,0.4]
          xtitle = 'log R_{23}'
          ytitle = 'log ([O II] \lambda3727 / H\alpha)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.R23 gt -9000.0) and $
               (regions.oii_h_alpha gt -9000.0))
             xregion = regions[good].R23
             yregion = regions[good].oii_h_alpha
          endif
             
       end

       23L: begin

          x = grids.nii_6584_h_alpha
          y = grids.oiii_5007_oii

          xrange = [-2.9,0.4]
          yrange = [-1.9,1.5]
          xtitle = 'log ([N II] \lambda6548 / H\alpha)'
          ytitle = 'log ([O III] \lambda5007 / [O II] \lambda3727)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.oiii_5007_oii gt -9000.0) and $
               (regions.nii_6584_h_alpha gt -9000.0))
             xregion = regions[good].nii_6584_h_alpha
             yregion = regions[good].oiii_5007_oii
          endif
             
       end

       24L: begin

          x = grids.oiii_5007_oii
          y = grids.oiii_5007_h_beta

          xrange = [-2.0,1.7]
          yrange = [-1.8,1.1]
          xtitle = 'log ([O III] \lambda5007 / [O II] \lambda3727)'
          ytitle = 'log ([O III] \lambda5007 / H\beta)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.oiii_5007_oii gt -9000.0) and $
               (regions.oiii_5007_h_beta gt -9000.0))
             xregion = regions[good].oiii_5007_oii
             yregion = regions[good].oiii_5007_h_beta
          endif
             
       end

       25L: begin

          x = grids.oiii_5007_oii
          y = grids.oi_6300_h_alpha

          xrange = [-1.7,1.7]
          yrange = [-3.0,-0.5]
          xtitle = 'log ([O III] \lambda5007 / [O II] \lambda3727)'
          ytitle = 'log ([O I] \lambda5007 / H\alpha)'
          
       end 

       26L: begin

          x = grids.oi_6300_h_alpha
          y = grids.sii_h_alpha

          xrange = [-3.3,-0.5]
          yrange = [-1.5,0.2]

          xtitle = 'log ([O I] \lambda6300 / H\alpha)'
          ytitle = 'log ([S II] \lambda\lambda6716,6731 / H\alpha)'
          
          if keyword_set(hiiregions) then begin
             splog, 'No HII region data on [O I].'
          endif

       end

       27L: begin

          x = grids.nii_6584_h_alpha
          y = grids.sii_h_alpha

          xrange = [-2.7,0.4]
          yrange = [-1.5,0.2]

          xtitle = 'log ([N II] \lambda6548 / H\alpha)'
          ytitle = 'log ([S II] \lambda\lambda6716,6731 / H\alpha)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.sii_h_alpha gt -9000.0) and $
               (regions.nii_6584_h_alpha gt -9000.0))
             xregion = regions[good].nii_6584_h_alpha
             yregion = regions[good].sii_h_alpha
          endif
             
       end

       28L: begin

          x = grids.O32
          xrange = [-2.4,1.6]
          xtitle = 'log O_{32}'

          y = alog10(Z # (fltarr(nU)+1)) + Zsun
          yrange = [7.4,9.5]
          ytitle = '12 + log (O/H)'
          
       end

       29L: begin ; inverse of plot 9

          x = grids.oiii_5007_h_beta
          y = grids.oii_h_beta

          xrange = [-1.8,1.5]
          yrange = [-2.0,1.2]
          xtitle = 'log ([O III] \lambda5007 / H\beta)'
          ytitle = 'log ([O II] \lambda3727 / H\beta)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.oii_h_beta gt -9000.0) and $
               (regions.oiii_5007_h_beta gt -9000.0))
             xregion = regions[good].oiii_5007_h_beta
             yregion = regions[good].oii_h_beta
          endif
             
       end

       30L: begin

          x = -grids.oiii_5007_oii
          y = grids.nii_6584_h_alpha

          xrange = [-1.5,1.5]
          yrange = [-2,1]
          xtitle = 'log ([O II] \lambda3727 / [O III] \lambda5007)'
          ytitle = 'log ([N II] \lambda6548 / H\alpha)'
          
          if keyword_set(hiiregions) then begin
             good = where((regions.oii_oiii_5007 gt -9000.0) and $
               (regions.nii_6584_h_alpha gt -9000.0))
             xregion = regions[good].oii_oiii_5007
             yregion = regions[good].nii_6584_h_alpha
          endif
             
       end

       else: begin
          splog, 'Unknown plot number!'
          return
       end

    endcase

    if keyword_set(siioffset) then x = x - siiratio ; this is NOT!!! general 
    
    if size(extra,/type) ne 8L then extra = {dummy: ''}

    postthick = 8.0

    if keyword_set(postscript) then thick = postthick else thick = 2.0
    if keyword_set(postscript) then xthick = postthick else xthick = 2.0
    if keyword_set(postscript) then ythick = postthick else ythick = 2.0
    if keyword_set(postscript) then charthick = postthick else charthick = 2.0
    if keyword_set(postscript) then lcharthick = postthick else lcharthick = 2.0

    if tag_exist(extra,'ULINESTYLE') then Ulinestyle = extra.Ulinestyle else Ulinestyle = 2
    if tag_exist(extra,'ZLINESTYLE') then Zlinestyle = extra.Zlinestyle else Zlinestyle = 0
    if tag_exist(extra,'XRANGE') then xrange = extra.xrange
    if tag_exist(extra,'yRANGE') then yrange = extra.yrange
    if tag_exist(extra,'XSTYLE') then xstyle = extra.xstyle else xstyle = 3
    if tag_exist(extra,'YSTYLE') then ystyle = extra.ystyle else ystyle = 3
    if tag_exist(extra,'CHARSIZE') then charsize = extra.charsize else charsize = 2.0
    if tag_exist(extra,'CHARTHICK') then charthick = extra.charthick ; overwrite
    if tag_exist(extra,'THICK') then thick = extra.thick ; overwrite

;   if tag_exist(extra,'THICK') then thick = extra.thick else thick = 2.0
;   if tag_exist(extra,'CHARSIZE') then charsize = extra.charsize else charsize = 2.0
;   if tag_exist(extra,'CHARTHICK') then charthick = extra.charthick else charthick = 2.0

    if not keyword_set(overplot) then begin
       if (not keyword_set(postscript)) and (!d.window ne 0L) then window, 0, xs=500, ys=500
       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
         charsize=charsize, charthick=charthick, xthick=xthick, ythick=xthick, $
         xtitle=xtitle, ytitle=ytitle, xstyle=xstyle, ystyle=ystyle, xlog=xlog, $
         _extra=extra
    endif

    if keyword_set(hiiregions) and (n_elements(xregion) ne 0L) then begin
       plotsym, 8, 0.5, /fill
       djs_oplot, xregion, yregion, ps=8, color='blue'
    endif
    
    if not keyword_set(noUgrid) then for i = 0L, nU-1L do $ ; fixed ionization parameter, varying metallicity
      djs_oplot, x[*,i], y[*,i], thick=thick, color=Ugridcolor, $
        linestyle=Ulinestyle
    
    if not keyword_set(noZgrid) then for i = 0L, nZ-1L do $ ; fixed metallicity, varying ionization parameter
      djs_oplot, x[i,*], y[i,*], thick=thick, color=Zgridcolor, $
        linestyle=Zlinestyle

; label the first and last grid points

    if (n_elements(labeltype) ne 0L) then begin

       case labeltype of

          0L: begin ; no labels
             
          end
          1L: begin ; label all the metallicity and ionization parameter beginning points
             for i = 0L, nZ-1L do $
               if (x[i,0] gt !x.crange[0]) and (x[i,0] lt !x.crange[1]) and $
                 (y[i,0] gt !y.crange[0]) and (y[i,0] lt !y.crange[1]) then $
               xyouts, x[i,0], y[i,0], string(10^z[i],format='(F3.1)'), /data, $
                 align=0.0, charsize=Zcharsize, charthick=charthick, color=djs_icolor(Zlabelcolor)
             for i = 0L, nU-1L do $
               if (x[0,i] gt !x.crange[0]) and (x[0,i] lt !x.crange[1]) and $
                 (y[0,i] gt !y.crange[0]) and (y[0,i] lt !y.crange[1]) then $
               xyouts, x[0,i], y[0,i], string(U[i],format='(F4.1)'), /data, $
                 align=1.0, charsize=Ucharsize, charthick=charthick, color=djs_icolor(Ulabelcolor)
          end
          2L: begin ; label all the metallicity endpoints at the highest ionization parameter
             for i = 0L, nZ-1L do $
               if (x[i,nU-1] gt !x.crange[0]) and (x[i,nU-1] lt !x.crange[1]) and $
                 (y[i,nU-1] gt !y.crange[0]) and (y[i,nU-1] lt !y.crange[1]) then $
               xyouts, x[i,nU-1], y[i,nU-1], string(10^z[i],format='(F3.1)'), /data, $
                 align=0.0, charsize=Zcharsize, charthick=charthick, color=djs_icolor(Zlabelcolor)
          end
          3L: begin ; label all the metallicity endpoints at the lowest ionization parameter
             for i = 0L, nZ-1L do $
               if (x[i,0] gt !x.crange[0]) and (x[i,0] lt !x.crange[1]) and $
                 (y[i,0] gt !y.crange[0]) and (y[i,0] lt !y.crange[1]) then $
               xyouts, x[i,0], y[i,0], string(10^z[i],format='(F3.1)'), /data, $
                 align=1.0, charsize=Zcharsize, charthick=charthick, color=djs_icolor(Zlabelcolor)
          end
          4L: begin ; label all the ionization parameter endpoints at the highest metallicity
             for i = 0L, nU-1L do $
               if (x[nZ-1,i] gt !x.crange[0]) and (x[nZ-1,i] lt !x.crange[1]) and $
                 (y[nZ-1,i] gt !y.crange[0]) and (y[nZ-1,i] lt !y.crange[1]) then $
               xyouts, x[nZ-1,i], y[nZ-1,i], string(U[i],format='(F4.1)'), /data, $
                 align=0.0, charsize=Ucharsize, charthick=charthick, color=djs_icolor(Ulabelcolor)
          end
          5L: begin ; label all the ionization parameter endpoints at the lowest metallicity
             for i = 0L, nU-1L do $
               if (x[0,i] gt !x.crange[0]) and (x[0,i] lt !x.crange[1]) and $
                 (y[0,i] gt !y.crange[0]) and (y[0,i] lt !y.crange[1]) then $
               xyouts, x[0,i], y[0,i], string(U[i],format='(F4.1)'), /data, $
                 align=1.0, charsize=Ucharsize, charthick=charthick, color=djs_icolor(Ulabelcolor)
          end

          else: 
      endcase

    endif

return
end    

