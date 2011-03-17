pro hcn_bpt
; jm09jul28nyu - customized version of Fig 1, the BPT diagram

    common hcn_bpt, sdss
    
; read the samples    
    
    path = hcn_path()

    sfsym = 108 & sfcolor = 'red' & sfpsize = 1.5 & sffill = 1
    sfsym_ul = 114 & sfcolor_ul = 'red' & sfpsize_ul = 4.0 & sffill_ul = 1

    agnsym = 105 & agncolor = 'blue' & agnpsize = 1.8 & agnfill = 1
    agnsym_ul = 114 & agncolor_ul = 'blue' & agnpsize_ul = 4.0 & agnfill_ul = 1

    agnsfsym = 106 & agnsfcolor = 'green' & agnsfpsize = 1.8 & agnsffill = 1
    agnsfsym_ul = 114 & agnsfcolor_ul = 'green' & agnsfpsize_ul = 4.0 & agnsffill_ul = 1

    sdsssym = 108 & sdsscolor = 'grey' & sdsspsize = 0.1 & sdssfill = 1

; ---------------------------------------------------------------------------    
; SDSS/BPT diagram

    im_plotconfig, 0, pos, psfile=path+'hcn_sdss_bpt.eps', xmargin=[1.5,0.4]

; ATLAS/HCN
    atlas_hcn = mrdfits(path+'atlas_new_sample_090612.fits.gz',1)

    indx = where((atlas_hcn.oiii_hb gt -900.0) and (atlas_hcn.nii_ha gt -900.0) and $
      (strmatch(atlas_hcn.galaxy,'*ngc7469*',/fold) eq 0),nindx)
    x = atlas_hcn[indx].nii_ha
    y = atlas_hcn[indx].oiii_hb
    xerr = atlas_hcn[indx].nii_ha_err
    yerr = atlas_hcn[indx].oiii_hb_err
    
    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

; SDSS
    if (n_elements(sdss) eq 0) then begin
       sdss1 = read_sdss_vagc_mpa(/ispecline)
    
       snrcut = 5.0
       indx = where($
         (sdss1.nii_6584[0]/sdss1.nii_6584[1] gt snrcut) and $
         (sdss1.oiii_5007[0]/sdss1.oiii_5007[1] gt snrcut) and $
         (sdss1.h_beta[0]/sdss1.h_beta[1] gt snrcut) and $
         (sdss1.h_alpha[0]/sdss1.h_alpha[1] gt snrcut),nsdss)
       sdss = sdss1[indx]
    endif

    xsdss = alog10(sdss.nii_6584[0]/sdss.h_alpha[0])
    ysdss = alog10(sdss.oiii_5007[0]/sdss.h_beta[0])

; make the plot    
    
    xtitle = 'log ([N II] \lambda6584/H\alpha)'
    ytitle = 'log ([O III] \lambda5007/H\beta)'

    xrange = [-1.5,0.7]
    yrange = [-1.1,1.2]

;   im_symbols, sdsssym, thick=postthick1, psize=sdsspsize, fill=sdssfill, color=sdsscolor
    hogg_scatterplot, xsdss, ysdss, /outliers, label=0, outpsym=3, outsymsize=1.0, $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=2.0, position=pos[*,0], outcolor=djs_icolor(sdsscolor), cthick=2.0
    im_symbols, sfsym, psize=sfpsize*0.8, fill=sffill, color=sfcolor
    oploterror, x[sf], y[sf], xerr[sf], yerr[sf], psym=8, $
      errcolor=djs_icolor(sfcolor)
    im_symbols, agnsym, psize=agnpsize, fill=agnfill, color=agncolor
    oploterror, x[agn], y[agn], xerr[agn], yerr[agn], psym=8, $
      errcolor=djs_icolor(agncolor)
    im_symbols, agnsfsym, psize=agnsfpsize*0.7, fill=agnsffill, color=agnsfcolor
    oploterror, x[agnsf], y[agnsf], xerr[agnsf], yerr[agnsf], psym=8, $
      errcolor=djs_icolor(agnsfcolor)
    
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)
    sjgalaxy = strcompress(atlas_hcn[indx].sj_galaxy,/remove)
    for ii = 0L, nindx-1L do begin
       align = 0.0 & xoff = 0.05 & yoff = 0.0
       if strmatch(galaxy[ii],'*ngc6240*',/fold) then begin
          align = 0.0 & xoff = +0.06 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*ic1623*',/fold) then begin
          align = 1.0 & xoff = -0.04 & yoff = 0.02
       endif
       if strmatch(galaxy[ii],'*ugc05101*',/fold) then begin
          align = 1.0 & xoff = -0.06 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*iras05189*',/fold) then begin
          align = 1.0 & xoff = -0.06 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*ugc08696*',/fold) then begin
          align = 0.0 & xoff = +0.06 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*iras17208*',/fold) then begin
          align = 1.0 & xoff = -0.07 & yoff = -0.0
       endif
       if strmatch(galaxy[ii],'*ngc3628*',/fold) then begin
          align = 1.0 & xoff = -0.04 & yoff = -0.05
       endif
       if strmatch(galaxy[ii],'*ngc0660*',/fold) then begin
          align = 0.0 & xoff = +0.03 & yoff = -0.04
       endif
       if strmatch(galaxy[ii],'*ngc3079*',/fold) then begin
          align = 1.0 & xoff = -0.03 & yoff = +0.02
       endif
       if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin
          align = 1.0 & xoff = -0.02 & yoff = -0.07
       endif
       if strmatch(galaxy[ii],'*ngc1144*',/fold) then begin
          align = 1.0 & xoff = -0.06 & yoff = 0.0
       endif
       if strmatch(galaxy[ii],'*ngc3893*',/fold) then begin
          align = 1.0 & xoff = -0.05 & yoff = 0.0
       endif
       if strmatch(galaxy[ii],'*ngc0695*',/fold) then begin
          align = 0.0 & xoff = 0.01 & yoff = 0.02
       endif
       if strmatch(galaxy[ii],'*ngc2903*',/fold) then begin
          align = 1.0 & xoff = -0.05 & yoff = -0.1
       endif
       if strmatch(galaxy[ii],'*ugc04881*',/fold) then begin
          align = 0.0 & xoff = +0.05 & yoff = -0.06
       endif
       if strmatch(galaxy[ii],'*ic5179*',/fold) then begin
          align = 1.0 & xoff = -0.05 & yoff = -0.05
       endif
       if strmatch(galaxy[ii],'*ngc7130*',/fold) then begin
          align = 0.5 & xoff = +0.0 & yoff = +0.05
       endif
       if strmatch(galaxy[ii],'*ngc7771*',/fold) then begin
          align = 0.0 & xoff = +0.04 & yoff = +0.01
       endif
       if strmatch(galaxy[ii],'*ic0883*',/fold) then begin
          align = 0.0 & xoff = +0.06 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*mrk0331*',/fold) then begin
          align = 0.0 & xoff = +0.06 & yoff = -0.03
       endif
       if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
          align = 0.0 & xoff = +0.03 & yoff = 0.0
       endif
       if strmatch(galaxy[ii],'*ngc6701*',/fold) then begin
          align = 0.0 & xoff = +0.05 & yoff = +0.01
       endif
       if strmatch(galaxy[ii],'*ngc4414*',/fold) then begin
          align = 1.0 & xoff = +0.0 & yoff = +0.04
       endif
       if strmatch(galaxy[ii],'*ngc0520*',/fold) then begin
          align = 0.0 & xoff = +0.05 & yoff = -0.03
       endif
;      if strmatch(galaxy[ii],'*ngc7469*',/fold) then begin
;         align = 0.0 & xoff = +0.06 & yoff = -0.02
;      endif
       if strmatch(galaxy[ii],'*ngc2146*',/fold) then begin
          align = 1.0 & xoff = -0.02 & yoff = +0.02
       endif
       if strmatch(galaxy[ii],'*ngc3893*',/fold) then begin
          align = 1.0 & xoff = -0.03 & yoff = -0.06
       endif
       if strmatch(galaxy[ii],'*ngc5194*',/fold) then begin
          align = 0.0 & xoff = +0.02 & yoff = -0.09
       endif
       xyouts, x[ii]+xoff, y[ii]+yoff, sjgalaxy[ii], $
         charsize=1.0, /data, align=align, charthick=2.5
    endfor

; overplot the mixing lines

    models = kewley_bpt_lines(/kauffmann,_extra=extra)
    oplot, models.x_nii, models.y_nii, line=0
    models = kewley_bpt_lines(_extra=extra)
    oplot, models.x_nii, models.y_nii, line=2

    im_plotconfig, /psclose

stop    
    
return
end
    
