pro Zgto
; jm04jun30uofa
; compute electron-temperature abundances for the GTO starbursts 

    oratio = 2.98D
    ocor = 1.0 + 1.0/oratio

    oiii_indx = 5L     ; index of the [O III] 5007 transition
    oii_3726_indx = 0L ; [O II] 3726
    oii_3728_indx = 1L ; [O II] 3728
    nii_indx = 5L      ; [N II] 6584

; read the data file

    readcol, 'Zgto.dat', galaxy, hiiregion, reference, oii_3727, oii_3727_err, $
      oiii_4363, oiii_4363_err, oiii_5007, oiii_5007_err, $
      nii_6584, nii_6584_err, sii_6716, sii_6716_err, sii_6731, sii_6731_err, $
      log12oh_lit, comment='#', /silent, format='A,A,A,D,D,D,D,D,D,D,D,D,D,D,D,D'
    nobject = n_elements(galaxy)
    
    data = {$
      galaxy:                 '', $
      hiiregion:              '', $
      reference:              '', $
      oii_3727:      [0.0D,0.0D], $
      oiii_4363:     [0.0D,0.0D], $
      oiii_5007:     [0.0D,0.0D], $
      nii_6584:      [0.0D,0.0D], $
      sii_6716:      [0.0D,0.0D], $
      sii_6731:      [0.0D,0.0D], $
      n:                    0.0D, $ ; density
      T_oiii:               0.0D, $
      T_oiii_lower:         0.0D, $
      T_oiii_upper:         0.0D, $
      T_oii:                0.0D, $
      oh_pp:                0.0D, $ ; O++ abundance
      oh_pp_lower:          0.0D, $
      oh_pp_upper:          0.0D, $
      oh_p:                 0.0D, $ ; O+ abundance
      oh_p_lower:           0.0D, $
      oh_p_upper:           0.0D, $
      log12oh:              0.0D, $ ; total oxygen abundance
      log12oh_lower:        0.0D, $
      log12oh_upper:        0.0D, $
      log12oh_lit:          0.0D, $ ; literature oxygen abundance
      log_no:            -999.0D}   ; log N/O ratio
    data = replicate(data,nobject)

    data.galaxy    = galaxy
    data.hiiregion = hiiregion
    data.reference = reference
    data.log12oh_lit = log12oh_lit

    for iobj = 0L, nobject-1L do begin

       splog, 'Considering galaxy '+data[iobj].galaxy+' - '+data[iobj].hiiregion+'.'

       data[iobj].oii_3727  = [oii_3727[iobj],oii_3727_err[iobj]]
       data[iobj].oiii_4363 = [oiii_4363[iobj],oiii_4363_err[iobj]]
       data[iobj].oiii_5007 = [oiii_5007[iobj],oiii_5007_err[iobj]]
       data[iobj].nii_6584  = [nii_6584[iobj],nii_6584_err[iobj]]
       data[iobj].sii_6716  = [sii_6716[iobj],sii_6716_err[iobj]]
       data[iobj].sii_6731  = [sii_6731[iobj],sii_6731_err[iobj]]

       oiii_ratio     = data[iobj].oiii_5007[0]*ocor/data[iobj].oiii_4363[0]
       oiii_ratio_err = im_compute_error(data[iobj].oiii_5007[0]*ocor,data[iobj].oiii_5007[1],$
         data[iobj].oiii_4363[0],data[iobj].oiii_4363[1],/quotient)

       sii_ratio     = data[iobj].sii_6716[0]/data[iobj].sii_6731[0]
       sii_ratio_err = im_compute_error(data[iobj].sii_6716[0],data[iobj].sii_6716[1],$
         data[iobj].sii_6731[0],data[iobj].sii_6731[1],/quotient)

; compute the temperature assuming a constant density, then compute
; the density given the temperature, and finally re-compute the
; temperature with Monte Carlo errors    
       
       fivel = idl_fivel(1,6,lineratio=oiii_ratio)
       fivel = idl_fivel(1,11,lineratio=sii_ratio,temperature=fivel.temperature)
       fivel = idl_fivel(1,6,temperature=fivel.temperature,density=fivel.density,$
         lineratio=oiii_ratio,err_lineratio=oiii_ratio_err,/montecarlo)

       data[iobj].t_oiii = fivel.temperature
       data[iobj].t_oiii_lower = fivel.temperature_lower
       data[iobj].t_oiii_upper = fivel.temperature_upper

       data[iobj].n = fivel.density

; adopt the Garnett (1992) photoionization relations to predict the
; temperature of the O+ zone from the O++ zone

       data[iobj].t_oii = 0.70*data[iobj].t_oiii + 3000.0
;      data[iobj].t_oii = 2D4/(1D4/data[iobj].t_oiii + 0.8) ; Stasinska (1990)

; compute the emissivities of the [O III] and [O II] lines; finally
; compute the partial-ionization zone abundances and the total
; abundance 

       fivel = idl_fivel(2,7,temperature=data[iobj].t_oiii,density=data[iobj].n)
       data[iobj].oh_pp = data[iobj].oiii_5007[0] * fivel.jhb / fivel.emissivity[oiii_indx]
       
       fivel = idl_fivel(2,6,temperature=data[iobj].t_oii,density=data[iobj].n)
       oii_3728 = data[iobj].oii_3727[0] / (1.0 + fivel.emissivity[oii_3726_indx]/fivel.emissivity[oii_3728_indx])
       data[iobj].oh_p = oii_3728 * fivel.jhb / fivel.emissivity[oii_3728_indx]

       data[iobj].log12oh = alog10(data[iobj].oh_pp + data[iobj].oh_p) + 12.0

; finally compute the log N/O ratio

       if (data[iobj].nii_6584[1] gt 0.0) then begin
          
          fivel = idl_fivel(2,3,temperature=data[iobj].t_oii,density=data[iobj].n)
          data[iobj].log_no = alog10((data[iobj].nii_6584[0] * fivel.jhb / fivel.emissivity[nii_indx]) / data[iobj].oh_p)

       endif
          
       print

    endfor

    out = struct_trimtags(data,select=['GALAXY','HIIREGION','REFERENCE','N','T_OIII',$
      'T_OII','LOG12OH','LOG12OH_LIT','LOG_NO'])
    struct_print, out

stop
    
return
end
    
