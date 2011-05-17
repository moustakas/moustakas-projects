function mzlz_grab_info, stroh, stranc, strmass, flux=flux, $
  t04=t04, m91=m91, kk04=kk04, mpajhu=mpajhu, nolimits=nolimits, $
  ambig=ambig, zmin=zmin, zmax=zmax, errcut=errcut
; jm09mar19nyu - grab the specified abundances and global properties

    if (n_elements(zmin) eq 0) then zmin = min(stranc.z)
    if (n_elements(zmax) eq 0) then zmax = max(stranc.z)

    id = lindgen(n_elements(stroh))
    masstag = tag_indx(strmass,'mass_50')
;   masstag = tag_indx(strmass,'mass_avg')
    masserrtag = tag_indx(strmass,'mass_err')
    
; special case for the SDSS - grab the MPA-JHU metallicities
    if keyword_set(mpajhu) then begin
       indx = where((stranc.oh_avg gt 0.0),nindx)
       oh = stranc[indx].oh_avg
       oh_err = (stranc[indx].oh_p84-stranc[indx].oh_p16)/2.0
       ohlimit = oh*0
    endif else begin

; grab T04 by default    
       if (keyword_set(t04) eq 0) and (keyword_set(m91) eq 0) and $
         (keyword_set(kk01) eq 0) then t04 = 1
       if keyword_set(t04) then calib = 't04'
       if keyword_set(m91) then calib = 'm91'
       if keyword_set(kk04) then calib = 'kk04'

; grab EWs by default, unless /FLUX
       if keyword_set(flux) then begin
          prefix = '' 
          suffix = '' 
       endif else begin
          prefix = 'ew_'
          suffix = '_ew'
       endelse

; optionally return objects with ambiguous abundances
       if keyword_set(ambig) then branch = 'A' else branch = 'U'
       
       branchtag = tag_indx(stroh,'r23branch'+'_'+calib+suffix)
       ohtag = tag_indx(stroh,prefix+'log12oh_'+calib)
       oherrtag = tag_indx(stroh,prefix+'log12oh_'+calib+'_err')

       if keyword_set(nolimits) then begin ; no upper limits
          indx = where((stroh.(ohtag) gt -900) and $
            (stroh.ohlimit eq 0) and $
            (stroh.(oherrtag) gt -900) and $
            (strtrim(stroh.(branchtag),2) eq branch) and $
            (stranc.z ge zmin) and (stranc.z le zmax),nindx)
       endif else begin         ; allow limits
          indx = where((stroh.(ohtag) gt -900) and $
            (stroh.(oherrtag) gt -900) and $
            (strtrim(stroh.(branchtag),2) eq branch) and $
            (stranc.z ge zmin) and (stranc.z le zmax),nindx)
       endelse

       if keyword_set(errcut) then begin
          ecut = weighted_quantile(stroh[indx].(oherrtag),quant=0.9)
          indx = indx[where(stroh[indx].(oherrtag) lt ecut)]
          nindx = n_elements(indx)
       endif
       
       if (nindx eq 0L) then begin
          splog, 'No objects found!'
          return, -1
       endif
       oh = stroh[indx].(ohtag)
       oh_err = stroh[indx].(oherrtag)
       ohlimit = stroh[indx].ohlimit
    endelse

; pack into a structure and return
    if tag_exist(stranc,'FINAL_WEIGHT') then $
      weight = stranc[indx].final_weight else $
        weight = oh*0.0+1.0

    if tag_exist(stranc,'OH_AVG') then $
      tremonti_oh = stranc[indx].oh_avg else $
        tremonti_oh = oh*0.0
    if tag_exist(stranc,'MASS_MEDIAN') then $
      kauffmann_mass = stranc[indx].mass_median else $
        kauffmann_mass = oh*0.0
       
;    ubvri_v2ab = k_vega2ab(filterlist=bessell_filterlist(),/silent,/kurucz)
;    zobj = strmass[indx].z

    out = {$
      id:      id[indx],$
      indx:        indx,$
      ngal:       nindx,$
      zmin:        zmin,$
      zmax:        zmax,$
      z:      stranc[indx].z,$
      oh:            oh,$
      oh_err:    oh_err,$
      weight:    weight,$
;     mb_vega:   strmass[indx].k_ubvrijhk_absmag_00[1]-ubvri_v2ab[1],$ ; Vega
;     mv_vega:   strmass[indx].k_ubvrijhk_absmag_00[2]-ubvri_v2ab[2],$ ; Vega
;     mr_vega:   strmass[indx].k_ubvrijhk_absmag_00[3]-ubvri_v2ab[3],$ ; Vega
;     mu_ab:     strmass[indx].k_ubvrijhk_absmag_00[0],$ ; AB
      mb_ab:     stranc[indx].k_ubvrijhk_absmag_00[1],$ ; AB
;     mg_ab:     strmass[indx].k_ugriz_absmag_01[1],$ ; AB
;     mr_ab:     strmass[indx].k_ugriz_absmag_01[2],$ ; AB
;     mi_ab:     strmass[indx].k_ugriz_absmag_01[3],$ ; AB
;     ub_ab:     fltarr(nindx),$
;     gr_ab:     fltarr(nindx),$
      mass:      strmass[indx].(masstag),$
      mass_err:  strmass[indx].(masserrtag),$
      ohlimit:   ohlimit,$
      tremonti_oh: tremonti_oh,$
      kauffmann_mass: kauffmann_mass}

; some colors    
;   out.ub_ab = out.mu_ab-out.mb_ab
;   out.gr_ab = out.mg_ab-out.mr_ab
    
return, out
end

