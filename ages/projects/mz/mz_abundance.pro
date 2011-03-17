;+
; NAME:
;   MZ_ABUNDANCE()
;
; PURPOSE:
;   Compute empirical abundances and various relevant line-ratios
;   from strong emission lines. 
;
; INPUTS:
;   line   - input structure of de-reddened line fluxes
;
; OPTIONAL INPUTS:
;   snrcut_abundance - require that individual emission lines have
;                      a signal-to-noise ratio greater than
;                      SNRCUT_ABUNDANCE 
;   ewalpha          - alpha parameter for Kobulnicky & Phillips
;                      (2003) EW abundances
;   nmonte           - number of Monte Carlo realizations; note
;                      that if NMONTE=0 then the abundance
;                      uncertainties for some calibrations will
;                      not be computed 
;
; KEYWORD PARAMETERS:
;   use_4959  - use the measured 4959 emission-line flux; the
;               default is to assume 4959 = 5007/2.984; note that
;               if LINE has not been reddening-corrected, then the
;               assumption of the intrinsic ratio will be wrong by
;               a little bit 
;   getdensity - compute electron densities (slow!)
;   silent    - do not print messages to STDOUT
;
; OUTPUTS:
;   abund - output data structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2004 Feb 8, U of A, written
;   jm04jul06uofa - added empirical Te-based abundances 
;   jm04jul22uofa - added SILENT and COMBINATION keywords 
;   jm04jul27uofa - added a warning flag if R23>1.1; added
;                   ELECTRONDENSITY keyword
;   jm04nov03uofa - documentation updated; HBETA keyword replaced
;                   with HbHg keyword for consistency with
;                   IUNRED_LINEDUST() 
;   jm04dec08uofa - general development of various empirical
;                   diagnostics; remove objects with log R23 > 1 
;   jm05jan01uofa - compute EW line-ratios and abundances
;   jm05jan06uofa - do not require reddening-corrected line
;                   fluxes
;   jm05feb08uofa - added NOTSTRICT keyword
;   jm05mar31uofa - cleaned up the structure field names
;   jm05sep05uofa - removed NOTSTRICT keyword in favor of
;                   R23STRICT keyword
;   jm06mar10uofa - added EWALPHA and NMONTE parameters 
;   jm07dec10nyu  - added S23
;   jm08feb06nyu  - somewhat major overhaul; cleaned up and
;                   streamlined; removed NORMALIZE and
;                   ELECTRONDENSITY keywords; added JUSTEW and
;                   JUSTFLUX keywords; KK04, M91, and PT05
;                   abundances now computed formally using their
;                   own routines
;   jm08mar20nyu - added NODENSITY and USE_4959 keywords; allow
;                  NMONTE=0 
;   jm08aug21nyu - compute the T04 calibration using EWs
;   jm08oct01nyu - added JUSTRATIOS keyword
;   jm10jul23ucsd - NODENSITY keyword changed to GETDENSITY; CL01 and
;     KD02/combined abundances removed since I never use them
;
; Copyright (C) 2004-2008, 2010, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function get_r23lines, line, ew=ew
; support routine for when we call, e.g., MONTE_LOG12OH_KK04
    r23lines = replicate({oii: [-999.0,-999.0], oiii: [-999.0,-999.0], $
      hbeta: [-999.0,-999.0]},n_elements(line))
    if keyword_set(ew) then begin
       r23lines.oii = line.oii_3727_ew
       r23lines.oiii = line.oiii_4959_5007_ew
       r23lines.hbeta = line.h_beta_ew
    endif else begin
       r23lines.oii = line.oii_3727
       r23lines.oiii = line.oiii_4959_5007
       r23lines.hbeta = line.h_beta
    endelse
return, r23lines
end    

function mz_abundance, line, snrcut_abundance=snrcut_abundance, $
  ewalpha=ewalpha1, justew=justew, justflux=justflux, justratios=justratios, $
  nmonte=nmonte, getdensity=getdensity, use_4959=use_4959, silent=silent, $
  debug=debug
    
    nspec = n_elements(line)
    if (nspec eq 0L) then begin
       doc_library, 'mz_abundance'
       return, -1L
    endif

    oratio = 2.984 ; intrinsic 5007/4959 ratio
    ocor = 1.0 + 1.0/oratio
    
; default: compute abundances based on both EWS and FLUXES 
    
    if (not keyword_set(justew)) and (not keyword_set(justflux)) then begin
       justew = 1
       justflux = 1
    endif
    
    if (n_elements(nmonte) eq 0) then nmonte = 250
    if (n_elements(snrcut_abundance) eq 0) then snrcut_abundance = 5.0
    if (keyword_set(silent) eq 0) then splog, 'S/N > '+$
      string(snrcut_abundance,format='(G0.0)')+', '+$
      'NMONTE = '+string(nmonte,format='(G0.0)')

; initialize the output data structure    
    
    abund = {$
; flux-ratios of interest
      zstrong_r23:                       -999.0, $ ; ([O II] + [O III] 4959,5007) / H-beta
      zstrong_r23_err:                   -999.0, $
      zstrong_o32:                       -999.0, $ ; [O III] 4959,5007 / [O II] 3727
      zstrong_o32_err:                   -999.0, $
      zstrong_p:                         -999.0, $ ; [O III] 4959,5007 / ([O II] + [O III] 4959,5007) 
      zstrong_p_err:                     -999.0, $
;      zstrong_oiiinii:                   -999.0, $ ; ([O III] 5007 / H-beta) / ([N II] 6584 / H-alpha)
;      zstrong_oiiinii_err:               -999.0, $
      zstrong_niiha:                     -999.0, $ ; [N II] 6584 / H-alpha
      zstrong_niiha_err:                 -999.0, $
      zstrong_niioii:                    -999.0, $ ; [N II] 6584 / [O II] 3727
      zstrong_niioii_err:                -999.0, $
;      zstrong_niisii:                    -999.0, $ ; [N II] 6584 / [S II] 6716,6731
;      zstrong_niisii_err:                -999.0, $
;      zstrong_oiioiii:                   -999.0, $ ; [O II] 3727 / [O III] 5007
;      zstrong_oiioiii_err:               -999.0, $
;      zstrong_oiiihb:                    -999.0, $
;      zstrong_oiiihb_err:                -999.0, $

;; density-sensitive ratios and densities
;      zstrong_sii:                       -999.0, $ ; [S II] 6716 / [S II] 6731
;      zstrong_sii_err:                   -999.0, $ 
;      zstrong_oii:                       -999.0, $ ; [O II] 3726 / [O II] 3729
;      zstrong_oii_err:                   -999.0, $ 
;      zstrong_sii_dens:                  -999.0, $ ; [S II] electron density (default 100 cm-3)
;      zstrong_sii_dens_err:              -999.0, $
;      zstrong_oii_dens:                  -999.0, $ ; [O II] electron density (default 100 cm-3)
;      zstrong_oii_dens_err:              -999.0, $

; R23-based O/H calibrations
;     zstrong_12oh_m91_frac:             -999.0, $ ; McGaugh (1991)
      zstrong_12oh_m91_avg:              -999.0, $
      zstrong_12oh_m91_avg_err:          -999.0, $
      zstrong_12oh_m91_upper:            -999.0, $
      zstrong_12oh_m91_upper_err:        -999.0, $
      zstrong_12oh_m91_lower:            -999.0, $
      zstrong_12oh_m91_lower_err:        -999.0, $

;     zstrong_12oh_kk04_frac:            -999.0, $ ; Kobulnicky & Kewley (2004)
      zstrong_12oh_kk04_avg:             -999.0, $
      zstrong_12oh_kk04_avg_err:         -999.0, $
      zstrong_12oh_kk04_upper:           -999.0, $
      zstrong_12oh_kk04_upper_err:       -999.0, $
      zstrong_12oh_kk04_lower:           -999.0, $
      zstrong_12oh_kk04_lower_err:       -999.0, $
      zstrong_logu_kk04_avg:             -999.0, $
      zstrong_logu_kk04_avg_err:         -999.0, $
      zstrong_logu_kk04_upper:           -999.0, $
      zstrong_logu_kk04_upper_err:       -999.0, $
      zstrong_logu_kk04_lower:           -999.0, $
      zstrong_logu_kk04_lower_err:       -999.0, $
      zstrong_converge_kk04_lower:            0, $
      zstrong_converge_kk04_upper:            0, $
                                         
;      zstrong_12oh_zkh94:                -999.0, $ ; ZKH94 (upper branch only)
;      zstrong_12oh_zkh94_err:            -999.0, $
      zstrong_12oh_t04:                  -999.0, $ ; Tremonti et al. 2004 (upper branch only)
      zstrong_12oh_t04_err:              -999.0, $

;; other O/H calibrations
;      zstrong_12oh_niiha_pettini:        -999.0, $ ; Pettini & Pagel (2004) [N II]/Ha calibration
;      zstrong_12oh_niiha_pettini_err:    -999.0, $ 
;      zstrong_12oh_oiiinii_pettini:      -999.0, $ ; empirical abundance based on the Pettini & Pagel (2004) calibration
;      zstrong_12oh_oiiinii_pettini_err:  -999.0, $ 

; EW-ratios of interest
      ewalpha:                              1.0, $
      zstrong_ew_r23:                    -999.0, $ ; {ew([O II]) + ew([O III] 4959,5007)} / ew(H-beta)
      zstrong_ew_r23_err:                -999.0, $
      zstrong_ew_o32:                    -999.0, $ ; ew([O III] 4959,5007) / ew([O II] 3727)
      zstrong_ew_o32_err:                -999.0, $
      zstrong_ew_p:                      -999.0, $ ; [O III] 4959,5007 / ([O II] + [O III] 4959,5007) 
      zstrong_ew_p_err:                  -999.0, $

; EW(R23)-based O/H calibrations
;     zstrong_ew_12oh_m91_frac:          -999.0, $ ; McGaugh (1991)
      zstrong_ew_12oh_m91_avg:           -999.0, $
      zstrong_ew_12oh_m91_avg_err:       -999.0, $
      zstrong_ew_12oh_m91_upper:         -999.0, $
      zstrong_ew_12oh_m91_upper_err:     -999.0, $
      zstrong_ew_12oh_m91_lower:         -999.0, $
      zstrong_ew_12oh_m91_lower_err:     -999.0, $
                                         
;     zstrong_ew_12oh_kk04_frac:         -999.0, $ ; Kobulnicky & Kewley (2004)
      zstrong_ew_12oh_kk04_avg:          -999.0, $
      zstrong_ew_12oh_kk04_avg_err:      -999.0, $
      zstrong_ew_12oh_kk04_upper:        -999.0, $
      zstrong_ew_12oh_kk04_upper_err:    -999.0, $
      zstrong_ew_12oh_kk04_lower:        -999.0, $
      zstrong_ew_12oh_kk04_lower_err:    -999.0, $
      zstrong_ew_logu_kk04_avg:          -999.0, $
      zstrong_ew_logu_kk04_avg_err:      -999.0, $
      zstrong_ew_logu_kk04_upper:        -999.0, $
      zstrong_ew_logu_kk04_upper_err:    -999.0, $
      zstrong_ew_logu_kk04_lower:        -999.0, $
      zstrong_ew_logu_kk04_lower_err:    -999.0, $
      zstrong_ew_converge_kk04_lower:         0, $
      zstrong_ew_converge_kk04_upper:         0, $
                                         
      zstrong_ew_12oh_t04:               -999.0, $ ; Tremonti et al. 2004 (upper branch only)
      zstrong_ew_12oh_t04_err:           -999.0}

    abund = replicate(abund,nspec)

; ###########################################################################
; compute EW-ratios and abundances
; ###########################################################################

    if keyword_set(justew) then begin

       if (not keyword_set(silent)) then begin
          if keyword_set(justratios) then splog, 'Computing EW line-ratios' else $
            splog, 'Computing line-ratios and abundances from EWs:'
       endif
       
; choose EWALPHA and then "correct" EW([O II]) to account for it
       
       if (n_elements(ewalpha1) eq 0L) then begin
          splog, 'Using alpha = 1.0 for all galaxies'
          ewalpha = replicate(1.0,nspec) 
       endif else begin
          if (n_elements(ewalpha1) eq 1L) then begin
             splog, 'Using alpha = '+string(ewalpha1,format='(G0.0)')+' for all galaxies'
             ewalpha = replicate(ewalpha1,nspec) 
          endif else begin
             if (nspec ne n_elements(ewalpha1)) then begin
                splog, 'Dimensions of LINE and EWALPHA do not agree'
                return, abund
             endif else begin
                splog, 'Using individual alpha value for each galaxy'
                ewalpha = ewalpha1
             endelse
          endelse
       endelse
       abund.ewalpha = ewalpha
       
       ewline = line
       if (tag_exist(ewline[0],'OII_3727_EW')) then begin
          doit = where((ewline.oii_3727_ew[1] gt 0.0),ndoit)
          if (ndoit ne 0L) then begin
             ewline[doit].oii_3727_ew[0] = ewalpha*ewline[doit].oii_3727_ew[0]
             ewline[doit].oii_3727_ew[1] = ewalpha*ewline[doit].oii_3727_ew[1]
          endif
       endif

; if /USE_4959 then use the measured [O III] 4959 line, otherwise use
; the intrinsic ratio, 2.984; note that if LINE has not been
; reddening-corrected, then the assumption of the intrinsic ratio will
; be wrong by a little bit (although it doesn't really matter for EWs)
       ewline = struct_addtags(temporary(ewline),replicate({oiii_4959_5007_ew: [0.0,-2.0]},nspec))
       if keyword_set(use_4959) then begin
          if tag_exist(ewline[0],'OIII_4959_EW') and tag_exist(ewline[0],'OIII_5007_EW') then begin
             doit = where((ewline.oiii_4959_ew[1] gt 0.0) and (ewline.oiii_5007_ew[1] gt 0.0),ndoit)
             if (ndoit ne 0L) then begin
                ewline[doit].oiii_4959_5007_ew[0] = ewline[doit].oiii_4959_ew[0]+ewline[doit].oiii_5007_ew[0]
                ewline[doit].oiii_4959_5007_ew[1] = sqrt(ewline[doit].oiii_4959_ew[1]^2+ewline[doit].oiii_5007_ew[1]^2)
             endif
          endif 
       endif else begin
          if (tag_exist(ewline[0],'OIII_5007_EW')) then begin
             doit = where((ewline.oiii_5007_ew[1] gt 0.0),ndoit)
             if (ndoit ne 0L) then ewline[doit].oiii_4959_5007_ew = ocor*ewline[doit].oiii_5007_ew
          endif 
       endelse
       
; EW(R23)
       if (tag_exist(ewline[0],'OII_3727_EW') and tag_exist(ewline[0],'OIII_4959_5007_EW') and $
         tag_exist(ewline[0],'H_BETA_EW')) then begin
          lineratio, ewline, ['OII_3727_EW','OIII_4959_5007_EW'], 'H_BETA_EW', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_ew_r23 = x
             abund[index].zstrong_ew_r23_err = xerr
          endif
       endif

; EW(O32)
       if (tag_exist(ewline[0],'OII_3727_EW') and tag_exist(ewline[0],'OIII_4959_5007_EW')) then begin
          lineratio, ewline, 'OIII_4959_5007_EW', 'OII_3727_EW', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_ew_o32 = x
             abund[index].zstrong_ew_o32_err = xerr
          endif
       endif

; EW(P)-parameter (Pilyugin 2001)
       if tag_exist(ewline[0],'OIII_4959_5007_EW') and tag_exist(ewline[0],'OII_3727_EW') then begin
          lineratio, ewline, 'OIII_4959_5007_EW', ['OII_3727_EW','OIII_4959_5007_EW'], $
            '', '', x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_ew_p = x
             abund[index].zstrong_ew_p_err = xerr
          endif
       endif

; save memory!
       ewr23lines = get_r23lines(ewline,/ew)
       ewline = 0

       if (keyword_set(justratios) eq 0) then begin ; compute abundances
          
; ---------------------------------------------------------------------------    
; McGaugh (1991) as fitted in Kobulnicky et al. (1999) - based on EWs
          if (not keyword_set(silent)) then splog, '   M91...'
          index = where((abund.zstrong_ew_o32 gt -900.0) and (abund.zstrong_ew_r23 gt -900.0),nindex)
          if (nindex ne 0L) then begin

             result = monte_log12oh_m91(abund[index].zstrong_ew_r23,$
               abund[index].zstrong_ew_o32,line=ewr23lines[index],$
               nmonte=nmonte)

;            abund[index].zstrong_ew_12oh_m91_frac = result.frac
             abund[index].zstrong_ew_12oh_m91_avg    = result.log12oh_avg
             abund[index].zstrong_ew_12oh_m91_upper  = result.log12oh_upper
             abund[index].zstrong_ew_12oh_m91_lower  = result.log12oh_lower

             abund[index].zstrong_ew_12oh_m91_avg_err    = result.log12oh_avg_err
             abund[index].zstrong_ew_12oh_m91_upper_err  = result.log12oh_upper_err
             abund[index].zstrong_ew_12oh_m91_lower_err  = result.log12oh_lower_err

          endif

; ---------------------------------------------------------------------------    
; Kobulnicky & Kewley (2004) - EW/R23
          if (not keyword_set(silent)) then splog, '   KK04...'
          index = where((abund.zstrong_ew_o32 gt -900.0) and (abund.zstrong_ew_r23 gt -900.0),nindex)
          if (nindex ne 0L) then begin

             result = monte_log12oh_kk04(abund[index].zstrong_ew_r23,$
               abund[index].zstrong_ew_o32,line=ewr23lines[index],$
               nmonte=nmonte,debug=debug)

;            abund[index].zstrong_ew_12oh_kk04_frac = result.frac
             abund[index].zstrong_ew_12oh_kk04_avg   = result.log12oh_avg
             abund[index].zstrong_ew_12oh_kk04_upper = result.log12oh_upper
             abund[index].zstrong_ew_12oh_kk04_lower = result.log12oh_lower
             abund[index].zstrong_ew_logu_kk04_avg   = result.logu_avg
             abund[index].zstrong_ew_logu_kk04_upper = result.logu_upper
             abund[index].zstrong_ew_logu_kk04_lower = result.logu_lower

             abund[index].zstrong_ew_converge_kk04_upper = result.converge_upper
             abund[index].zstrong_ew_converge_kk04_lower = result.converge_lower
             
             abund[index].zstrong_ew_12oh_kk04_avg_err   = result.log12oh_avg_err
             abund[index].zstrong_ew_12oh_kk04_upper_err = result.log12oh_upper_err
             abund[index].zstrong_ew_12oh_kk04_lower_err = result.log12oh_lower_err
             abund[index].zstrong_ew_logu_kk04_avg_err   = result.logu_avg_err
             abund[index].zstrong_ew_logu_kk04_upper_err = result.logu_upper_err
             abund[index].zstrong_ew_logu_kk04_lower_err = result.logu_lower_err
          endif

; ---------------------------------------------------------------------------    
; Tremonti et al. (2004)
          if (keyword_set(silent) eq 0) then splog, '   T04...'
          index = where(abund.zstrong_ew_r23 gt -900.0,nindex)
          if (nindex ne 0L) then begin
             result = monte_log12oh_t04(abund[index].zstrong_ew_r23,$
               line=ewr23lines[index],nmonte=nmonte)
             abund[index].zstrong_ew_12oh_t04     = result.log12oh
             abund[index].zstrong_ew_12oh_t04_err = result.log12oh_err
          endif
       endif                    ; close the JUSTRATIOS case
    endif                       ; close the JUSTEW case

; ###########################################################################
; compute flux-ratios and abundances
    if keyword_set(justflux) then begin
       if (keyword_set(silent) eq 0) then begin
          if keyword_set(justratios) then splog, 'Computing line-flux ratios' else $
            splog, 'Computing line-ratios and abundances from fluxes:'
       endif

; if /USE_4959 then use the measured [O III] 4959 line, otherwise use
; the intrinsic ratio, 2.984; note that if LINE has not been
; reddening-corrected, then the assumption of the intrinsic ratio will
; be wrong by a little bit

       fluxline = struct_addtags(line,replicate({oiii_4959_5007: [0.0,-2.0]},nspec))
       if keyword_set(use_4959) then begin
          if tag_exist(fluxline[0],'OIII_4959') and tag_exist(fluxline[0],'OIII_5007') then begin
             doit = where((fluxline.oiii_4959[1] gt 0.0) and (fluxline.oiii_5007[1] gt 0.0),ndoit)
             if (ndoit ne 0L) then begin
                fluxline[doit].oiii_4959_5007[0] = fluxline[doit].oiii_4959[0]+fluxline[doit].oiii_5007[0]
                fluxline[doit].oiii_4959_5007[1] = sqrt(fluxline[doit].oiii_4959[1]^2+fluxline[doit].oiii_5007[1]^2)
             endif
          endif 
       endif else begin
          if (tag_exist(fluxline[0],'OIII_5007')) then begin
             doit = where((fluxline.oiii_5007[1] gt 0.0),ndoit)
             if (ndoit ne 0L) then fluxline[doit].oiii_4959_5007 = ocor*fluxline[doit].oiii_5007
          endif 
       endelse

; R23
       if (tag_exist(fluxline[0],'OII_3727') and tag_exist(fluxline[0],'OIII_4959_5007') and $
         tag_exist(fluxline[0],'H_BETA')) then begin
          lineratio, fluxline, ['OII_3727','OIII_4959_5007'], 'H_BETA', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_r23 = x
             abund[index].zstrong_r23_err = xerr
          endif
       endif 

; O32
       if (tag_exist(fluxline[0],'OII_3727') and tag_exist(fluxline[0],'OIII_4959_5007')) then begin
          lineratio, fluxline, 'OIII_4959_5007', 'OII_3727', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_o32 = x
             abund[index].zstrong_o32_err = xerr
          endif
       endif 

; P-parameter (Pilyugin 2001)
       if tag_exist(fluxline[0],'OIII_4959_5007') and tag_exist(fluxline[0],'OII_3727') then begin
          lineratio, fluxline, 'OIII_4959_5007', ['OII_3727','OIII_4959_5007'], $
            '', '', x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_p = x
             abund[index].zstrong_p_err = xerr
          endif
       endif 

;; ([O III]/Hb)/([N II]/Ha) - Pettini & Pagel (2004)
;       if (tag_exist(fluxline[0],'NII_6584') and tag_exist(fluxline[0],'OIII_5007') and $
;         tag_exist(fluxline[0],'H_BETA') and tag_exist(fluxline[0],'H_ALPHA')) then begin
;          lineratio, fluxline, 'OIII_5007', 'H_BETA', 'NII_6584', 'H_ALPHA', $
;            x1, x1err, x2, x2err, index=index, nindex=nindex, snrcut=snrcut_abundance
;          if (nindex ne 0L) then begin
;             abund[index].zstrong_oiiinii = x1 - x2
;             abund[index].zstrong_oiiinii_err = sqrt(x1err^2 + x2err^2)
;          endif
;       endif

; [N II]/Ha    
       if (tag_exist(fluxline[0],'NII_6584') and tag_exist(fluxline[0],'H_ALPHA')) then begin
          lineratio, fluxline, 'NII_6584', 'H_ALPHA', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
          if (nindex ne 0L) then begin
             abund[index].zstrong_niiha = x
             abund[index].zstrong_niiha_err = xerr
          endif
       endif

; [N II]/[O II]
       if (tag_exist(fluxline[0],'NII_6584') and tag_exist(fluxline[0],'OII_3727')) then begin
          lineratio, fluxline, 'NII_6584', 'OII_3727', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
          if (nindex ne 0L) then begin
             abund[index].zstrong_niioii = x
             abund[index].zstrong_niioii_err = xerr
          endif
       endif

;; [N II]/[S II]
;       if (tag_exist(fluxline[0],'NII_6584') and tag_exist(fluxline[0],'SII_6716') and $
;         tag_exist(fluxline[0],'SII_6731')) then begin
;          lineratio, fluxline, 'NII_6584', ['SII_6716','SII_6731'], '', '', $
;            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
;          if (nindex ne 0L) then begin
;             abund[index].zstrong_niisii = x
;             abund[index].zstrong_niisii_err = xerr
;          endif
;       endif

;; [O II]/[O III]
;       if (tag_exist(fluxline[0],'OII_3727') and tag_exist(fluxline[0],'OIII_5007')) then begin
;          lineratio, fluxline, 'OII_3727', 'OIII_5007', '', '', $
;            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
;          if (nindex ne 0L) then begin
;             abund[index].zstrong_oiioiii = x
;             abund[index].zstrong_oiioiii_err = xerr
;          endif
;       endif

;; [O III]/H-beta
;       if (tag_exist(fluxline[0],'OIII_5007') and tag_exist(fluxline[0],'H_BETA')) then begin
;          lineratio, fluxline, 'OIII_5007', 'H_BETA', '', '', $
;            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
;          if (nindex ne 0L) then begin
;             abund[index].zstrong_oiiihb = x
;             abund[index].zstrong_oiiihb_err = xerr
;          endif
;       endif

;; S23
;       if (tag_exist(fluxline[0],'SII_6716') and tag_exist(fluxline[0],'SII_6731') and $
;         tag_exist(fluxline[0],'SIII_9069') and tag_exist(fluxline[0],'SIII_9532') and $
;         tag_exist(fluxline[0],'H_BETA')) then begin
;          lineratio, fluxline, ['SII_6716','SII_6731','SIII_9069','SIII_9532'], 'H_BETA', '', '', $
;            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
;          if (nindex ne 0L) then begin
;             abund[index].zstrong_S23 = x
;             abund[index].zstrong_S23_err = xerr
;          endif
;       endif

;; [S II] 6716 / [S II] 6731
;       if (tag_exist(fluxline[0],'SII_6716') and tag_exist(fluxline[0],'SII_6731')) then begin
;          lineratio, fluxline, 'SII_6716', 'SII_6731', '', '', $
;            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
;          if (nindex ne 0L) then begin
;             abund[index].zstrong_sii = x
;             abund[index].zstrong_sii_err = xerr
;          endif
;       endif
;
;; [O II] 3726 / [O II] 3729
;       if (tag_exist(fluxline[0],'OII_3726') and tag_exist(fluxline[0],'OII_3729')) then begin
;          lineratio, fluxline, 'OII_3726', 'OII_3729', '', '', $
;            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
;          if (nindex ne 0L) then begin
;             abund[index].zstrong_oii = x
;             abund[index].zstrong_oii_err = xerr
;          endif
;       endif

       r23lines = get_r23lines(fluxline)
       fluxline = 0 ; save memory
       
       if (keyword_set(justratios) eq 0) then begin ; compute abundances

; ---------------------------------------------------------------------------    
; McGaugh (1991) as fitted in Kobulnicky et al. (1999)
          if (not keyword_set(silent)) then splog, '   M91...'
          index = where((abund.zstrong_o32 gt -900.0) and (abund.zstrong_r23 gt -900.0),nindex)
          if (nindex ne 0L) then begin

             result = monte_log12oh_m91(abund[index].zstrong_r23,$
               abund[index].zstrong_o32,line=r23lines[index],$
               nmonte=nmonte)

;            abund[index].zstrong_12oh_m91_frac = result.frac
             abund[index].zstrong_12oh_m91_avg    = result.log12oh_avg
             abund[index].zstrong_12oh_m91_upper  = result.log12oh_upper
             abund[index].zstrong_12oh_m91_lower  = result.log12oh_lower

             abund[index].zstrong_12oh_m91_avg_err    = result.log12oh_avg_err
             abund[index].zstrong_12oh_m91_upper_err  = result.log12oh_upper_err
             abund[index].zstrong_12oh_m91_lower_err  = result.log12oh_lower_err

          endif

; ---------------------------------------------------------------------------    
; Kobulnicky & Kewley (2004) - Flux/R23
          if (keyword_set(silent) eq 0) then splog, '   KK04...'
          index = where((abund.zstrong_o32 gt -900.0) and (abund.zstrong_r23 gt -900.0),nindex)
          if (nindex ne 0L) then begin

             result = monte_log12oh_kk04(abund[index].zstrong_r23,$
               abund[index].zstrong_o32,line=r23lines[index],$
               nmonte=nmonte)
             
;            abund[index].zstrong_12oh_kk04_frac = result.frac
             abund[index].zstrong_12oh_kk04_avg   = result.log12oh_avg
             abund[index].zstrong_12oh_kk04_upper = result.log12oh_upper
             abund[index].zstrong_12oh_kk04_lower = result.log12oh_lower
             abund[index].zstrong_logu_kk04_avg   = result.logu_avg
             abund[index].zstrong_logu_kk04_upper = result.logu_upper
             abund[index].zstrong_logu_kk04_lower = result.logu_lower
             
             abund[index].zstrong_converge_kk04_upper = result.converge_upper
             abund[index].zstrong_converge_kk04_lower = result.converge_lower

             abund[index].zstrong_12oh_kk04_avg_err   = result.log12oh_avg_err
             abund[index].zstrong_12oh_kk04_upper_err = result.log12oh_upper_err
             abund[index].zstrong_12oh_kk04_lower_err = result.log12oh_lower_err
             abund[index].zstrong_logu_kk04_avg_err   = result.logu_avg_err
             abund[index].zstrong_logu_kk04_upper_err = result.logu_upper_err
             abund[index].zstrong_logu_kk04_lower_err = result.logu_lower_err
          endif

;; ---------------------------------------------------------------------------    
;; Zaritsky et al. (1994)
;          index = where((abund.zstrong_r23 gt -900.0),nindex)
;          if (nindex ne 0L) then begin
;             result = monte_log12oh_zkh94(abund[index].zstrong_r23,$
;               line=r23lines[index],nmonte=nmonte)
;             abund[index].zstrong_12oh_zkh94     = result.log12oh
;             abund[index].zstrong_12oh_zkh94_err = result.log12oh_err
;          endif

; ---------------------------------------------------------------------------    
; Tremonti et al. (2004)
          index = where(abund.zstrong_r23 gt -900.0,nindex)
          if (nindex ne 0L) then begin
             result = monte_log12oh_t04(abund[index].zstrong_r23,$
               line=r23lines[index],nmonte=nmonte)
             abund[index].zstrong_12oh_t04     = result.log12oh
             abund[index].zstrong_12oh_t04_err = result.log12oh_err
          endif
          
;; ---------------------------------------------------------------------------    
;; Pettini & Pagel (2004) [N II]/Ha calibration
;          if (not keyword_set(silent)) then splog, '   N2...'
;          index = where((abund.zstrong_niiha ge -2.5) and (abund.zstrong_niiha lt -0.3),nindex)
;          if (nindex ne 0L) then begin
;
;             c = [8.90D,0.57D]
;;            cerr = [0.18D,0.0D]
;             cerr = [0.0D,0.0D]
;             
;             x = abund[index].zstrong_niiha
;             xerr = abund[index].zstrong_niiha_err
;             
;             abund[index].zstrong_12oh_niiha_pettini = c[0] + c[1] * x
;             abund[index].zstrong_12oh_niiha_pettini_err = sqrt( cerr[0]^2 + $
;               (c[1] * xerr)^2.0 + (cerr[1] * x)^2.0 )
;
;             x = 0 & xerr = 0   ; save memory
;
;          endif
;          
;; ---------------------------------------------------------------------------    
;; Pettini & Pagel (2004) ([O III]/Hb)/([N II]/Ha) calibration
;          if (not keyword_set(silent)) then splog, '   O3N2...'
;          index = where((abund.zstrong_oiiinii ge -1.0) and (abund.zstrong_oiiinii le 1.9),nindex)
;          if (nindex ne 0L) then begin
;
;             c = [8.73D,-0.32D]
;             cerr = [0.0D,0.0D]
;
;             x = abund[index].zstrong_oiiinii
;             xerr = abund[index].zstrong_oiiinii_err
;
;             abund[index].zstrong_12oh_oiiinii_pettini = c[0] + c[1] * x
;             abund[index].zstrong_12oh_oiiinii_pettini_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )
;
;             x = 0 & xerr = 0   ; save memory
;
;          endif

;; [S II] electron density [cm-3]
;          index = where((abund.zstrong_sii gt -900.0),nindex)
;          if (nindex ne 0L) and keyword_set(getdensity) then begin
;             if (not keyword_set(silent)) then splog, 'Computing [S II] electron densities...'
;             temden = im_temden('s_ii',abund[index].zstrong_sii,$
;               ratio_err=abund[index].zstrong_sii_err,nmonte=nmonte)
;             abund[index].zstrong_sii_dens     = temden.dens
;             abund[index].zstrong_sii_dens_err = temden.dens_err
;          endif
;
;; [O II] electron density [cm-3]
;          index = where((abund.zstrong_oii gt -900.0),nindex)
;          if (nindex ne 0L) and keyword_set(getdensity) then begin
;             if (not keyword_set(silent)) then splog, 'Computing [O II] electron densities...'
;             temden = im_temden('o_ii_dens',abund[index].zstrong_oii,$
;               ratio_err=abund[index].zstrong_oii_err,nmonte=nmonte)
;             abund[index].zstrong_oii_dens     = temden.dens
;             abund[index].zstrong_oii_dens_err = temden.dens_err
;          endif
       endif                    ; close the JUSTRATIOS case
    endif                       ; close the JUSTFLUX case
       
return, abund
end
