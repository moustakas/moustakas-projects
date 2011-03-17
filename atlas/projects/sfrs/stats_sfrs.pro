pro stats_sfrs, atlasdust, atlasnodust, nfgsdust, nfgsnodust, latex=latex
; jm05feb01uofa

    path = atlas_path(/projects)+'sfrs/'
    
    if (n_elements(atlasdust) eq 0L) then atlasdust = read_integrated(linefitnodust=atlasnodust,/snrcuts,/hiionly)
    if (n_elements(nfgsdust) eq 0L) then nfgsdust = read_nfgs(linefitnodust=nfgsnodust,/snrcuts,/hiionly)
    
    if n_elements(snrcut) eq 0L then snrcut = 3.0

    irconst = alog10(4.5D-44)   ; K98 conversion L(IR) --> SFR(IR)
    haconst = alog10(7.9D-42)   ; K98 conversion L(IR) --> SFR(IR)

    result = {$

      LHaobs_LIR_median:               '', $
      LHaobs_LIR_error:                '', $
      LHacor_LIR_median:              '', $
      LHacor_LIR_error:               '', $

; ----------      
; H-beta
; ----------      
      
      LHbobs_LHaobs_NoAbs_median: '', $
      LHbobs_LHaobs_NoAbs_error:  '', $
      LHbobs_LHaobs_Abs_median: '', $
      LHbobs_LHaobs_Abs_error:  '', $
      LHbobs_LHacor_Abs_median: '', $
      LHbobs_LHacor_Abs_error:  '', $
      LHbobs_LHacor_NoAbs_median: '', $
      LHbobs_LHacor_NoAbs_error:  '', $

      SFR_LHbobs_NoAbs_median: '', $ ; no figure for this
      SFR_LHbobs_NoAbs_error:  '', $
      SFR_LHbobs_Abs_median: '', $
      SFR_LHbobs_Abs_error:  '', $

; ----------      
; [O II]
; ----------      

      Loiiobs_LHaobs_median: '', $
      Loiiobs_LHaobs_error:  '', $
      Loiiobs_LHacor_median: '', $
      Loiiobs_LHacor_error:  '', $
      Loiicor_LHacor_median: '', $
      Loiicor_LHacor_error:  '', $

      SFR_Loiiobs_median: '', $
      SFR_Loiiobs_error:  '', $
      SFR_Loiicor_median: '', $
      SFR_Loiicor_error:  '', $

      Loiicor_LHacor_lowZ_median: '', $
      Loiicor_LHacor_lowZ_error:  '', $
      Loiicor_LHacor_midZ_median: '', $
      Loiicor_LHacor_midZ_error:  '', $
      Loiicor_LHacor_highZ_median: '', $
      Loiicor_LHacor_highZ_error:  '', $

      SFR_Loiicor_lowZ_median: '', $
      SFR_Loiicor_lowZ_error:  '', $
      SFR_Loiicor_midZ_median: '', $
      SFR_Loiicor_midZ_error:  '', $
      SFR_Loiicor_highZ_median: '', $
      SFR_Loiicor_highZ_error:  '', $

; ----------      
; U-band
; ----------      

      LUobs_LHaobs_median: '', $
      LUobs_LHaobs_error:  '', $
      LUobs_LHacor_median: '', $
      LUobs_LHacor_error:  '', $
      LUcor_LHacor_median: '', $
      LUcor_LHacor_error:  '', $

      SFR_LUobs_median: '', $
      SFR_LUobs_error:  '', $
      SFR_LUcor_median: '', $
      SFR_LUcor_error:  ''}

; ---------------------------------------------------------------------------    
; SFR(Ha)
; ---------------------------------------------------------------------------    

    indx = where((atlasnodust.ir_flux gt -900) and (atlasnodust.ebv_hahb_err gt 0.0),nindx)
    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]

    indxnfgs = where((nfgsdust.ir_flux gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)
    lirnfgs = nfgsdust[indxnfgs].ir_flux

; ##################################################    
; Observed
; ##################################################    
    
    ha = atlasdust[indx].h_alpha[0]
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]

    y = alog10(ha/lir)
    ynfgs = alog10(hanfgs/lirnfgs)

    stats = im_stats([y,ynfgs])
    result.LHaobs_LIR_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.LHaobs_LIR_error = strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
; ##################################################    
; Corrected
; ##################################################    

    ha = atlasnodust[indx].h_alpha[0]
    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]

    y = alog10(ha/lir)
    ynfgs = alog10(hanfgs/lirnfgs)

    stats = im_stats([y,ynfgs])
    result.LHacor_LIR_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.LHacor_LIR_error = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ---------------------------------------------------------------------------    
; SFR(Hb)
; ---------------------------------------------------------------------------    

    ewcut = 1.0
    
    indx = where((atlasdust.rc3_B_lum gt -900) and (atlasnodust.ebv_hahb_err gt 0.0) and $
      (atlasdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    indxnfgs = where((nfgsdust.rc3_B_lum gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0) and $
      (nfgsdust.h_beta_ew_uncor[0] ge ewcut),nindxnfgs)

; ##################################################    
; 1 - No Abs, Hb Observed, Ha Observed
; ##################################################    

    hb = atlasdust[indx].h_beta_uncor[0]
    ha = atlasdust[indx].h_alpha[0]
    y = alog10(hb/ha)
    
    hbnfgs = nfgsdust[indxnfgs].h_beta_uncor[0]
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    ynfgs = alog10(hbnfgs/hanfgs)

    stats = im_stats([y,ynfgs])
    result.LHbobs_LHaobs_NoAbs_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.LHbobs_LHaobs_NoAbs_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ##################################################    
; 2 - Abs, Hb observed, Ha observed
; ##################################################    

    hb = atlasdust[indx].h_beta[0]
    ha = atlasdust[indx].h_alpha[0]
    y = alog10(hb/ha)
    
    hbnfgs = nfgsdust[indxnfgs].h_beta[0]
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    ynfgs = alog10(hbnfgs/hanfgs)

    stats = im_stats([y,ynfgs])
    result.LHbobs_LHaobs_Abs_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.LHbobs_LHaobs_Abs_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ##################################################    
; 3 - Abs, Hb observed, Ha corrected
; ##################################################    

    hb = atlasdust[indx].h_beta[0]
    ha = atlasnodust[indx].h_alpha[0]
    y = alog10(hb/ha)
    
    hbnfgs = nfgsdust[indxnfgs].h_beta[0]
    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    ynfgs = alog10(hbnfgs/hanfgs)
    
    stats = im_stats([y,ynfgs])

    result.LHbobs_LHacor_Abs_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.LHbobs_LHacor_Abs_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    result.SFR_LHbobs_Abs_median    = strtrim(string(haconst-stats.median,format='(F12.2)'),2)
    result.SFR_LHbobs_Abs_error     = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ##################################################    
; 4 - Abs, Hb observed, Ha corrected - No Figure
; ##################################################    

    hb = atlasdust[indx].h_beta_uncor[0]
    ha = atlasnodust[indx].h_alpha[0]
    y = alog10(hb/ha)
    
    hbnfgs = nfgsdust[indxnfgs].h_beta_uncor[0]
    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    ynfgs = alog10(hbnfgs/hanfgs)
    
    stats = im_stats([y,ynfgs])

    result.LHbobs_LHacor_NoAbs_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.LHbobs_LHacor_NoAbs_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    result.SFR_LHbobs_NoAbs_median    = strtrim(string(haconst-stats.median,format='(F12.2)'),2)
    result.SFR_LHbobs_NoAbs_error     = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ---------------------------------------------------------------------------    
; SFR(OII)
; ---------------------------------------------------------------------------    

    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.rc3_B_lum gt -900.0) and (atlasnodust.ehbha_err gt 0.0),nindx)

    indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
      (nfgsnodust.rc3_B_lum gt -900.0) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

; ##########################################################
; 1 - [O II] observed, Ha observed
; ##########################################################

    y1 = atlasdust[indx].oii_3727[0]
    y2 = atlasdust[indx].h_alpha[0]
    y = alog10(y1/y2)

    y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
    y2nfgs = nfgsdust[indxnfgs].h_alpha[0]
    ynfgs = alog10(y1nfgs/y2nfgs)

    stats = im_stats([y,ynfgs])
    result.Loiiobs_LHaobs_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.Loiiobs_LHaobs_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ##########################################################
; 2 - [O II] observed, Ha corrected
; ##########################################################

    y1 = atlasdust[indx].oii_3727[0]
    y2 = atlasnodust[indx].h_alpha[0]
    y = alog10(y1/y2)

    y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    ynfgs = alog10(y1nfgs/y2nfgs)

    stats = im_stats([y,ynfgs])

    result.Loiiobs_LHacor_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.Loiiobs_LHacor_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    result.SFR_Loiiobs_median = strtrim(string(haconst-stats.median,format='(F12.2)'),2)
    result.SFR_Loiiobs_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ##########################################################
; 3 - [O II] corrected, Ha corrected
; ##########################################################

    y1 = atlasnodust[indx].oii_3727[0]
    y2 = atlasnodust[indx].h_alpha[0]
    y = alog10(y1/y2)

    y1nfgs = nfgsnodust[indxnfgs].oii_3727[0]
    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    ynfgs = alog10(y1nfgs/y2nfgs)

    stats = im_stats([y,ynfgs])

    result.Loiicor_LHacor_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.Loiicor_LHacor_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    result.SFR_Loiicor_median = strtrim(string(haconst-stats.median,format='(F12.2)'),2)
    result.SFR_Loiicor_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ----------
; Metallicity Sensitivity; use the O3N2 metallicity calibration, and
; the N2 calibration for objects below the range of the O3N2
; calibration
; ----------

    w1 = where((atlasnodust.z_12oh_o3n2_pettini gt -900),nw1)
    w2 = where((atlasnodust.z_12oh_o3n2_pettini lt -900) and (atlasnodust.z_12oh_n2_denicolo gt -900),nw2)
    if (nw2 ne 0L) then w = [w1,w2]
    newindx = cmset_op(indx,'AND',w)

    w1nfgs = where((nfgsnodust.z_12oh_o3n2_pettini gt -900),nw1nfgs)
    w2nfgs = where((nfgsnodust.z_12oh_o3n2_pettini lt -900) and (nfgsnodust.z_12oh_n2_denicolo gt -900),nw2nfgs)
    if (nw2nfgs ne 0L) then wnfgs = [w1nfgs,w2nfgs]
    newindxnfgs = cmset_op(indxnfgs,'AND',wnfgs)

    y1 = atlasnodust[newindx].oii_3727[0]
    y2 = atlasnodust[newindx].h_alpha[0]
    y = alog10(y1/y2)
    Z = atlasnodust[newindx].z_12oh_o3n2_pettini
    badZ = where(Z lt -900,nbadZ)
    if (nbadZ ne 0L) then Z[badZ] = atlasnodust[newindx[badZ]].z_12oh_n2_denicolo

    y1nfgs = nfgsnodust[newindxnfgs].oii_3727[0]
    y2nfgs = nfgsnodust[newindxnfgs].h_alpha[0]
    ynfgs = alog10(y1nfgs/y2nfgs)
    Znfgs = nfgsnodust[newindxnfgs].z_12oh_o3n2_pettini
    badZnfgs = where(Znfgs lt -900,nbadZnfgs)
    if (nbadZnfgs ne 0L) then Znfgs[badZnfgs] = nfgsnodust[newindxnfgs[badZnfgs]].z_12oh_n2_denicolo
    
; ##########################################################
; 1 - [O II] corrected, Ha corrected, low-Z
; ##########################################################

    lowZ = where(Z lt 8.14)
    lowZnfgs = where(Znfgs lt 8.14) ; <-- None!

    if (lowZnfgs[0] ne -1L) then stats = im_stats([y[lowZ],ynfgs[lowZnfgs]]) else $
      stats = im_stats(y[lowZ])

    result.Loiicor_LHacor_lowZ_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.Loiicor_LHacor_lowZ_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    result.SFR_Loiicor_lowZ_median = strtrim(string(haconst-stats.median,format='(F12.2)'),2)
    result.SFR_Loiicor_lowZ_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ##########################################################
; 2 - [O II] corrected, Ha corrected, mid-Z
; ##########################################################

    MidZ = where((Z gt 8.14) and (Z lt 8.68))
    MidZnfgs = where((Znfgs gt 8.14) and (Znfgs lt 8.68))
    
    stats = im_stats([y[MidZ],ynfgs[MidZnfgs]])

    result.Loiicor_LHacor_MidZ_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.Loiicor_LHacor_MidZ_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    result.SFR_Loiicor_MidZ_median = strtrim(string(haconst-stats.median,format='(F12.2)'),2)
    result.SFR_Loiicor_MidZ_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ##########################################################
; 3 - [O II] corrected, Ha corrected, high-Z
; ##########################################################

    HighZ = where(Z gt 8.68)
    HighZnfgs = where(Znfgs gt 8.68)
    
    stats = im_stats([y[HighZ],ynfgs[HighZnfgs]])

    result.Loiicor_LHacor_HighZ_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.Loiicor_LHacor_HighZ_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    result.SFR_Loiicor_HighZ_median = strtrim(string(haconst-stats.median,format='(F12.2)'),2)
    result.SFR_Loiicor_HighZ_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ---------------------------------------------------------------------------    
; SFR(U)
; ---------------------------------------------------------------------------    

    Uinfo = im_filterspecs(filterlist='bessell_U.par')
    Uconstant = Uinfo.weff*Uinfo.vega_flam

    indx = where((atlasdust.d4000_narrow_model[1] gt 0.0) and (atlasdust.U gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0),nindx)
    indxnfgs = where((nfgsdust.d4000_narrow_model[1] gt 0.0) and (nfgsdust.U gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0),nindxnfgs)

; ##########################################################
; 1 - U observed, Ha observed
; ##########################################################

    y1 = Uconstant*10^(-0.4*atlasdust[indx].U)
    y2 = atlasdust[indx].h_alpha[0]
    y = alog10(y1/y2)

    y1nfgs = Uconstant*10^(-0.4*nfgsdust[indxnfgs].U)
    y2nfgs = nfgsdust[indxnfgs].h_alpha[0]
    ynfgs = alog10(y1nfgs/y2nfgs)
    
    stats = im_stats([y,ynfgs])
    result.LUobs_LHaobs_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.LUobs_LHaobs_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ##########################################################
; 2 - U observed, Ha corrected
; ##########################################################

    y1 = Uconstant*10^(-0.4*atlasdust[indx].U)
    y2 = atlasnodust[indx].h_alpha[0]
    y = alog10(y1/y2)

    y1nfgs = Uconstant*10^(-0.4*nfgsdust[indxnfgs].U)
    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    ynfgs = alog10(y1nfgs/y2nfgs)

    stats = im_stats([y,ynfgs])

    result.LUobs_LHacor_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.LUobs_LHacor_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    result.SFR_LUobs_median = strtrim(string(haconst-stats.median,format='(F12.2)'),2)
    result.SFR_LUobs_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ##########################################################
; 3 - U corrected, Ha corrected
; ##########################################################

    dust_fraction = 0.5

    AU = atlasnodust[indx].ebv_hahb*k_lambda(Uinfo.weff,/charlot)*dust_fraction
    AUnfgs = nfgsnodust[indxnfgs].ebv_hahb*k_lambda(Uinfo.weff,/charlot)*dust_fraction

    y1 = Uconstant*10^(-0.4*atlasdust[indx].U)*10^(0.4*AU)
    y2 = atlasnodust[indx].h_alpha[0]
    y = alog10(y1/y2)

    y1nfgs = Uconstant*10^(-0.4*nfgsdust[indxnfgs].U)*10^(0.4*AUnfgs)
    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    ynfgs = alog10(y1nfgs/y2nfgs)

    stats = im_stats([y,ynfgs])

    result.LUcor_LHacor_median = strtrim(string(stats.median,format='(F12.2)'),2)
    result.LUcor_LHacor_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    result.SFR_LUcor_median = strtrim(string(haconst-stats.median,format='(F12.2)'),2)
    result.SFR_LUcor_error  = strtrim(string(stats.sig68mean,format='(F12.2)'),2)

; ---------------------------------------------------------------------------
; generate a latex table
; ---------------------------------------------------------------------------

    paperpath = atlas_path(/papers)+'sfrs/'
    if keyword_set(latex) then begin
    
       openw, lun, paperpath+'sfr_conversions.tex', /get_lun
       printf, lun, '\begin{deluxetable}{ccccccl}'
       printf, lun, '\tabletypesize{\scriptsize}'
       printf, lun, '\tablecolumns{7}'
       printf, lun, '\tablecaption{Star-Formation Rate Indicator Conversions}'
       printf, lun, '\tablewidth{0in}'
       printf, lun, '\tablehead{'
       printf, lun, '\colhead{Indicator} & '
       printf, lun, '\colhead{$\lambda$} & '
       printf, lun, '\colhead{Extinction} & '
       printf, lun, '\colhead{$\log\, [L(\lambda)/L(\ha)]$\tablenotemark{b}} & '
       printf, lun, '\colhead{$\log\, [\psi(\ha)/L(\lambda)]$\tablenotemark{c}} & '
       printf, lun, '\colhead{Figure} & '
       printf, lun, '\colhead{Remarks} \\ '
       printf, lun, '\colhead{} & '
       printf, lun, '\colhead{(\AA)} & '
       printf, lun, '\colhead{Corrected?\tablenotemark{a}} & '
       printf, lun, '\colhead{} & '
       printf, lun, '\colhead{$(\sfrunits)/(\lunits)$} & '
       printf, lun, '\colhead{Number} & '
       printf, lun, '\colhead{}'
       printf, lun, '}'
       printf, lun, '\startdata'

; H-beta

       printf, lun, '\hb & 4861 & No & '+$
         '$'+result.LHbobs_LHacor_Abs_median+'\pm'+result.LHbobs_LHacor_Abs_error+'$ & '+$
         '$'+result.SFR_LHbobs_Abs_median+'\pm'+result.SFR_LHbobs_Abs_error+'$ & '+$
         ' \ref{fig:LB_hbha}\emph{c} & Absorption-Corrected \\'
;      printf, lun, ' &  & No & '+$
;        '$'+result.LHbobs_LHacor_NoAbs_median+'\pm'+result.LHbobs_LHacor_NoAbs_error+'$ & '+$
;        '$'+result.SFR_LHbobs_NoAbs_median+'\pm'+result.SFR_LHbobs_NoAbs_error+'$ & '+$
;        '   & No Absorption Correction \\'

; [O II] 3727

       printf, lun, '\oii & 3727 & No & '+$
         '$'+result.Loiiobs_LHacor_median+'\pm'+result.Loiiobs_LHacor_error+'$ & '+$
         '$'+result.SFR_Loiiobs_median+'\pm'+result.SFR_Loiiobs_error+'$ & '+$
         ' \ref{fig:LB_oiiha}\emph{b} & Full Sample \\'
       printf, lun, '     &      & Yes & '+$
         '$'+result.Loiicor_LHacor_median+'\pm'+result.Loiicor_LHacor_error+'$ & '+$
         '$'+result.SFR_Loiicor_median+'\pm'+result.SFR_Loiicor_error+'$ & '+$
         ' \ref{fig:LB_oiiha}\emph{c} & Full Sample \\'
       printf, lun, '     &      & Yes & '+$
         '$'+result.Loiicor_LHacor_lowz_median+'\pm'+result.Loiicor_LHacor_lowz_error+'$ & '+$
         '$'+result.SFR_Loiicor_lowz_median+'\pm'+result.SFR_Loiicor_lowz_error+'$ & '+$
         ' \ref{fig:OH12_oiiha} & $\logoh\lesssim8.14$ \\'
       printf, lun, '     &      & Yes & '+$
         '$'+result.Loiicor_LHacor_midz_median+'\pm'+result.Loiicor_LHacor_midz_error+'$ & '+$
         '$'+result.SFR_Loiicor_midz_median+'\pm'+result.SFR_Loiicor_midz_error+'$ & '+$
         ' \ref{fig:OH12_oiiha} & $8.14\lesssim\logoh\lesssim8.68$ \\'
       printf, lun, '     &      & Yes & '+$
         '$'+result.Loiicor_LHacor_highz_median+'\pm'+result.Loiicor_LHacor_highz_error+'$ & '+$
         '$'+result.SFR_Loiicor_highz_median+'\pm'+result.SFR_Loiicor_highz_error+'$ & '+$
         ' \ref{fig:OH12_oiiha} & $\logoh\gtrsim8.68$ \\'

; U-band

       printf, lun, '$U$-band & 3585 & No & '+$
         '$'+result.LUobs_LHacor_median+'\pm'+result.LUobs_LHacor_error+'$ & '+$
         '$'+result.SFR_LUobs_median+'\pm'+result.SFR_LUobs_error+'$ & '+$
         ' \ref{fig:D4000_UHa}\emph{b} &  \\'
       printf, lun, '     &      & Yes & '+$
         '$'+result.LUcor_LHacor_median+'\pm'+result.LUcor_LHacor_error+'$ & '+$
         '$'+result.SFR_LUcor_median+'\pm'+result.SFR_LUcor_error+'$ & '+$
         ' \ref{fig:D4000_UHa}\emph{c} & $\eta=0.5$ (see text) \\'

       printf, lun, '\enddata'
       printf, lun, '\label{table:sfr_conversions}'
       printf, lun, '\tablenotetext{a}{This column indicates whether or not the indicator has been '+$
         'corrected for extinction using the Balmer decrement.}'
       printf, lun, '\tablenotetext{b}{Median luminosity ratio and $1\sigma$ scatter.}'
       printf, lun, '\tablenotetext{c}{Scalar factor and uncertainty to convert from luminosity '+$
         'to \ha{} star-formation rate, $\psi(\ha)$.}'
       printf, lun, '%\tablecomments{Comments}'
       printf, lun, '%\tablerefs{}'
       printf, lun, '\end{deluxetable}'
       free_lun, lun

    endif

    mwrfits, result, path+'stats_sfrs.fits', /create
    
;   help, result, /str

stop    
    
return
end
    
