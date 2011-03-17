pro write_atlas_hcn_sample, result, write=write, thumbs=thumbs, gao=gao
; jm06jul08uofa - match my atlas with Gao & Solomon (2004)
; jm07sep24nyu - total re-write and clean up; see below for detailed
;                comments 
; jm08apr07nyu - minor revisions and clean-up
    
; THUMBS - generate thumbnail visualizations

; read the integrated and nuclear spectral atlas

    hcnpath = atlas_path(/projects)+'hcn/' ; CUSTOMIZE!
    atlas = read_integrated()
    nuclear = read_nuclear()

; set up some path names and constants

    gao = 1L ; use either Gao or Baan
    baan = 0L
    
    mpc2cm = 3.086D24     ; [cm/Mpc]
    light = 2.99792458D18 ; speed of light [A/s]
    lsun = 3.826D33       ; bolometric solar luminosity [erg/s]

    dust_snrcut = 2.0
    class_snrcut = 1.0

    kenn_sfr_ha_const = 7.9D-42 ; Salpeter IMF
    kenn_sfr_ir_const = 4.5D-44 ; Salpeter IMF
    
    kenn_factor = 0.031D ; 0.038D
    kenn_factor_err = 0.006D

    lhcnerr = 0.1 ; NOTE!
    lcoerr = 0.1  ; NOTE!
    
; ###########################################################################    
; match against the integrated spectral atlas and write out 
; ###########################################################################    

; read the Baan et al. or the Gao & Solomon sample

    if keyword_set(baan) then begin
    
       hcn_sample = 'Baan et al.'
       hcn1 = mrdfits(hcnpath+'baan.ned.fits.gz',1,/silent)

       ra = 15.0*im_hms2dec(hcn1.ra)
       de = im_hms2dec(hcn1.dec)

       indx = speclinefit_locate(atlas,'IIZW096')
       atlas[indx[[1,2]]].ra = '-99:00:00'
       atlas[indx[[1,2]]].dec = '-99:00:00'

       raref = 15.0*im_hms2dec(atlas.ra)
       deref = im_hms2dec(atlas.dec)

       searchrad = 3.0

    endif
    
    if keyword_set(gao) then begin

       hcn_sample = 'Gao & Solomon'
       hcn1 = mrdfits(hcnpath+'gao.ned.fits.gz',1,/silent)

       ra = 15.0*im_hms2dec(hcn1.ra)
       de = im_hms2dec(hcn1.dec)

       gao_data = rsex(hcnpath+'04gao_apj.dat')
       hcn1 = struct_addtags(hcn1,gao_data)
       
       raref = 15.0*im_hms2dec(atlas.ra)
       deref = im_hms2dec(atlas.dec)

       nuclear_raref = 15.0*im_hms2dec(nuclear.ra)
       nuclear_deref = im_hms2dec(nuclear.dec)

       searchrad = 3.0          ; 17.0

    endif

    splog, 'Matching the integrated spectral atlas to the '+hcn_sample+' sample.'
    spherematch, raref, deref, ra, de, searchrad/3600.0, atlas_match, $
      hcn_atlas_match, dist12, maxmatch=1
    nmatch = n_elements(atlas_match)
    splog, 'Found '+string(nmatch,format='(I0)')+' galaxies in common.'
    
    splog, 'Matching the integrated spectral atlas to the '+hcn_sample+' sample.'
    spherematch, nuclear_raref, nuclear_deref, ra, de, searchrad/3600.0, $
      nuclear_match, hcn_nuclear_match, maxmatch=1
    nnuclear_match = n_elements(nuclear_match)
    splog, 'Found '+string(nnuclear_match,format='(I0)')+' galaxies in common.'

    srt = sort(atlas[atlas_match].ra)
    nuclear_srt = sort(nuclear[nuclear_match].ra)

; print the matching statistics       

    splog, '###########################################################################'
    splog, 'Integrated:'
    splog, '###########################################################################'
    for i = 0L, nmatch-1L do print, strtrim(hcn1[hcn_atlas_match[srt[i]]].galaxy,2), $
      strtrim(atlas[atlas_match[srt[i]]].galaxy,2), strtrim(hcn1[hcn_atlas_match[srt[i]]].nedgalaxy,2), $
      strtrim(atlas[atlas_match[srt[i]]].ned_galaxy,2), dist12[srt[i]]*3600.0, $
      strtrim(atlas[atlas_match[srt[i]]].drift_comments,2), format='(A20,A20,A25,A25,F12.3,A45)'

    splog, '###########################################################################'
    splog, 'Nuclear:'
    splog, '###########################################################################'
    for i = 0L, nnuclear_match-1L do print, strtrim(hcn1[hcn_nuclear_match[nuclear_srt[i]]].galaxy,2), $
      strtrim(nuclear[nuclear_match[nuclear_srt[i]]].galaxy,2), strtrim(hcn1[hcn_nuclear_match[nuclear_srt[i]]].nedgalaxy,2), $
      strtrim(nuclear[nuclear_match[nuclear_srt[i]]].ned_galaxy,2), dist12[nuclear_srt[i]]*3600.0, $
      strtrim(nuclear[nuclear_match[nuclear_srt[i]]].drift_comments,2), format='(A20,A20,A25,A25,F12.3,A45)'

; subscript the ATLAS structure and compute reddening-corrected
; quantities
    
    hcn_atlas = atlas[atlas_match[srt]]
    hcn_nuclear = nuclear[nuclear_match[nuclear_srt]]
    hcn = hcn1[hcn_atlas_match[srt]]
;   struct_print, hcn

    hcn_atlas_nodust = iunred_linedust(hcn_atlas,snrcut=dust_snrcut,/silent)
    dust_indx = where(hcn_atlas_nodust.ehbha_err gt 0.0)
;   niceprint, hcn_atlas.galaxy, hcn_atlas_nodust.ebv_hahb, hcn_atlas_nodust.ebv_hahb_err
    
; which of the Baan galaxies are *not* in Gao & Solomon?

;      gao = mrdfits(hcnpath+'gao.ned.fits.gz',1,/silent)
;      jj = cmset_op(strtrim(hcn.nedgalaxy,2),'and',/not2,strtrim(gao.nedgalaxy,2),/index)
;      niceprint, hcn[jj].galaxy, hcn_atlas[jj].galaxy

; initialize the output data structure       
    
    result = {$
      galaxy:                 '-', $ 
      alt_galaxy:             '-', $
      gao_galaxy:             '-', $
      ra_J2000:               '-', $
      dec_J2000:              '-', $
      morph:                  '-', $
      distance:            -999.0, $
      distance_method:        '-', $
      distance_ref:           '-', $
      d25_maj:             -999.0, $
      d25_min:             -999.0, $
      nuclear:                  0, $
      
      gao_distance:        -999.0, $
      gao_lir:             -999.0, $
      gao_lco:             -999.0, $
      gao_lco_err:         -999.0, $
      gao_lhcn:            -999.0, $
      gao_lhcn_err:        -999.0, $
      
      class:                    '',$
      oiii_hb:             -999.0, $
      oiii_hb_err:         -999.0, $
      nii_ha:              -999.0, $
      nii_ha_err:          -999.0, $

      iras_12:             -999.0, $
      iras_25:             -999.0, $
      iras_60:             -999.0, $
      iras_100:            -999.0, $
      iras_12_err:         -999.0, $
      iras_25_err:         -999.0, $
      iras_60_err:         -999.0, $
      iras_100_err:        -999.0, $
      iras_ref:             '...', $
      
      photflag:               '?', $
      h_alpha:             -999.0, $
      h_alpha_err:         -999.0, $
      h_beta:              -999.0, $
      h_beta_err:          -999.0, $
      lha_obs:             -999.0, $
      lha_obs_err:         -999.0, $
      lha_cor:             -999.0, $
      lha_cor_err:         -999.0, $
      lir:                 -999.0, $
      lir_err:             -999.0, $
      l25:                 -999.0, $
      l25_err:             -999.0, $
      l24:                 -999.0, $
      l24_err:             -999.0, $
      iras_60_100:         -999.0, $
      filter_f25_f24:      -999.0, $
      l24_lha_obs:         -999.0, $
      l24_lha_obs_err:     -999.0, $
      
      ehbha:               -999.0, $
      ehbha_err:           -999.0, $
      a_ha:                -999.0, $
      a_ha_err:            -999.0, $
      a_ir:                -999.0, $
      a_ir_err:            -999.0, $
      
      sfr_lir:             -999.0, $
      sfr_lir_err:         -999.0, $
      sfr_lha_obs:         -999.0, $
      sfr_lha_obs_err:     -999.0, $
      sfr_lha_cor:         -999.0, $
      sfr_lha_cor_err:     -999.0, $
      sfr_l24_lha_obs:     -999.0, $
      sfr_l24_lha_obs_err: -999.0}

    result = replicate(result,nmatch)

; basic galaxy properties       
    
    result.galaxy          = strtrim(hcn_atlas.galaxy,2)
    result.gao_galaxy      = strtrim(hcn.gao_galaxy,2)
    result.ra_J2000        = strtrim(hcn_atlas.ra,2)
    result.dec_J2000       = strtrim(hcn_atlas.dec,2)
    result.morph           = strtrim(hcn_atlas.lit_type,2)

    result.nuclear         = hcn_atlas.nuclear
    
    result.iras_12         = hcn_atlas.iras_12
    result.iras_25         = hcn_atlas.iras_25
    result.iras_60         = hcn_atlas.iras_60
    result.iras_100        = hcn_atlas.iras_100
    result.iras_12_err     = hcn_atlas.iras_12_err
    result.iras_25_err     = hcn_atlas.iras_25_err
    result.iras_60_err     = hcn_atlas.iras_60_err
    result.iras_100_err    = hcn_atlas.iras_100_err
    result.iras_ref        = hcn_atlas.iras_100_ref

    result.distance        = hcn_atlas.distance
    result.distance_ref    = hcn_atlas.distance_ref
    result.distance_method = hcn_atlas.distance_method
    result.d25_maj         = hcn_atlas.d25_maj
    result.d25_min         = hcn_atlas.d25_min

    result.alt_galaxy = strcompress(hcn_atlas.alt_galaxy,/remove)
    blank = where(strcompress(result.alt_galaxy,/remove) eq '',nblank)
    if (nblank ne 0L) then result[blank].alt_galaxy = '...'

; quantities from Gao & Solomon (2004); account for the different
; distance scales used 
    
    result.gao_distance = hcn.distance                                                          ; [Mpc]
    result.gao_lir      = alog10(hcn.lir*1D10) + 2.0*alog10(hcn_atlas.distance/hcn.distance)    ; log [L_sun]
    result.gao_lco      = alog10(hcn.lco*1D8) + 2.0*alog10(hcn_atlas.distance/hcn.distance)     ; log [K*km/s/pc^2]
    result.gao_lco_err  = lcoerr
    result.gao_lhcn     = alog10(abs(hcn.lhcn)*1D8) + 2.0*alog10(hcn_atlas.distance/hcn.distance) ; log [K*km/s/pc^2]
    result.gao_lhcn_err = lhcnerr

    neg = where(hcn.lhcn lt 0.0,nneg)
    if (nneg ne 0L) then result[neg].gao_lhcn = -result[neg].gao_lhcn ; upper limits are negative
    neg = where(hcn.lco lt 0.0,nneg)
    if (nneg ne 0L) then result[neg].gao_lco = -result[neg].gao_lco ; upper limits are negative

; classify; IRAS23365+3604 is clearly an AGN, but it fails the S/N cut
; on [NII], so use the [OI]/Ha classification stored in BPT_CLASS;
; also compute the emission-line ratios we're going to need for the
; BPT diagram in HCN_PLOTS; note that IRAS23365+3604 failes the S/N
; cut because [NII] 6584 was not measured

    hcn_class = iclassification(hcn_atlas,ratios=hcn_ratios,$
      snrcut_class=0.0,silent=0,doplot=1,chi2cut=500.0,sigma=1000)
;   hcn_class = iclassification(hcn_atlas,ratios=hcn_ratios,$
;     snrcut_class=class_snrcut,/silent,doplot=1)
stop    
;   niceprint, hcn_class.bpt_galaxy, hcn_class.bpt_nii_mixture_class, hcn_atlas.bpt_class
    result.class = hcn_class.bpt_nii_mixture_class

    iras = where(strmatch(hcn_atlas.galaxy,'*iras23365+3604*',/fold))
    result[iras].class = hcn_class[iras].bpt_class

    lineratio, hcn_atlas, 'NII_6584', 'H_ALPHA', 'OIII_5007', 'H_BETA', $
      nii_ha, nii_ha_err, oiii_hb, oiii_hb_err, index=indx, snrcut=class_snrcut
    result[indx].oiii_hb     = oiii_hb
    result[indx].oiii_hb_err = oiii_hb_err
    result[indx].nii_ha      = nii_ha
    result[indx].nii_ha_err  = nii_ha_err

; store the spectral atlas fluxes and luminosities

    result.photflag    = strtrim(hcn_atlas.drift_photflag,2)                      ; photometric or not?
    result.h_alpha     = hcn_atlas.h_alpha[0]                                     ; [erg/s/cm2]
    result.h_alpha_err = hcn_atlas.h_alpha[1]                                     ; [erg/s/cm2]
    result.h_beta      = hcn_atlas.h_beta[0]                                      ; [erg/s/cm2]
    result.h_beta_err  = hcn_atlas.h_beta[1]                                      ; [erg/s/cm2]
    result.lir         = hcn_atlas.ir_lum                                         ; total L(IR)=L(8-1000) [log L_sun]
    result.lir_err     = hcn_atlas.ir_lum_err                                     ; [log L_sun]
    result.lha_obs     = hcn_atlas.h_alpha_lum[0]                                 ; observed, [log L_sun]
    result.lha_obs_err = hcn_atlas.h_alpha_lum[1]                                 ; observed, [log L_sun]
    result[dust_indx].lha_cor     = hcn_atlas_nodust[dust_indx].h_alpha_lum[0]    ; reddening-corrected, [log L_sun]
    result[dust_indx].lha_cor_err = hcn_atlas_nodust[dust_indx].h_alpha_lum[1]    ; reddening-corrected, [log L_sun]

; compute L(25); recall: 1 Jy = 1D-23 erg/s/cm2/Hz; also compute L(24)
; based on the filter transmission correction function given to me by
; D. Dale (which is a function of the 60/100 flux ratio); once we have
; L(24) we can compute all the optical+mid-IR diagnostics

    indx = where((hcn_atlas.iras_25 gt -900.0) and (hcn_atlas.iras_60 gt -900.0) and $
      (hcn_atlas.iras_100 gt -900.0),nindx)
    if (nindx ne 0L) then begin

       area = 4.0*!dpi*(hcn_atlas[indx].distance*mpc2cm)^2.0
       f25     = hcn_atlas[indx].iras_25*25D4*1D-23*light/(25.0D4)^2.0        ; [erg/s/cm2]
       f25_err = hcn_atlas[indx].iras_25_err*25D4*1D-23*light/(25.0D4)^2.0    ; [erg/s/cm2]

       l25     = f25*area/lsun        ; [L_sun]
       l25_err = f25_err*area/lsun    ; [L_sun]

; filter correction
       
       splog, 'Reading '+hcnpath+'25_24.dat'
       readcol, hcnpath+'25_24.dat', alpha, f25_24, f60_100, format='F,F,F', /silent

       iras_60_100 = hcn_atlas[indx].iras_60/hcn_atlas[indx].iras_100
       filter_f25_f24 = interpol(f25_24,f60_100,iras_60_100)

       l24     = l25/filter_f25_f24
       l24_err = l25_err/filter_f25_f24

       result[indx].l25            = alog10(l25)
       result[indx].l25_err        = l25_err/l25/alog(10.0)
       result[indx].l24            = alog10(l24)
       result[indx].l24_err        = l24_err/l24/alog(10.0)
       result[indx].iras_60_100    = iras_60_100
       result[indx].filter_f25_f24 = filter_f25_f24

    endif

; compute and store some other SFR measures

    indx = where((result.lir gt -900.0),nindx)
    if (nindx ne 0L) then begin
       result.sfr_lir         = result.lir + alog10(lsun) + alog10(kenn_sfr_ir_const)
       result.sfr_lir_err     = result.lir_err
    endif

    result.sfr_lha_obs     = result.lha_obs + alog10(lsun) + alog10(kenn_sfr_ha_const)
    result.sfr_lha_obs_err = result.lha_obs_err
    result[dust_indx].sfr_lha_cor     = result[dust_indx].lha_cor + alog10(lsun) + alog10(kenn_sfr_ha_const)
    result[dust_indx].sfr_lha_cor_err = result[dust_indx].lha_cor_err

; store the reddening values

    result[dust_indx].ehbha     = hcn_atlas_nodust[dust_indx].ehbha
    result[dust_indx].ehbha_err = hcn_atlas_nodust[dust_indx].ehbha_err

    kl_factor = k_lambda(6563.0,/odon)/(k_lambda(4861.0,/odon)-k_lambda(6563.0,/odon))
    result[dust_indx].a_ha     = hcn_atlas_nodust[dust_indx].ehbha*kl_factor
    result[dust_indx].a_ha_err = hcn_atlas_nodust[dust_indx].ehbha_err*kl_factor
;      niceprint, result.aha, -2.5*alog10(10^(result.lha_obs-result.lha_cor)) ; NOTE equivalence!

; store the composite L(Ha)_obs + L(24) luminosity and the
; optical+mid-IR SFR estimator proposed by Calzetti et al. and
; Kennicutt et al.; for simplicity, we first "undo" the logarithms;
; also compute the reddening implied by the L(24)/L(Ha)_obs ratio

    indx = where((result.lha_obs gt -900.0) and (result.l24 gt -900.0),nindx)
    if (nindx ne 0L) then begin

       lha_obs     = lsun*10.0^result[indx].lha_obs 
       lha_obs_err = result[indx].lha_obs_err*lha_obs*alog(10.0)
       l24         = lsun*10.0^result[indx].l24
       l24_err     = result[indx].l24_err*l24*alog(10.0)

       l24_lha_obs     = lha_obs+kenn_factor*l24
       l24_lha_obs_err = sqrt((lha_obs_err)^2.0+(kenn_factor*l24_err)^2.0+(kenn_factor_err*l24)^2.0)

       result[indx].l24_lha_obs     = alog10(l24_lha_obs/lsun)
       result[indx].l24_lha_obs_err = l24_lha_obs_err/l24_lha_obs/alog(10.0)

       sfr_l24_lha_obs     = kenn_sfr_ha_const*l24_lha_obs
       sfr_l24_lha_obs_err = kenn_sfr_ha_const*l24_lha_obs_err
       
       result[indx].sfr_l24_lha_obs = alog10(sfr_l24_lha_obs)
       result[indx].sfr_l24_lha_obs_err = sfr_l24_lha_obs_err/sfr_l24_lha_obs/alog(10.0)

       argument = 1.0+kenn_factor*l24/lha_obs
       argument_err = sqrt((kenn_factor*(l24_err/lha_obs-l24*lha_obs_err/lha_obs^2.0))^2.0 + $
         (kenn_factor_err*(l24/lha_obs))^2.0)
       
       result[indx].a_ir = 2.5*alog10(argument)
       result[indx].a_ir_err = 2.5*argument_err/argument/alog(10.0)

    endif

; write out

    if keyword_set(write) then begin

;      struct_print, struct_trimtags(result,select=['galaxy','class','nuclear'])
;      struct_print, struct_trimtags(result,select=['galaxy','lha*','*l24*'])

       if file_test(hcnpath+'atlas_hcn_sample.fits.gz',/regular) then begin
          splog, 'Press Y to overwrite existing file '+hcnpath+'atlas_hcn_sample.fits.gz'
          cc = get_kbrd(1)
          if (strmatch(cc,'y',/fold) eq 0B) then return
       endif

       splog, 'Writing '+hcnpath+'atlas_hcn_sample.fits.gz'
       mwrfits, result, hcnpath+'atlas_hcn_sample.fits', /create
       spawn, 'gzip -f '+hcnpath+'atlas_hcn_sample.fits', /sh

    endif
       
; ###########################################################################    
; make thumbnail DSS visualizations
; ###########################################################################    

    if keyword_set(thumbs) then begin
       
       dsspath = atlas_path(/dss)
       dssfits = dsspath+strlowcase(strtrim(hcn_atlas.galaxy,2))+'.fits.gz'

       if keyword_set(write) then begin
          dfpsplot, hcnpath+'atlas_hcn_sample.ps', /square, /color
          postthick = 5.0
       endif else begin
          im_window, 0, xratio=0.4, /square
          postthick = 2.0
       endelse

       for i = 0L, nmatch-1L do begin
          
          dssimage = readfits(dssfits[i],hdss,/silent)
          gsssextast, hdss, astr

          imsize = size(dssimage,/dimension)
          xsize = imsize[0] & xcen = xsize/2.0 & ysize = imsize[1] & ycen = ysize/2.0

          xpixscale = (astr.pltscl*1E-3*astr.xsz)/60 ; [arcmin/pixel]
          ypixscale = (astr.pltscl*1E-3*astr.ysz)/60 ; [arcmin/pixel]

          xaxis = (findgen(xsize)-xcen)*xpixscale ; centered on the image [arcsec]
          yaxis = (findgen(ysize)-ycen)*ypixscale ; centered on the image [arcsec]

          img = logscl(dssimage,exponent=1.0,negative=keyword_set(write),omin=35,omax=255)
          
          plotimage, img, /preserve_aspect, position=pos, /normal, imgxrange=minmax(xaxis), $
            imgyrange=minmax(yaxis), charsize=1.8, charthick=postthick, xthick=postthick, $
            ythick=postthick, xtitle=textoidl('\Delta\alpha [arcmin]'), $
            ytitle=textoidl('\Delta\delta [arcmin]'), title=name ;, color=djs_icolor('white')

          gsssadxy, astr, 15.0*im_hms2dec(hcn_atlas[i].ra), im_hms2dec(hcn_atlas[i].dec), xrc3, yrc3
          xrc3 = (xrc3 - xcen)*xpixscale & yrc3 = (yrc3 - ycen)*ypixscale

          if hcn_atlas[i].drift_agnflag then broad = 'Type 1 AGN' else broad = ''
          if (result[i].gao_lhcn lt 0.0) then limit = '>' else limit = ''
          if (result[i].lir lt -900.0) then lir = '?' else lir = strtrim(string(result[i].lir,format='(F5.2)'),2)

          if (strtrim(result[i].galaxy,2) ne strtrim(result[i].gao_galaxy,2)) then $
            gal = strtrim(result[i].galaxy,2)+'='+strtrim(result[i].gao_galaxy,2) else $
              gal = strtrim(result[i].galaxy,2)
          label = [gal,$
            'L(IR) = '+lir,$
            'L(HCN) = '+limit+strtrim(string(result[i].gao_lhcn,format='(F5.2)'),2),$
            'L(CO) = '+strtrim(string(result[i].gao_lco,format='(F5.2)'),2)]
          legend, label, /left, /top, box=0, charthick=postthick, charsize=1.5
          legend, broad, /right, /top, box=0, charthick=postthick, charsize=1.5
          tvcircle, 1.1/2.0, 0.0, 0.0, /data, thick=3.0
;         tvcircle, 0.5, !x.crange[0], !y.crange[0], /data, thick=3.0
          
          if (not keyword_set(write)) then cc = get_kbrd(1)

       endfor
       
       if keyword_set(write) then begin
          dfpsclose
          spawn, ['gzip -f '+hcnpath+'atlas_hcn_sample.ps'], /sh
       endif

    endif
       
return
end

;; ---------------------------------------------------------------------------
;       n3690_ne = where(strtrim(atlas.galaxy,2) eq 'NGC3690NE')
;       n3690_sw = where(strtrim(atlas.galaxy,2) eq 'NGC3690SW')
;       n3690 = where(strtrim(atlas.galaxy,2) eq 'NGC3690')
;;      print, atlas[n3690].h_alpha_lum[0], atlas[n3690_ne].h_alpha_lum[0], atlas[n3690_sw].h_alpha_lum[0]
;
;       atlas[n3690_ne].ra           = atlas[n3690].ra
;       atlas[n3690_ne].dec          = atlas[n3690].dec
;       atlas[n3690_ne].iras_12      = atlas[n3690].iras_12
;       atlas[n3690_ne].iras_25      = atlas[n3690].iras_25
;       atlas[n3690_ne].iras_60      = atlas[n3690].iras_60
;       atlas[n3690_ne].iras_100     = atlas[n3690].iras_100
;       atlas[n3690_ne].iras_12_err  = atlas[n3690].iras_12_err
;       atlas[n3690_ne].iras_25_err  = atlas[n3690].iras_25_err
;       atlas[n3690_ne].iras_60_err  = atlas[n3690].iras_60_err
;       atlas[n3690_ne].iras_100_err = atlas[n3690].iras_100_err
;       atlas[n3690_ne].ir_lum       = atlas[n3690].ir_lum
;       atlas[n3690_ne].ir_lum_err   = atlas[n3690].ir_lum_err
;
;       atlas[n3690].ra    = '-99:00:00' & atlas[n3690].dec    = '-99:00:00'
;       atlas[n3690_sw].ra = '-99:00:00' & atlas[n3690_sw].dec = '-99:00:00'
;;      atlas[n3690_ne].ra = '-99:00:00' & atlas[n3690_ne].dec = '-99:00:00'
;
;; ---------------------------------------------------------------------------
       
; mess with the coordinates for NGC3690/ARP299 and IC1623 to get the
; right matching objects; see comments for 04GAO_APJ.DAT; a search
; radius of 17" ensures that the matching is not contaminated by
; NGC1143 or ARP118=NGC1143/NGC1144; precess the coordinates in Gao &
; Solomon

;         hcn1 = rsex(path+'04gao_apj.dat')
;         hcn1 = hcn1[sort(hcn1.ra)]
;         ra = 15.0*im_hms2dec(hcn1.ra_b1950)
;         de = im_hms2dec(hcn1.dec_b1950)
;         precess, ra, de, 1950.0, 2000.0

;         indx = speclinefit_locate(atlas,'NGC3690')
;         atlas[indx[[1,2]]].ra = '-99:00:00'
;         atlas[indx[[1,2]]].dec = '-99:00:00'
;
;         indx = speclinefit_locate(atlas,'IC1623')
;         atlas[indx[[1,2]]].ra = '-99:00:00'
;         atlas[indx[[1,2]]].dec = '-99:00:00'

