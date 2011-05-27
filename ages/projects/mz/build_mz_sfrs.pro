function init_sfrs, ngal, sdss=sdss
; initialize the data table
    if keyword_set(sdss) then begin
       pre = {object_position: 0L, vagc_object_position: 0L}
    endif else begin
       pre = {ages_id: 0L}
    endelse
    sfrs = struct_addtags(pre,{$
      ra:            0.0D,$
      dec:           0.0D,$
      z:           -999.0,$
      agn:             -1,$     ; 1=yes, 0=no, -1=no classification
      ohsample:         0,$     ; part of the metallicity sample?
      weight:         1.0,$
      
      mips:             0,$     ; MIPS source? 1=detection, 0=limit
      f24:         -999.0,$     ; 24-micron flux and error
      f24_err:     -999.0,$ 
      sfr24:       -999.0,$     ; 24-micron SFR
      sfr24_err:   -999.0,$ 

      mass_50:     -999.0,$
      mass_avg:    -999.0,$
      mass_err:    -999.0,$
      massmpa_50:  -999.0,$
      massmpa_avg: -999.0,$
      massmpa_err: -999.0,$
      
      sfr_50:  -999.0,$         ; from iSEDfit; instantaneous
      sfr_avg: -999.0,$
      sfr_err: -999.0,$
        
      sfr100_50:  -999.0,$      ; from iSEDfit; 100 Myr
      sfr100_avg: -999.0,$
      sfr100_err: -999.0,$
        
      sfrmpa_flag:     0,$      ; 0=good
      sfrmpa_50:  -999.0,$      ; from Brinchmann+
      sfrmpa_avg: -999.0,$
      sfrmpa_err: -999.0,$

      sfrmpa_dr4_50:  -999.0,$  ; from Brinchmann+
      sfrmpa_dr4_avg: -999.0,$
      sfrmpa_dr4_err: -999.0})
    sfrs = replicate(sfrs,ngal)
return, sfrs
end

pro build_mz_sfrs, sfrs, sdss=sdss, clobber=clobber
; jm11may14ucsd - build the AGES and SDSS SFR measurements table;
; identify objects that are part of the metallicity sample

    mzpath = mz_path()
    catpath = ages_path(/catalogs)

    parent = read_mz_sample(sdss=sdss,/parent)
    mass = read_mz_sample(sdss=sdss,/mass)
    oh = read_mz_sample(sdss=sdss,/mzhii_log12oh)
    ispec = read_mz_sample(sdss=sdss,/mzhii_ispec)
    ngal = n_elements(parent)

    sfrs = init_sfrs(ngal,sdss=sdss)
    sfrs.z = parent.z
    sfrs.ra = parent.ra
    sfrs.dec = parent.dec
    sfrs.weight = parent.final_weight

    sfrs.mass_50 = mass.mass_50
    sfrs.mass_avg = mass.mass_avg
    sfrs.mass_err = mass.mass_err

    sfrs.sfr_50 = mass.sfr_50
    sfrs.sfr_avg = mass.sfr_avg
    sfrs.sfr_err = mass.sfr_err

    sfrs.sfr100_50 = mass.sfr100_50
    sfrs.sfr100_avg = mass.sfr100_avg
    sfrs.sfr100_err = mass.sfr100_err

    if keyword_set(sdss) then begin ; SDSS
       splog, file=mzpath+'build_mz_sfrs.sdss.log'
       splog, '#########################'
       splog, 'Building SDSS SFRs'

       sfrs.object_position = parent.object_position
       sfrs.vagc_object_position = parent.vagc_object_position

; identify AGN
       vagc = mz_get_vagc(sample=sample,letter=letter,poststr=poststr)
       allispec = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/ispec)
       allispec = allispec[parent.vagc_object_position]

       class = iclassification(allispec,snrcut=3.0,ratios=iratios)
       sfrs.agn = strmatch(iratios.final_class,'*AGN*',/fold)
       sfrs[where(strtrim(iratios.final_class,2) eq 'SF')].agn = 0
       sfrs[where(strtrim(iratios.final_class,2) eq 'UNKNOWN')].agn = -1

; in the metallicity sample?       
       spherematch, parent.ra, parent.dec, ispec.ra, ispec.dec, 1.0/3600.0, m1, m2
       sfrs[m1].ohsample = 1

; get the MPA SFRs, including, for comparison, the DR4 estimates 
       sfrs.massmpa_50 = parent.mass_median-0.07 ; Kroupa-->Chabrier
       sfrs.massmpa_avg = parent.mass_avg-0.07
       sfrs.massmpa_err = (parent.mass_p84-parent.mass_p16)/2.0

       sfrs.sfrmpa_50 = parent.sfr_median-0.07 ; Kroupa-->Chabrier
       sfrs.sfrmpa_avg = parent.sfr_avg-0.07
       sfrs.sfrmpa_err = (parent.sfr_p84-parent.sfr_p16)/2.0

       dr4 = mrdfits(sdss_path(/mpa_dr4)+'mpamassoh_dr4_v5_1b.fits.gz',1)
       gd = where(dr4.ra gt 0.0)
       spherematch, parent.ra, parent.dec, dr4[gd].ra, dr4[gd].dec, $
         1.0/3600.0, m1, m2
       sfrs[m1].sfrmpa_dr4_50 = dr4[gd[m2]].sfr_median-0.07 ; Kroupa-->Chabrier
       sfrs[m1].sfrmpa_dr4_avg = dr4[gd[m2]].sfr_avg-0.07
       sfrs[m1].sfrmpa_dr4_err = (dr4[gd[m2]].sfr_p84-dr4[gd[m2]].sfr_p16)/2.0
       
       im_mwrfits, sfrs, mzpath+'sdss_parent_sfrs.fits', clobber=clobber

    endif else begin            ; AGES
       splog, file=mzpath+'build_mz_emline.sfrs.log'
       splog, '#########################'
       splog, 'Building AGES SFRs'

; identify AGN
       ppxf = read_ages(/ppxf)
       match, parent.ages_id, ppxf.ages_id, m1, m2
       ppxf = ppxf[m2]

       class = iclassification(ppxf,snrcut=3.0,ratios=iratios)
       sfrs.agn = strmatch(iratios.final_class,'*AGN*',/fold)
       sfrs[where(strtrim(iratios.final_class,2) eq 'SF')].agn = 0
       sfrs[where(strtrim(iratios.final_class,2) eq 'UNKNOWN')].agn = -1

; in the metallicity sample?
       sfrs.ages_id = parent.ages_id
       match, sfrs.ages_id, ispec.ages_id, m1, m2
       sfrs[m1].ohsample = 1

; get the SFR using the prescription from Rieke+09       
       mips = where(parent.phot_mips24 gt -900.0,nmips,comp=nomips)
       sfrs[mips].mips = 1

       sfrs.f24 = 0.27 ; 1-sigma limiting flux [mJy]
       sfrs[mips].f24 = parent[mips].phot_mips24
       sfrs[mips].f24_err = parent[mips].phot_mips24_err
;      lir = im_f24tolir(sfrs.z,sfrs.f24,/rieke)

       sfrs.sfr24 = sfr_09rieke(sfrs.z,sfrs.f24*1D-3)-0.07 ; log (Msun/yr); Kroupa02-->Chabrier
       sfrs.sfr24_err = 0.3 ; assume a factor of 2 [dex]

;; the statistical errors are tiny!       
;       nmonte = 100
;       for ii = 0, nmips-1 do begin
;          f24monte = sfrs[mips[ii]].f24 + randomn(seed,nmonte)*sfrs[mips[ii]].f24_err
;          sfr24monte = sfr_09rieke(fltarr(nmonte)+sfrs[mips[ii]].z,f24monte*1D-3)
;          sfrs[mips[ii]].sfr24_err = djsig(sfr24monte)
;       endfor

       im_mwrfits, sfrs, mzpath+'ages_parent_sfrs.fits', clobber=clobber
    endelse
    splog, /close
    
return
end
