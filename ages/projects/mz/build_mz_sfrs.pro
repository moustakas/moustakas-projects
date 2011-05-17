function init_sfrs, ngal, sdss=sdss
; initialize the data table
    if keyword_set(sdss) then begin
       pre = {object_position: 0L}
    endif else begin
       pre = {$
         ages_id:     0L,$
         mips:         0,$ ; MIPS source?
         f24:     -999.0,$ ; 24-micron flux and error
         f24_err: -999.0,$ 
         sfr24:   -999.0}  ; 24-micron SFR
    endelse
    sfrs = struct_addtags(pre,{$
      ra:           0.0D,$
      dec:          0.0D,$
      z:          -999.0,$
      agn:             0,$
      mass:       -999.0,$
      mass_err:   -999.0,$
      sfr100:     -999.0,$ ; from iSEDfit
      sfr100_err: -999.0})
    sfrs = replicate(sfrs,ngal)
return, sfrs
end

pro build_mz_sfrs, sdss=sdss, clobber=clobber, doplot=doplot
; jm11may14ucsd - build the AGES and SDSS SFR measurements table;
; identify objects that are part of the metallicity sample

    mzpath = mz_path()
    catpath = ages_path(/catalogs)

    parent = read_mz_sample(sdss=sdss,/parent)
    mass = read_mz_sample(sdss=sdss,/mass)
    oh = read_mz_sample(sdss=sdss,/mzhii_log12oh)
    ngal = n_elements(parent)

    sfrs = init_sfrs(ngal,sdss=sdss)
    sfrs.z = parent.z
    sfrs.ra = parent.ra
    sfrs.dec = parent.dec

    sfrs.mass = mass.mass_50
    sfrs.mass_err = mass.mass_err
    sfrs.sfr100 = mass.sfr100_avg
;   sfrs.sfr100 = mass.sfr100_50
    sfrs.sfr100_err = mass.sfr100_err
    
; identify AGN
;   sfrs.agn = sample.x_match or ages_irac_agn(sample)

    if keyword_set(sdss) then begin ; SDSS
       splog, file=mzpath+'build_mz_sfrs.sdss.log'
       splog, '#########################'
       splog, 'Building SDSS SFRs'

       sfrs.object_position = parent.object_position

       im_mwrfits, mz_parent, mzpath+'sdss_mz_sfrs.fits', clobber=clobber

    endif else begin            ; AGES
       splog, file=mzpath+'build_mz_emline.sfrs.log'
       splog, '#########################'
       splog, 'Building AGES SFRs'

       sfrs.ages_id = parent.ages_id

; compute L(IR) and L(24) using the observed 24-micron flux densities
; and the average of the Chary & Elbaz (2001) and Dale & Helou (2002)
; infrared templates
       mips = where(parent.phot_mips24 gt -900.0,nmips)
       sfrs[mips].mips = 1
       sfrs[mips].f24 = parent[mips].phot_mips24
       sfrs[mips].f24_err = parent[mips].phot_mips24_err
       sfrs[mips].sfr24 = sfr_09rieke(sfrs[mips].z,sfrs[mips].f24*1D-3)

;      plot, sfrs[mips].z, sfrs[mips].sfr24, ps=6, ysty=3
       plot, sfrs[mips].mass, alog10(sfrs[mips].sfr100)-sfrs[mips].sfr24, $
         ps=6, xr=[8,12], yr=[-2,2], sym=0.1
       ss = im_stats(alog10(sfrs[mips].sfr100)-sfrs[mips].sfr24,/ver) 

stop       
       
       im_mwrfits, mz_parent, mzpath+'ages_parent_sfrs.fits', clobber=clobber
    endelse
    splog, /close
    
return
end
