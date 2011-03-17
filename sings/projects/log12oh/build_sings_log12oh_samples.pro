pro build_sings_log12oh_samples, hii=hii, sdss=sdss, clobber=clobber
; jm10mar06ucsd - build the requisite samples for the SINGS/log12oh
; sample; rewritten version of WRITE_SINGS_LOG12OH_SAMPLES 

    version = sings_log12oh_version()
    outpath = sings_path(/projects)+'log12oh/'

; ###########################################################################
; HII regions
    if keyword_set(hii) then begin

; read and cross-match the list of SINGS galaxies against the
; HII-region database
       sings = sings_read_info()
       bighii = read_hiiregions()

       spherematch, 15.0*im_hms2dec(sings.ra), im_hms2dec(sings.dec), $
         15.0*im_hms2dec(bighii.galaxy_ra), im_hms2dec(bighii.galaxy_dec), $
         3.0/3600.0, sings_match, bighii_match, maxmatch=0 ; return all matches
       hii = bighii[bighii_match]

; require H-beta, [OII], and [OIII]
       good = where((hii.zstrong_p gt -900.0) and (hii.zstrong_r23 gt -900.0) and $
;        (hii.hii_rc3_rr25 gt 0.0) and $
         (hii.hii_rc3_rr25 lt 1.0),ngood)
       hii = hii[good]
       
;      doit = match_string(sings.ned_galaxy,bighii.ned_galaxy,index=index,/exact,/silent)
;      good = where(index ne -1L,ngood)
;      hii = bighii[index[good]]

       uindx = uniq(strtrim(hii.ned_galaxy,2),sort(strtrim(hii.ned_galaxy)))
       nmatch = n_elements(uindx)
       splog, 'Identified '+string(nmatch,format='(I0)')+' galaxies with fluxes for '+$
         string(ngood,format='(I0)')+' HII regions'
       
;      niceprint, ((sings[sings_match])[uindx]).galaxy, ((sings[sings_match])[uindx]).ned_galaxy, $
;        hii[uindx].ned_galaxy, hii[uindx].hii_galaxy

       ref = strtrim(hii.texref,2)
       uref = ref[uniq(ref,sort(ref))]
       splog, 'Unique references: '+string(n_elements(uref),format='(I0)')
       niceprint, uref
       
       hiifile = 'sings_hiiregions_'+version+'.fits'
       im_mwrfits, hii, hiifile, clobber=clobber

       return
    endif

; ###########################################################################
; read and write out an SDSS emission-line sample
    if keyword_set(sdss) then begin
       sample = 'dr72'
       letter = 'bsafe'
       poststr = '25' ; 0.01<z<0.25
       suffix = sample+letter+poststr

       lsspath = getenv('LSS_REDUX')+'/'+sample+'/'+letter+'/'+poststr+'/'

       postfile = lsspath+'post_catalog.'+suffix+'.fits'
       splog, 'Reading '+postfile
       post = mrdfits(postfile,1)

       ispecfile = lsspath+'ispecline.'+suffix+'.fits.gz'
       splog, 'Reading '+ispecfile
       ispec = mrdfits(ispecfile,1)

       kcorrfile = getenv('VAGC_REDUX')+'/kcorrect/kcorrect.nearest.model.fits.gz'
       splog, 'Reading '+kcorrfile
       kcorr = mrdfits(kcorrfile,1,rows=post.object_position)

; require a match to the VAGC with good stellar masses
       good = where((kcorr.mass gt -900.0) and (ispec.plateid gt -900),ngood)
       splog, 'SDSS sample: ', ngood

       out = struct_addtags(kcorr[good],struct_trimtags(ispec[good],$
         select=['linename',strtrim(ispec[0].linename,2)+'*']))
       outfile = outpath+'sdss_log12oh_'+version+'.fits'
       im_mwrfits, out, outfile, clobber=clobber

       return
    endif
    
; ###########################################################################
; SINGS data    
    
    snrcut1 = 2.0
    nmonte = 500

    class = mrdfits(outpath+'sings_class_'+version+'.fits',1,/silent)
    allnuclear = read_sings_gandalf(/nuclear)
    alldrift20 = read_sings_gandalf(/drift20)
    alldrift56 = read_sings_gandalf(/drift56)

; -------------------------
    splog, '##################################################'
    splog, 'DRIFT56:'
    splog, '##################################################'

    match, strtrim(alldrift56.ned_galaxy,2), strtrim(class.ned_galaxy,2), m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    nagn = fix(total(strtrim(class[m2].drift56_class,2) eq 'AGN'))

    drift56_keep = where((strtrim(class[m2].drift56_class,2) ne 'AGN') and $
      (alldrift56.h_alpha[0]/alldrift56.h_alpha[1] gt snrcut1) and $
      (alldrift56.h_beta[0]/alldrift56.h_beta[1] gt snrcut1) and $
      (alldrift56.oii_3727[0]/alldrift56.oii_3727[1] gt snrcut1) and $
      (alldrift56.oiii_5007[0]/alldrift56.oiii_5007[1] gt snrcut1),$
      ndrift56_keep,comp=drift56_fail,ncomp=ndrift56_fail)

; compute the reddening; NGC5474 has an unphysical decrement
    drift56 = alldrift56[drift56_keep]
    drift56_crap = alldrift56[drift56_fail]

    splog, 'Computing the reddening.'
    drift56_nodust = iunred_linedust(drift56,snrcut=1.0,use_hahb=2.86,$
      allow_r_min=1,silent=1,nopropagate=0)
    drift56_dust_good = where(drift56_nodust.ebv_hahb_err gt 0.0,ndrift56_dust_good,$
      comp=drift56_dust_crap,ncomp=ndrift56_dust_crap)

    if (ndrift56_dust_crap ne 0) then begin
       splog, 'Unphysical Balmer decrements:'
       for ii = 0L, ndrift56_dust_crap-1L do print, '   '+$
         strtrim(drift56_nodust[drift56_dust_crap[ii]].galaxy,2), $
         drift56_nodust[drift56_dust_crap[ii]].hahb, $
         drift56_nodust[drift56_dust_crap[ii]].hahb_err
       drift56_crap = struct_append(drift56_crap,drift56[drift56_dust_crap])
       drift56_crap = drift56_crap[sort(drift56_crap.galaxy)]
    endif

; index the objects that made the reddening cut    
    if (ndrift56_dust_good ne 0L) then begin
       drift56 = drift56[drift56_dust_good]
       drift56_nodust = drift56_nodust[drift56_dust_good]
    endif
    ndrift56_good = n_elements(drift56)
    ndrift56_crap = n_elements(drift56_crap)
    
    splog, string(ndrift56_good,format='(I0)')+'/'+string(n_elements(alldrift56),format='(I0)')+$
      ' ('+ string(100.0*ndrift56_good/float(n_elements(alldrift56)),format='(F4.1)')+$
      '%) galaxies made the S/N and reddening cuts, excluding '+string(nagn,format='(I0)')+' AGN'

; compute oxygen abundances
    splog, 'Computing oxygen abundances.'
    abund = im_abundance(drift56_nodust,snrcut=0.0,nmonte=nmonte,$
      silent=0,/justflux)
    drift56_nodust = struct_addtags(drift56_nodust,abund)

    niceprint, drift56_nodust.galaxy, drift56_nodust.ebv_hahb, $
      abund.zstrong_r23, abund.zstrong_r23_err, abund.zstrong_o32, $
      abund.zstrong_o32_err
    
;   splog, 'DRIFT56: Generating QA postscript files.'
;   sings_display_spectrum, drift56, labeltype=3, /drift56, /post, $
;     psname=outpath+'spectra_drift56_keep_'+version+'.ps'
;   sings_display_spectrum, drift56_crap, labeltype=3, /drift56, /post, $
;     psname=outpath+'spectra_drift56_crap_'+version+'.ps'

    outfile = outpath+'sings_log12oh_drift56_speclinefit_'+version+'.fits'
    outfile_nodust = repstr(outfile,'.fits','_nodust.fits')
    im_mwrfits, drift56, outfile, clobber=clobber
    im_mwrfits, drift56_nodust, outfile_nodust, clobber=clobber

; -------------------------
    splog, '##################################################'
    splog, 'DRIFT20:'
    splog, '##################################################'

    match, strtrim(alldrift20.ned_galaxy,2), strtrim(class.ned_galaxy,2), m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    nagn = fix(total(strtrim(class[m2].drift20_class,2) eq 'AGN'))
    
; possibly also reject NGC7331    
    drift20_keep = where((strtrim(class[m2].drift20_class,2) ne 'AGN') and $
      (alldrift20.h_alpha[0]/alldrift20.h_alpha[1] gt snrcut1) and $
      (alldrift20.h_beta[0]/alldrift20.h_beta[1] gt snrcut1) and $
      (alldrift20.oii_3727[0]/alldrift20.oii_3727[1] gt snrcut1) and $
      (alldrift20.oiii_5007[0]/alldrift20.oiii_5007[1] gt snrcut1),$
      ndrift20_keep,comp=drift20_fail,ncomp=ndrift20_fail)

; compute the reddening; NGC2915 has an unphysical Balmer decrement;
; use Ha/Hb=3.1 for AGN and 2.86 otherwise
    drift20 = alldrift20[drift20_keep]
    drift20_crap = alldrift20[drift20_fail]

; NGC2915 has an unphysical decrement
    splog, 'Computing the reddening.'
    drift20_nodust = iunred_linedust(drift20,snrcut=1.0,use_hahb=2.86,$
      allow_r_min=1,silent=1,nopropagate=0)
    drift20_dust_good = where(drift20_nodust.ebv_hahb_err gt 0.0,ndrift20_dust_good,$
      comp=drift20_dust_crap,ncomp=ndrift20_dust_crap)

    if (ndrift20_dust_crap ne 0) then begin
       splog, 'Unphysical Balmer decrements:'
       for ii = 0, ndrift20_crap-1 do print, '   '+$
         strtrim(drift20_nodust[drift20_dust_crap[ii]].galaxy,2), $
         drift20_nodust[drift20_dust_crap[ii]].hahb, $
         drift20_nodust[drift20_dust_crap[ii]].hahb_err, $
         drift20_nodust[drift20_dust_crap[ii]].ebv_hahb, $
         drift20_nodust[drift20_dust_crap[ii]].ebv_hahb_err
       drift20_crap = struct_append(drift20_crap,drift20[drift20_dust_crap])
       drift20_crap = drift20_crap[sort(drift20_crap.galaxy)]
    endif

; index the objects that made the reddening cut    
    if (ndrift20_dust_good ne 0) then begin
       drift20 = drift20[drift20_dust_good]
       drift20_nodust = drift20_nodust[drift20_dust_good]
    endif
    ndrift20_good = n_elements(drift20)
    ndrift20_crap = n_elements(drift20_crap)
    
    splog, string(ndrift20_good,format='(I0)')+'/'+string(n_elements(alldrift20),format='(I0)')+$
      ' ('+ string(100.0*ndrift20_good/float(n_elements(alldrift20)),format='(F4.1)')+$
      '%) galaxies made the S/N and reddening cuts, excluding '+string(nagn,format='(I0)')+' AGN'

; compute oxygen abundances
    splog, 'Computing oxygen abundances.'
    abund = im_abundance(drift20_nodust,snrcut=0.0,$
      nmonte=nmonte,silent=0,/justflux)
    drift20_nodust = struct_addtags(drift20_nodust,abund)

    niceprint, drift20_nodust.galaxy, drift20_nodust.ebv_hahb, $
      abund.zstrong_r23, abund.zstrong_r23_err, abund.zstrong_o32, $
      abund.zstrong_o32_err

;   splog, 'DRIFT20: Generating QA postscript files.'
;   sings_display_spectrum, drift20, labeltype=3, /drift20, /post, $
;     psname=outpath+'spectra_drift20_keep_'+version+'.ps'
;   sings_display_spectrum, drift20_crap, labeltype=3, /drift20, /post, $
;     psname=outpath+'spectra_drift20_crap_'+version+'.ps'
    
    outfile = outpath+'sings_log12oh_drift20_speclinefit_'+version+'.fits'
    outfile_nodust = repstr(outfile,'.fits','_nodust.fits')
    im_mwrfits, drift20, outfile, clobber=clobber
    im_mwrfits, drift20_nodust, outfile_nodust, clobber=clobber

; -------------------------
    splog, '##################################################'
    splog, 'NUCLEAR:'
    splog, '##################################################'

    match, strtrim(allnuclear.ned_galaxy,2), strtrim(class.ned_galaxy,2), m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    nagn = fix(total(strtrim(class[m2].nuclear_class,2) eq 'AGN'))
    
    nuclear_keep = where((strtrim(class[m2].nuclear_class,2) ne 'AGN') and $
      (allnuclear.h_alpha[0]/allnuclear.h_alpha[1] gt snrcut1) and $
      (allnuclear.h_beta[0]/allnuclear.h_beta[1] gt snrcut1) and $
      (allnuclear.oii_3727[0]/allnuclear.oii_3727[1] gt snrcut1) and $
      (allnuclear.oiii_5007[0]/allnuclear.oiii_5007[1] gt snrcut1),$
      nnuclear_keep,comp=nuclear_fail,ncomp=nnuclear_fail)

; compute the reddening; NGC0925, NGC3049, and NGC3351 have unphysical
; Balmer decrements, while the decrement for NGC3049 requires
; ALLOW_R_MIN=1; use Ha/Hb=3.1 for AGN and 2.86 otherwise
    nuclear = allnuclear[nuclear_keep]
    nuclear_crap = allnuclear[nuclear_fail]

; NGC5474 has a formally negative reddening value because H-beta has
; been over-subtracted, I think; but keep it because it has strong
; lines     
    splog, 'Computing the reddening.'
    nuclear_nodust = iunred_linedust(nuclear,snrcut=0.0,use_hahb=2.86,$
      allow_r_min=1,silent=1,nopropagate=0)
    nuclear_dust_good = where(nuclear_nodust.ebv_hahb_err gt 0.0,nnuclear_dust_good,$
      comp=nuclear_dust_crap,ncomp=nnuclear_dust_crap)
    if (nnuclear_dust_crap ne 0) then begin
       splog, 'Unphysical Balmer decrements:'
       for ii = 0L, nnuclear_dust_crap-1L do print, '   '+$
         strtrim(nuclear_nodust[nuclear_dust_crap[ii]].galaxy,2), $
         nuclear_nodust[nuclear_dust_crap[ii]].hahb, $
         nuclear_nodust[nuclear_dust_crap[ii]].hahb_err
       nuclear_crap = struct_append(nuclear_crap,nuclear[nuclear_dust_crap])
       nuclear_crap = nuclear_crap[sort(nuclear_crap.galaxy)]
    endif

; index the objects that made the reddening cut    
    if (nnuclear_dust_good ne 0) then begin
       nuclear = nuclear[nuclear_dust_good]
       nuclear_nodust = nuclear_nodust[nuclear_dust_good]
    endif
    nnuclear_good = n_elements(nuclear)
    nnuclear_crap = n_elements(nuclear_crap)
    
    splog, string(nnuclear_good,format='(I0)')+'/'+string(n_elements(allnuclear),format='(I0)')+$
      ' ('+ string(100.0*nnuclear_good/float(n_elements(allnuclear)),format='(F4.1)')+$
      '%) galaxies made the S/N and reddening cuts, excluding '+string(nagn,format='(I0)')+' AGN'

; compute oxygen abundances
    splog, 'Computing oxygen abundances.'
    abund = im_abundance(nuclear_nodust,snrcut=0.0,$
      nmonte=nmonte,silent=0,/justflux)
    nuclear_nodust = struct_addtags(nuclear_nodust,abund)

    niceprint, nuclear_nodust.galaxy, nuclear_nodust.ebv_hahb, $
      abund.zstrong_r23, abund.zstrong_r23_err, abund.zstrong_o32, $
      abund.zstrong_o32_err

;   splog, 'NUCLEAR: Generating QA postscript files.'
;   sings_display_spectrum, nuclear, labeltype=3, /nuclear, /post, $
;     psname=outpath+'spectra_nuclear_keep_'+version+'.ps'
;   sings_display_spectrum, nuclear_crap, labeltype=3, /nuclear, /post, $
;     psname=outpath+'spectra_nuclear_crap_'+version+'.ps'

    outfile = outpath+'sings_log12oh_nuclear_speclinefit_'+version+'.fits'
    outfile_nodust = repstr(outfile,'.fits','_nodust.fits')
    im_mwrfits, nuclear, outfile, clobber=clobber
    im_mwrfits, nuclear_nodust, outfile_nodust, clobber=clobber

return
end
