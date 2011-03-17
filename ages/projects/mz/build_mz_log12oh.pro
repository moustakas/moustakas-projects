function get_log12oh, abund, r23branch=r23branch, debug=debug, $
  justew=justew, justflux=justflux
;   abund = struct_trimtags(abund,except=['*pt05*',$
;     '*pettini*','*zkh*','*dens*','*s23*'])
; KK04
    oh_kk04 = mz_assign_r23branch(abund,/kk04,debug=debug,$
      ewalpha=ewalpha,justew=justew,justflux=justflux,$
      r23branch=r23branch)
    oh_kk04 = im_struct_trimtags(oh_kk04,select=tag_names(oh_kk04),$
      newtags=repstr(tag_names(oh_kk04),'R23BRANCH','R23BRANCH_KK04'))
; M91
    oh_m91 = mz_assign_r23branch(abund,/m91,debug=debug,$
      ewalpha=ewalpha,justew=justew,justflux=justflux,$
      r23branch=r23branch)
    oh_m91 = struct_trimtags(struct_trimtags(oh_m91,$
      select=['r23branch*','*m91*']),except='*LOGU*')
    oh_m91 = im_struct_trimtags(oh_m91,select=tag_names(oh_m91),$
      newtags=repstr(tag_names(oh_m91),'R23BRANCH','R23BRANCH_M91'))
; T04: the T04 abundances are a special case (only upper-branch); but
; reject abundances with R23>10
    m91tags = tag_names(oh_m91)
    oh_t04 = im_struct_trimtags(oh_m91,select=$
      m91tags,newtags=repstr(m91tags,'M91','T04'))
    oh_t04 = im_struct_assign(abund,oh_t04,/nozero)
    oh_t04.r23branch_t04 = 'U'
    oh_t04.r23branch_t04_ew = 'U'
    rej = where(abund.zstrong_r23 gt 10.0,nrej)
    if (nrej ne 0L) then begin
       oh_t04[rej].r23branch_t04 = 'Rejected'
       oh_t04[rej].zstrong_12oh_t04 = -999.0
       oh_t04[rej].zstrong_12oh_t04_err = -999.0
    endif
    ewrej = where(abund.zstrong_ew_r23 gt 10.0,newrej)
    if (newrej ne 0L) then begin
       oh_t04[ewrej].r23branch_t04_ew = 'Rejected'
       oh_t04[ewrej].zstrong_ew_12oh_t04 = -999.0
       oh_t04[ewrej].zstrong_ew_12oh_t04_err = -999.0
    endif
; rename the tags and merge everything
    oh = struct_addtags(struct_addtags(oh_kk04,oh_m91),oh_t04)
    oldtags = tag_names(oh)
    oh = im_struct_trimtags(oh,select=oldtags,$
      newtags=repstr(repstr(oldtags,'ZSTRONG_',''),'12OH','LOG12OH'))
    oh = struct_addtags(temporary(oh),struct_trimtags(abund,select=['ohlimit']))
return, oh
end

pro build_mz_log12oh, sdss=sdss, agn=agn, debug=debug, $
  nmonte=nmonte, clobber=clobber
; jm09mar18nyu - build the AGES and SDSS abundances samples (i.e., AGN
;   have been removed and the R23 abundances are well-behaved) 
; jm10jul29ucsd - near-total rewrite

    mzpath = ages_path(/projects)+'mz/'
    nmonte = 100
;   nmonte = 250 ; 500
    
    if keyword_set(sdss) then begin ; SDSS
       prefix = 'sdss'
       snrcut1 = 2.0
    endif else begin            ; AGES
       prefix = 'ages'
       snrcut1 = 2.0
    endelse

; optionally read the AGN sample (for testing with the SDSS sample)
    if keyword_set(agn) then begin
       mzagn_ispec = 1
       mzhii_ispec = 0
       mzagn_ancillary = 1
       mzhii_ancillary = 0
       rootfile = prefix+'_mz_agn_log12oh'
    endif else begin
       mzagn_ispec = 0
       mzhii_ispec = 1
       mzagn_ancillary = 0
       mzhii_ancillary = 1
       rootfile = prefix+'_mz_hii_log12oh'
    endelse

    dust = read_mz_sample(mzhii_ispec=mzhii_ispec,$
      mzagn_ispec=mzagn_ispec,sdss=sdss)
;   dust = dust[1.5E4:1.55E4]
;   dust = dust[0:100]
    ngal = n_elements(dust)

; determine R23 branches based on [NII]/Ha, assuming a default upper
; branch 
    splog, 'Assigning R23 branches'
    r23branch = replicate('U',ngal) ; default
    lodetect = (dust.bpt_nii_ha gt -900) and (dust.bpt_nii_ha lt -1.1)
    lolimit = (dust.bpt_nii_ha_limit gt -900) and (dust.bpt_nii_ha_limit lt -1.1)
    lo = where(lodetect or lolimit,nlo)
;   if tag_exist(dust,'bpt_nii_ha_limit') then begin
;      lolimit = (dust.bpt_nii_ha_limit gt -900) and (dust.bpt_nii_ha_limit lt -1.1)
;      lo = where(lodetect or lolimit,nlo)
;   endif else lo = where(lodetect,nlo)
    if (nlo ne 0L) then r23branch[lo] = 'L'

; compute abundances from observed fluxes and EWs assuming alpha=1;
; deal with upper limits on [OIII]: compute the error on the
; abundances, etc., for these objects by fixing [OIII] at the upper
; limit and at zero and take the difference
    adust = mz_abundance(dust,nmonte=nmonte,snrcut=0.0,ewalpha=1.0)

    splog, 'Computing upper limits on the metallicities'
    lim = where(dust.oiii_5007[1] eq -1.0,nlim)
    if (nlim ne 0L) then begin
       dust_lim = dust[lim]
; at the upper limit
       dust_lim.oiii_5007 = transpose([[dust_lim.oiii_5007_limit],[dust_lim.oiii_5007_limit*0.01]])
       dust_lim.oiii_5007_ew = transpose([[dust_lim.oiii_5007_ew_limit],[dust_lim.oiii_5007_limit*0.01]])
       adust_lim1 = mz_abundance(dust_lim,nmonte=0,snrcut=0.0,ewalpha=1.0)
; assume [OIII]=0 (actually, just a very small number)
       dust_lim.oiii_5007 = [1D-21,1D-22]
       dust_lim.oiii_5007_ew = [1D-5,1D-7]
       adust_lim2 = mz_abundance(dust_lim,nmonte=0,snrcut=0.0,ewalpha=1.0)

       tags = 'zstrong_'+['r23','o32','p','12oh_'+[$
         'm91_'+['avg','upper','lower'],$
         'kk04_'+['avg','upper','lower']],$
         'logu_kk04_'+['avg','upper','lower'],'12oh_t04']
       tags = [tags,repstr(tags,'zstrong_','zstrong_ew_')]
       for jj = 0, n_elements(tags)-1 do begin
          indx1 = tag_indx(adust_lim1,tags[jj])
          indx2 = tag_indx(adust_lim1,tags[jj]+'_err')
          adust_lim1.(indx2) = abs(adust_lim1.(indx1)-adust_lim2.(indx1))
;         niceprint, adust_lim1.(indx1), adust_lim2.(indx1)
       endfor
       adust[lim] = adust_lim1
    endif

; this has to go here    
    adust = struct_addtags(temporary(adust),replicate({ohlimit: 0},ngal))
    if (nlim ne 0L) then adust[lim].ohlimit = 1
    
    ohdust = get_log12oh(adust,r23branch=r23branch,debug=debug,$
      justew=justew,justflux=justflux)

; compute abundances from reddening-corrected fluxes; ignore upper
; limits here 
    nodust = iunred_linedust(dust,snrcut=snrcut1)
    anodust = mz_abundance(nodust,nmonte=nmonte,snrcut=0.0,/justflux)
    anodust = struct_addtags(temporary(anodust),replicate({ohlimit: 0},ngal))
    
    ohnodust = get_log12oh(anodust,r23branch=r23branch,debug=debug,/justflux)
    ohnodust = struct_addtags(struct_trimtags(nodust,select='*hahb*'),$
      struct_trimtags(ohnodust,except='*ew*'))

; for the SDSS also add the Tremonti+ metallicity to both output
; structures (for convenience)
    if keyword_set(sdss) then begin
       anc = read_mz_sample(mzhii_ancillary=mzhii_ancillary,$
         mzagn_ancillary=mzagn_ancillary,sdss=sdss)
       more = replicate({log12oh_tremonti: -999.0, $
         log12oh_tremonti_err: -999.0},ngal)
       good = where((anc.oh_median gt 0.0) and $
         (anc.oh_p84 gt 0.0) and (anc.oh_p16 gt 0.0),ngood)
       more[good].log12oh_tremonti = anc[good].oh_median
       more[good].log12oh_tremonti_err = (anc[good].oh_p84-anc[good].oh_p16)/2.0
       ohdust = struct_addtags(temporary(ohdust),more)
       ohnodust = struct_addtags(temporary(ohnodust),more)
    endif

; write out    
    im_mwrfits, ohdust, mzpath+rootfile+'.fits', clobber=clobber
    im_mwrfits, ohnodust, mzpath+rootfile+'_nodust.fits', clobber=clobber

return
end
