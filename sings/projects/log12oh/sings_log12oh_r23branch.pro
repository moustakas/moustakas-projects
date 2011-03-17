pro sings_log12oh_r23branch, branchinfo, clobber=clobber
; jm06feb10uofa - attempt to pick the correct R23 branch based on the
;   nuclear, circumnuclear, and radial strip spectra
; jm06oct19nyu - the procedure is now iterative; the final branches
;   are in the sextractor file, but the decision was based on
;   considering all three spectra and the B-band luminosity  
; jm08feb07nyu - updated
; jm10mar10ucsd - final version; note that this is an ancillary
;   support routine to help build the SINGS_R23_BRANCH.SEX file

    version = sings_log12oh_version()
    outpath = sings_path(/projects)+'log12oh/'

    sings = sings_read_info()
    ngal = n_elements(sings)
    
    r23info = rsex(outpath+'sings_r23_branch.sex')
    nuclear = read_sings_gandalf(/nuclear)
    drift20 = read_sings_gandalf(/drift20)
    drift56 = read_sings_gandalf(/drift56)

    nuclear = read_sings_log12oh_samples(/nodust_nuclear)
    drift20 = read_sings_log12oh_samples(/nodust_drift20)
    drift56 = read_sings_log12oh_samples(/nodust_drift56)

    branchinfo = {galaxy: '', m_b: -999.0, r23_branch: '?', $
      nuc_branch: '?', d20_branch: '?', d56_branch: '?', $
      nuc_reason: '-', d20_reason: '-', d56_reason: '-', $
      nuc_r23: -999.0, nuc_n2: -999.0, nuc_n2o2: -999.0, $
      d20_r23: -999.0, d20_n2: -999.0, d20_n2o2: -999.0, $
      d56_r23: -999.0, d56_n2: -999.0, d56_n2o2: -999.0}
    branchinfo = replicate(branchinfo,ngal)

    branchinfo.galaxy = sings.galaxy
    branchinfo.m_b = sings.rc3_ubv_absmag[1]
    branchinfo.r23_branch = r23info.r23_branch

    for k = 0L, ngal-1L do begin

       match = where(sings[k].ned_galaxy eq nuclear.ned_galaxy,nmatch)
       if (nmatch ne 0L) then begin
          branchinfo[k].nuc_r23           = nuclear[match].zstrong_r23
          branchinfo[k].nuc_n2            = nuclear[match].zstrong_niiha
          branchinfo[k].nuc_n2o2          = nuclear[match].zstrong_niioii
;         branchinfo[k].nuc_oh_o3n2       = nuclear[match].zstrong_12oh_oiiinii_niiha
;         branchinfo[k].nuc_oh_kk04_upper = nuclear_nodust[match].zstrong_12oh_kk04_r23_upper
;         branchinfo[k].nuc_oh_kk04_lower = nuclear_nodust[match].zstrong_12oh_kk04_r23_lower
       endif

       match = where(sings[k].ned_galaxy eq drift20.ned_galaxy,nmatch)
       if (nmatch ne 0L) then begin
          branchinfo[k].d20_r23           = drift20[match].zstrong_r23
          branchinfo[k].d20_n2            = drift20[match].zstrong_niiha
          branchinfo[k].d20_n2o2          = drift20[match].zstrong_niioii
;         branchinfo[k].d20_oh_o3n2       = drift20[match].zstrong_12oh_oiiinii_niiha
;         branchinfo[k].d20_oh_kk04_upper = drift20_nodust[match].zstrong_12oh_kk04_r23_upper
;         branchinfo[k].d20_oh_kk04_lower = drift20_nodust[match].zstrong_12oh_kk04_r23_lower
       endif

       match = where(sings[k].ned_galaxy eq drift56.ned_galaxy,nmatch)
       if (nmatch ne 0L) then begin
          branchinfo[k].d56_r23           = drift56[match].zstrong_r23
          branchinfo[k].d56_n2            = drift56[match].zstrong_niiha
          branchinfo[k].d56_n2o2          = drift56[match].zstrong_niioii
;         branchinfo[k].d56_oh_o3n2       = drift56[match].zstrong_12oh_oiiinii_niiha
;         branchinfo[k].d56_oh_kk04_upper = drift56_nodust[match].zstrong_12oh_kk04_r23_upper
;         branchinfo[k].d56_oh_kk04_lower = drift56_nodust[match].zstrong_12oh_kk04_r23_lower
       endif

    endfor

; nuclear    

    both = where((branchinfo.nuc_n2 gt -900.0) and (branchinfo.nuc_n2o2 gt -900.0),nboth)
    if (nboth ne 0L) then begin
       branchinfo[both].nuc_reason = 'N2,N2O2'
       lo = where((branchinfo[both].nuc_n2 lt -1.0) and (branchinfo[both].nuc_n2o2 lt -1.05),nlo)
       if (nlo ne 0L) then branchinfo[both[lo]].nuc_branch = 'L'
       up = where((branchinfo[both].nuc_n2 gt -1.0) and (branchinfo[both].nuc_n2o2 gt -0.8),nup)
       if (nup ne 0L) then branchinfo[both[up]].nuc_branch = 'U'
       ambig = where((branchinfo[both].nuc_n2o2 gt -1.05) and (branchinfo[both].nuc_n2o2 lt -0.8),nambig)
       if (nambig ne 0L) then branchinfo[both[ambig]].nuc_branch = 'A'
    endif

    single = where((branchinfo.nuc_n2 gt -900.0) and (branchinfo.nuc_n2o2 lt -900.0),nsingle)
    if (nsingle ne 0L) then begin
       branchinfo[single].nuc_reason = 'N2'
       lo = where((branchinfo[single].nuc_n2 lt -1.0),nlo)
       if (nlo ne 0L) then branchinfo[single[lo]].nuc_branch = 'L'
       up = where((branchinfo[single].nuc_n2 gt -1.0),nup)
       if (nup ne 0L) then branchinfo[single[up]].nuc_branch = 'U'
    endif

; drift20
    
    both = where((branchinfo.d20_n2 gt -900.0) and (branchinfo.d20_n2o2 gt -900.0),nboth)
    if (nboth ne 0L) then begin
       branchinfo[both].d20_reason = 'N2,N2O2'
       lo = where((branchinfo[both].d20_n2 lt -1.0) and (branchinfo[both].d20_n2o2 lt -1.05),nlo)
       if (nlo ne 0L) then branchinfo[both[lo]].d20_branch = 'L'
       up = where((branchinfo[both].d20_n2 gt -1.0) and (branchinfo[both].d20_n2o2 gt -0.8),nup)
       if (nup ne 0L) then branchinfo[both[up]].d20_branch = 'U'
       ambig = where((branchinfo[both].d20_n2o2 gt -1.05) and (branchinfo[both].d20_n2o2 lt -0.8),nambig)
       if (nambig ne 0L) then branchinfo[both[ambig]].d20_branch = 'A'
    endif

    single = where((branchinfo.d20_n2 gt -900.0) and (branchinfo.d20_n2o2 lt -900.0),nsingle)
    if (nsingle ne 0L) then begin
       branchinfo[single].d20_reason = 'N2'
       lo = where((branchinfo[single].d20_n2 lt -1.0),nlo)
       if (nlo ne 0L) then branchinfo[single[lo]].d20_branch = 'L'
       up = where((branchinfo[single].d20_n2 gt -1.0),nup)
       if (nup ne 0L) then branchinfo[single[up]].d20_branch = 'U'
    endif

; drift56    
    
    both = where((branchinfo.d56_n2 gt -900.0) and (branchinfo.d56_n2o2 gt -900.0),nboth)
    if (nboth ne 0L) then begin
       branchinfo[both].d56_reason = 'N2,N2O2'
       lo = where((branchinfo[both].d56_n2 lt -1.0) and (branchinfo[both].d56_n2o2 lt -1.05),nlo)
       if (nlo ne 0L) then branchinfo[both[lo]].d56_branch = 'L'
       up = where((branchinfo[both].d56_n2 gt -1.0) and (branchinfo[both].d56_n2o2 gt -0.8),nup)
       if (nup ne 0L) then branchinfo[both[up]].d56_branch = 'U'
       ambig = where((branchinfo[both].d56_n2o2 gt -1.05) and (branchinfo[both].d56_n2o2 lt -0.8),nambig)
       if (nambig ne 0L) then branchinfo[both[ambig]].d56_branch = 'A'
    endif

    single = where((branchinfo.d56_n2 gt -900.0) and (branchinfo.d56_n2o2 lt -900.0),nsingle)
    if (nsingle ne 0L) then begin
       branchinfo[single].d56_reason = 'N2'
       lo = where((branchinfo[single].d56_n2 lt -1.0),nlo)
       if (nlo ne 0L) then branchinfo[single[lo]].d56_branch = 'L'
       up = where((branchinfo[single].d56_n2 gt -1.0),nup)
       if (nup ne 0L) then branchinfo[single[up]].d56_branch = 'U'
    endif

; write out
    splog, 'Writing '+outpath+'sings_r23_branch_'+version+'.dat'
    splog, 'Writing '+outpath+'sings_r23_branch_'+version+'.fits'
    struct_print, branchinfo, file=outpath+'sings_r23_branch_'+version+'.dat'
    im_mwrfits, branchinfo, outpath+'sings_r23_branch_'+version+'.fits', clobber=clobber

; old code below here:

;   flag = where(((nuc.branch eq '?') or (nuc.branch eq 'A') or (strmatch(nuc.galaxy,'*tololo89*',/fold)) or $
;     (strmatch(nuc.galaxy,'*ngc1705*',/fold))) and (nuc.r23 gt -900.0),nflag)
;   if (nflag ne 0L) then begin
;      for i = 0L, nflag-1L do begin
;         match = where(strtrim(nuc[flag[i]].galaxy,2) eq strtrim(r23info.galaxy,2))
;         splog, 'Manually setting R23 branch for '+strtrim(nuc[flag[i]].galaxy,2)+': '+$
;           nuc[flag[i]].branch+'-->'+r23info[match].r23_branch+'.'
;         nuc[flag[i]].branch = r23info[match].r23_branch
;      endfor
;      print
;   endif
;   
;   flag = where(((d20.branch eq '?') or (d20.branch eq 'A') or (strmatch(d20.galaxy,'*tololo89*',/fold))) and (d20.r23 gt -900.0),nflag)
;   if (nflag ne 0L) then begin
;      for i = 0L, nflag-1L do begin
;         match = where(strtrim(d20[flag[i]].galaxy,2) eq strtrim(r23info.galaxy,2))
;         splog, 'Manually setting R23 branch for '+strtrim(d20[flag[i]].galaxy,2)+': '+$
;           d20[flag[i]].branch+'-->'+r23info[match].r23_branch+'.'
;         d20[flag[i]].branch = r23info[match].r23_branch
;      endfor
;      print
;   endif
;
;   flag = where(((d56.branch eq '?') or (d56.branch eq 'A') or (strmatch(d56.galaxy,'*tololo89*',/fold))) and (d56.r23 gt -900.0),nflag)
;   if (nflag ne 0L) then begin
;      for i = 0L, nflag-1L do begin
;         match = where(strtrim(d56[flag[i]].galaxy,2) eq strtrim(r23info.galaxy,2))
;         splog, 'Manually setting R23 branch for '+strtrim(d56[flag[i]].galaxy,2)+': '+$
;           d56[flag[i]].branch+'-->'+r23info[match].r23_branch+'.'
;         d56[flag[i]].branch = r23info[match].r23_branch
;      endfor
;      print
;   endif

return
end
    
