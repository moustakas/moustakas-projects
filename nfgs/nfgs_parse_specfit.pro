pro nfgs_parse_specfit, line, _extra=extra
; jm04jan22uofa
; jm04mar21uofa - updated
    
; noticeable second-order contamination redward of ~6700 A
    
;   redbad = [$
;     'A00113+3037','A00289+0556','A00389-0159','A00442+3224','A00570+1504', $
;     'A01344+2838','A02056+1444','A02464+1807','A10321+4649','A10337+1358', $
;     'A10465+0711','A11142+1804','A11310+3254','A12331+7230', $
;     'A15523+1645','A22306+0750','A23176+1541','A23264+1703','A23514+2813', $
;     'IC1124','IC1504','IC197','IC2520', 'NGC1029','NGC1298','NGC1552', $
;     'NGC193','NGC2692','NGC315','NGC3279','NGC3633','NGC382','NGC516', $
;     'NGC5425','NGC5541','NGC5888','NGC7194','NGC7328','NGC7360','NGC7436',$
;     'NGC7460','NGC7537','NGC825','NGC984']

    datapath = nfgs_path(/specfit)
    outpath = getenv('CATALOGS_DIR')+'/nfgs/'
    analysis_path = nfgs_path(/analysis)

    root = 'nfgs_int'
    outfile = root+'_speclinefit.fits'

; read the ancillary nfgs data, the synthesized magnitudes, and the
; stellar masses

    infofile = 'nfgs_info.fits.gz'
    if (file_test(analysis_path+infofile,/regular) eq 0L) then begin
       splog, 'Unable to find '+analysis_path+infofile+'.'
       return
    endif

    splog, 'Reading '+analysis_path+infofile+'.'
    nfgs = mrdfits(analysis_path+infofile,1,/silent)
    nfgs = struct_trimtags(nfgs,except=except)
 
    synthfile = root+'_synthmags.fits.gz'
    if file_test(datapath+synthfile,/regular) then begin
       splog, 'Reading and appending '+datapath+synthfile+'.'
       mags = mrdfits(datapath+synthfile,1,/silent)
       match, strtrim(nfgs.galaxy,2), strtrim(mags.galaxy,2), nfgsindx, magsindx
       nfgs = struct_addtags(nfgs[nfgsindx],struct_trimtags(mags[magsindx],except='GALAXY'))
    endif

    massfile = 'nfgs_stellar_mass_photo.fits.gz'
    if file_test(analysis_path+massfile,/regular)then begin
       splog, 'Reading and appending '+analysis_path+massfile+'.'
       mass = mrdfits(analysis_path+massfile,1,/silent)
       match, strtrim(nfgs.galaxy,2), strtrim(mass.galaxy,2), nfgsindx, massindx
       nfgs = struct_addtags(nfgs[nfgsindx],struct_trimtags(mass[massindx],except='GALAXY'))
    endif

    kcorrmassfile = 'nfgs_stellar_mass_kcorrect.fits.gz'
    if file_test(analysis_path+kcorrmassfile,/regular)then begin
       splog, 'Reading and appending '+analysis_path+kcorrmassfile+'.'
       kcorrmass = mrdfits(analysis_path+kcorrmassfile,1,/silent)
       match, strtrim(nfgs.galaxy,2), strtrim(kcorrmass.galaxy,2), nfgsindx, kcorrmassindx
       nfgs = struct_addtags(nfgs[nfgsindx],struct_trimtags(kcorrmass[kcorrmassindx],except='GALAXY'))
    endif

; ---------------------------------------------------------------------------    
; UBVR photometry from Jansen et al. (2000) has been tabulated in
; NFGS_INFO.  here, correct the photometry for emission lines but do
; not let the correction be less than zero; use synthesized photometry
; where no other photometry is available; we could also correct for
; internal extinction here, or to face-on inclination; (see also
; WRITE_00JANSEN)
; ---------------------------------------------------------------------------    

    if file_test(datapath+synthfile,/regular) then begin
       
; U

       cor = (nfgs.synth_u-nfgs.synth_u_obs) > 0.0 ; <-- NOTE!
       good = where(nfgs.u_obs gt -900.0,ngood,comp=synth,ncomp=nsynth)
       if (ngood ne 0L) then begin
          nfgs[good].u         = nfgs[good].u_obs + cor[good]
          nfgs[good].u_err     = nfgs[good].u_obs_err
          nfgs[good].M_u       = nfgs[good].M_u_obs + cor[good]
          nfgs[good].M_u_err   = nfgs[good].M_u_obs_err
          nfgs[good].u_lum     = nfgs[good].u_lum_obs - 0.4*cor[good]
          nfgs[good].u_lum_err = nfgs[good].u_lum_obs_err
       endif
       if (nsynth ne 0L) then begin
          nfgs[synth].u_obs         = nfgs[synth].synth_u_obs
          nfgs[synth].u_obs_err     = nfgs[synth].synth_u_obs_err
          nfgs[synth].M_u_obs       = nfgs[synth].synth_M_u + cor[synth]
          nfgs[synth].M_u_obs_err   = nfgs[synth].synth_M_u_err
          nfgs[synth].u_lum_obs     = nfgs[synth].synth_u_lum - 0.4*cor[synth]
          nfgs[synth].u_lum_obs_err = nfgs[synth].synth_u_lum_err

          nfgs[synth].u         = nfgs[synth].synth_u
          nfgs[synth].u_err     = nfgs[synth].synth_u_err
          nfgs[synth].M_u       = nfgs[synth].synth_M_u
          nfgs[synth].M_u_err   = nfgs[synth].synth_M_u_err
          nfgs[synth].u_lum     = nfgs[synth].synth_u_lum
          nfgs[synth].u_lum_err = nfgs[synth].synth_u_lum_err
       endif

; B
       
       cor = (nfgs.synth_b-nfgs.synth_b_obs) > 0.0 ; <-- NOTE!
       good = where(nfgs.b_obs gt -900.0,ngood,comp=synth,ncomp=nsynth)
       if (ngood ne 0L) then begin
          nfgs[good].b         = nfgs[good].b_obs + cor[good]
          nfgs[good].b_err     = nfgs[good].b_obs_err
          nfgs[good].M_b       = nfgs[good].M_b_obs + cor[good]
          nfgs[good].M_b_err   = nfgs[good].M_b_obs_err
          nfgs[good].b_lum     = nfgs[good].b_lum_obs - 0.4*cor[good]
          nfgs[good].b_lum_err = nfgs[good].b_lum_obs_err
       endif
       if (nsynth ne 0L) then begin
          nfgs[synth].b_obs         = nfgs[synth].synth_b_obs
          nfgs[synth].b_obs_err     = nfgs[synth].synth_b_obs_err
          nfgs[synth].M_b_obs       = nfgs[synth].synth_M_b + cor[synth]
          nfgs[synth].M_b_obs_err   = nfgs[synth].synth_M_b_err
          nfgs[synth].b_lum_obs     = nfgs[synth].synth_b_lum - 0.4*cor[synth]
          nfgs[synth].b_lum_obs_err = nfgs[synth].synth_b_lum_err

          nfgs[synth].b         = nfgs[synth].synth_b
          nfgs[synth].b_err     = nfgs[synth].synth_b_err
          nfgs[synth].M_b       = nfgs[synth].synth_M_b
          nfgs[synth].M_b_err   = nfgs[synth].synth_M_b_err
          nfgs[synth].b_lum     = nfgs[synth].synth_b_lum
          nfgs[synth].b_lum_err = nfgs[synth].synth_b_lum_err
       endif

; V    
       
       cor = (nfgs.synth_v-nfgs.synth_v_obs) > 0.0 ; <-- NOTE!
       good = where(nfgs.v_obs gt -900.0,ngood,comp=synth,ncomp=nsynth)
       if (ngood ne 0L) then begin
          nfgs[good].v         = nfgs[good].v_obs + cor[good]
          nfgs[good].v_err     = nfgs[good].v_obs_err
          nfgs[good].M_v       = nfgs[good].M_v_obs + cor[good]
          nfgs[good].M_v_err   = nfgs[good].M_v_obs_err
          nfgs[good].v_lum     = nfgs[good].v_lum_obs - 0.4*cor[good]
          nfgs[good].v_lum_err = nfgs[good].v_lum_obs_err
       endif
       if (nsynth ne 0L) then begin
          nfgs[synth].v_obs         = nfgs[synth].synth_v_obs
          nfgs[synth].v_obs_err     = nfgs[synth].synth_v_obs_err
          nfgs[synth].M_v_obs       = nfgs[synth].synth_M_v + cor[synth]
          nfgs[synth].M_v_obs_err   = nfgs[synth].synth_M_v_err
          nfgs[synth].v_lum_obs     = nfgs[synth].synth_v_lum - 0.4*cor[synth]
          nfgs[synth].v_lum_obs_err = nfgs[synth].synth_v_lum_err

          nfgs[synth].v         = nfgs[synth].synth_v
          nfgs[synth].v_err     = nfgs[synth].synth_v_err
          nfgs[synth].M_v       = nfgs[synth].synth_M_v
          nfgs[synth].M_v_err   = nfgs[synth].synth_M_v_err
          nfgs[synth].v_lum     = nfgs[synth].synth_v_lum
          nfgs[synth].v_lum_err = nfgs[synth].synth_v_lum_err
       endif

; R    
       
       cor = (nfgs.synth_r-nfgs.synth_r_obs) > 0.0 ; <-- NOTE!
       good = where(nfgs.r_obs gt -900.0,ngood,comp=synth,ncomp=nsynth)
       if (ngood ne 0L) then begin
          nfgs[good].r         = nfgs[good].r_obs + cor[good]
          nfgs[good].r_err     = nfgs[good].r_obs_err
          nfgs[good].M_r       = nfgs[good].M_r_obs + cor[good]
          nfgs[good].M_r_err   = nfgs[good].M_r_obs_err
          nfgs[good].r_lum     = nfgs[good].r_lum_obs - 0.4*cor[good]
          nfgs[good].r_lum_err = nfgs[good].r_lum_obs_err
       endif
       if (nsynth ne 0L) then begin
          nfgs[synth].r_obs         = nfgs[synth].synth_r_obs
          nfgs[synth].r_obs_err     = nfgs[synth].synth_r_obs_err
          nfgs[synth].M_r_obs       = nfgs[synth].synth_M_r + cor[synth]
          nfgs[synth].M_r_obs_err   = nfgs[synth].synth_M_r_err
          nfgs[synth].r_lum_obs     = nfgs[synth].synth_r_lum - 0.4*cor[synth]
          nfgs[synth].r_lum_obs_err = nfgs[synth].synth_r_lum_err

          nfgs[synth].r         = nfgs[synth].synth_r
          nfgs[synth].r_err     = nfgs[synth].synth_r_err
          nfgs[synth].M_r       = nfgs[synth].synth_M_r
          nfgs[synth].M_r_err   = nfgs[synth].synth_M_r_err
          nfgs[synth].r_lum     = nfgs[synth].synth_r_lum
          nfgs[synth].r_lum_err = nfgs[synth].synth_r_lum_err
       endif

    endif

; write out    
    
    select_lines = ['OII_3727','H_BETA','OIII_5007','H_ALPHA']
    disttag = 'DISTANCE'
    disterrtag = 'DISTANCE_ERR'

    snrcut_linedust = 3.0
    snrcut_abundance = 3.0

    syserr = 0.0 ; DO NOT USE!!
;   syserr = 6.0 ; systematic error [%]

; parse everything 
    
    line = parse_ispeclinefit(datapath=datapath,outpath=outpath,prepend=nfgs,$
      root=root,trimtags='GALAXY',select_lines=select_lines,syserr=syserr,$
      disttag=disttag,disterrtag=disterrtag,photerrtag=photerrtag,$
      outfile=outfile,electrondensity=0,snrcut_linedust=snrcut_linedust,$
      snrcut_abundance=snrcut_abundance,/match,/kauffmann,/write,$
      /odonnell,/nopropagate,_extra=extra)

return
end
