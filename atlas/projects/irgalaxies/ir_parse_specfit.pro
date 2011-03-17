pro ir_parse_specfit, line, _extra=extra
; jm05aug15uofa - parse 

    datapath = atlas_path(/projects)+'irgalaxies/specfit/'
    atlaspath = atlas_path(/specfit)
    analysis_path = atlas_path(/analysis)

    root = 'ir_04'
    outfile = root+'_speclinefit.fits'

; read the ancillary atlas data, the synthesized magnitudes, and the
; stellar masses

    infofile = 'atlas1d_info.fits.gz'
    if (file_test(analysis_path+infofile,/regular) eq 0L) then begin
       splog, 'Unable to find '+analysis_path+infofile+'.'
       return
    endif

    splog, 'Reading '+analysis_path+infofile+'.'
    atlas = mrdfits(analysis_path+infofile,1,/silent)
    atlas = struct_trimtags(atlas,except=except)

    synthfile = 'integrated_atlas_synthmags.fits.gz' ; NOTE CHANGE IN FILE/PATH!
    if file_test(atlaspath+synthfile,/regular) then begin 
       splog, 'Reading and appending '+atlaspath+synthfile+'.'
       mags = mrdfits(atlaspath+synthfile,1,/silent)
       match, strtrim(atlas.galaxy,2), strtrim(mags.galaxy,2), atlasindx, magsindx
       atlas = struct_addtags(atlas[atlasindx],struct_trimtags(mags[magsindx],except='GALAXY'))
    endif

    massfile = 'atlas_stellar_mass_photo.fits.gz'
    if file_test(analysis_path+massfile,/regular)then begin
       splog, 'Reading and appending '+analysis_path+massfile+'.'
       mass = mrdfits(analysis_path+massfile,1,/silent)
       match, strtrim(atlas.galaxy,2), strtrim(mass.galaxy,2), atlasindx, massindx
       atlas = struct_addtags(atlas[atlasindx],struct_trimtags(mass[massindx],except='GALAXY'))
    endif

; ---------------------------------------------------------------------------    
; finalize the UBV photometry; correct the RC3 photometry for emission
; lines but do not let the correction be less than zero; use
; synthesized photometry where no RC3 photometry is available; we
; could also correct for internal extinction here, or to face-on
; inclination
; ---------------------------------------------------------------------------    

    if file_test(atlaspath+synthfile,/regular) then begin
       
; U

       good = where(atlas.rc3_u gt -900.0,ngood,comp=synth,ncomp=nsynth)
       if (ngood ne 0L) then begin
          cor = (atlas.synth_u-atlas.synth_u_obs) > 0.0 ; <-- NOTE!
          atlas[good].u         = atlas[good].rc3_u + cor[good]
          atlas[good].u_err     = atlas[good].rc3_u_err
          atlas[good].M_u       = atlas[good].rc3_M_u + cor[good]
          atlas[good].M_u_err   = atlas[good].rc3_M_u_err
          atlas[good].u_lum     = atlas[good].rc3_u_lum - 0.4*cor[good]
          atlas[good].u_lum_err = atlas[good].rc3_u_lum_err
       endif
       if (nsynth ne 0L) then begin
          atlas[synth].u         = atlas[synth].synth_u
          atlas[synth].u_err     = atlas[synth].synth_u_err
          atlas[synth].M_u       = atlas[synth].synth_M_u
          atlas[synth].M_u_err   = atlas[synth].synth_M_u_err
          atlas[synth].u_lum     = atlas[synth].synth_u_lum
          atlas[synth].u_lum_err = atlas[synth].synth_u_lum_err
       endif

; B
       
       good = where(atlas.rc3_b gt -900.0,ngood,comp=synth,ncomp=nsynth)
       if (ngood ne 0L) then begin
          cor = (atlas.synth_b-atlas.synth_b_obs) > 0.0 ; <-- NOTE!
          atlas[good].b         = atlas[good].rc3_b + cor[good]
          atlas[good].b_err     = atlas[good].rc3_b_err
          atlas[good].M_b       = atlas[good].rc3_M_b + cor[good]
          atlas[good].M_b_err   = atlas[good].rc3_M_b_err
          atlas[good].b_lum     = atlas[good].rc3_b_lum - 0.4*cor[good]
          atlas[good].b_lum_err = atlas[good].rc3_b_lum_err
       endif
       if (nsynth ne 0L) then begin
          atlas[synth].b         = atlas[synth].synth_b
          atlas[synth].b_err     = atlas[synth].synth_b_err
          atlas[synth].M_b       = atlas[synth].synth_M_b
          atlas[synth].M_b_err   = atlas[synth].synth_M_b_err
          atlas[synth].b_lum     = atlas[synth].synth_b_lum
          atlas[synth].b_lum_err = atlas[synth].synth_b_lum_err
       endif

; V    
       
       good = where(atlas.rc3_v gt -900.0,ngood,comp=synth,ncomp=nsynth)
       if (ngood ne 0L) then begin
          cor = (atlas.synth_v-atlas.synth_v_obs) > 0.0 ; <-- NOTE!
          atlas[good].v         = atlas[good].rc3_v + cor[good]
          atlas[good].v_err     = atlas[good].rc3_v_err
          atlas[good].M_v       = atlas[good].rc3_M_v + cor[good]
          atlas[good].M_v_err   = atlas[good].rc3_M_v_err
          atlas[good].v_lum     = atlas[good].rc3_v_lum - 0.4*cor[good]
          atlas[good].v_lum_err = atlas[good].rc3_v_lum_err
       endif
       if (nsynth ne 0L) then begin
          atlas[synth].v         = atlas[synth].synth_v
          atlas[synth].v_err     = atlas[synth].synth_v_err
          atlas[synth].M_v       = atlas[synth].synth_M_v
          atlas[synth].M_v_err   = atlas[synth].synth_M_v_err
          atlas[synth].v_lum     = atlas[synth].synth_v_lum
          atlas[synth].v_lum_err = atlas[synth].synth_v_lum_err
       endif

    endif
       
; write out    
    
    select_lines = ['OII_3727','H_BETA','OIII_5007','H_ALPHA']
    disttag = 'DISTANCE'
    disterrtag = 'DISTANCE_ERR'
;   photerrtag = 'DRIFTABSERROR'

    snrcut_linedust = 3.0
    snrcut_abundance = 3.0

; parse everything 
    
    line = parse_ispeclinefit(datapath=datapath,prepend=atlas,root=root,$
      trimtags='GALAXY',select_lines=select_lines,$
      disttag=disttag,disterrtag=disterrtag,photerrtag=photerrtag,$
      outfile=outfile,/electrondensity,snrcut_linedust=snrcut_linedust,$
      snrcut_abundance=snrcut_abundance,/match,/kauffmann,/write,$
      /odonnell,/nopropagate,_extra=extra)

return
end
