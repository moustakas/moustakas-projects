;+
; NAME:
;       ATLAS_PARSE_SPECFIT
;
; PURPOSE:
;       Parse the output from fitting the atlas.
;
; CALLING SEQUENCE:
;       atlas_parse_specfit, datapath=, /nuclear, _extra=extra
;
; INPUTS:
;       None.
;
; OPTIONAL INPUTS:
;       extra - keywords for PARSE_ISPECLINEFIT()
;
; KEYWORD PARAMETERS:
;       nuclear - parse the nuclear spectral fitting
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       line - output from PARSE_ISPECLINEFIT()
;
; PROCEDURES USED:
;       STRUCT_TRIMTAGS(), PARSE_ISPECLINEFIT() 
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Jun 13, U of A - written
;       jm03dec04uofa - generalized
;       jm04mar11uofa - streamlined and now calls PARSE_ISPECLINEFIT() 
;       jm08apr14nyu - remove the corrections for emission lines (the
;                      requisite tags were removed from
;                      INIT_ANCILLARY_DATA); also relegated
;                      concatenating the stellar masses and
;                      synthesized magnitudes
;-

pro atlas_parse_specfit, atlas_final, nuclear=nuclear

    version = atlas_version(/specfit)
    
    datapath = atlas_path(/specfit)
    outpath = atlas_path(/spectral_atlas)
    analysis_path = atlas_path(/analysis)
    
    if keyword_set(nuclear) then begin
       root = 'nuclear_atlas'
    endif else begin
       root = 'integrated_atlas'
    endelse

    snrcut_linedust = 3.0
    outfile = root+'_speclinefit_'+version+'.fits'
    outfile_nodust = repstr(outfile,'.fits','_nodust.fits')

    atlas = atlas_read_info()

; read the latest specdatafile    
    
    specfitfiles = file_search(djs_filepath('?????_'+root+'_specdata.fits.gz'),count=fcount)
    specfitfile = specfitfiles[(reverse(sort(specfitfiles)))[0]]

    splog, 'Reading '+specfitfile
    line = mrdfits(specfitfile,1)

    if keyword_set(nuclear) then begin

       match, strtrim(atlas.nuclear_file,2), strtrim(line.specfile,2), m1, m2
       atlas_final = struct_addtags(atlas[m1],struct_trimtags(line[m2],except='GALAXY'))
       atlas_final_nodust = iunred_linedust(atlas_final,snrcut=snrcut_linedust)

       splog, 'Writing '+outpath+outfile
       mwrfits, atlas_final, outpath+outfile, /create
       spawn, 'gzip -f '+outpath+outfile
       
       splog, 'Writing '+outpath+outfile_nodust
       mwrfits, atlas_final_nodust, outpath+outfile_nodust, /create
       spawn, 'gzip -f '+outpath+outfile_nodust

    endif else begin

       match, strtrim(atlas.drift_file,2), strtrim(line.specfile,2), m1, m2
       atlas_final = struct_addtags(atlas[m1],struct_trimtags(line[m2],except='GALAXY'))
       atlas_final_nodust = iunred_linedust(atlas_final,snrcut=snrcut_linedust)
       
       splog, 'Writing '+outpath+outfile
       mwrfits, atlas_final, outpath+outfile, /create
       spawn, 'gzip -f '+outpath+outfile
       
       splog, 'Writing '+outpath+outfile_nodust
       mwrfits, atlas_final_nodust, outpath+outfile_nodust, /create
       spawn, 'gzip -f '+outpath+outfile_nodust
       
    endelse
    
return
end

; ALL CODE BELOW HERE RELEGATED FOR NOW!! jm08apr14nyu    
    
;   magfloor_err = 0.05
    
;;    synthfile = root+'_synthmags.fits.gz'
;;    if file_test(datapath+synthfile,/regular) then begin
;;       splog, 'Reading and appending '+datapath+synthfile+'.'
;;       mags = mrdfits(datapath+synthfile,1,/silent)
;;       match, strtrim(atlas.galaxy,2), strtrim(mags.galaxy,2), atlasindx, magsindx
;;       atlas = struct_addtags(atlas[atlasindx],struct_trimtags(mags[magsindx],except=['GALAXY','DISTANCE','DISTANCE_ERR']))
;;    endif
;;
;;    massfile = 'atlas_stellar_mass_photo.fits.gz'
;;    if file_test(analysis_path+massfile,/regular)then begin
;;       splog, 'Reading and appending '+analysis_path+massfile+'.'
;;       mass = mrdfits(analysis_path+massfile,1,/silent)
;;       match, strtrim(atlas.galaxy,2), strtrim(mass.galaxy,2), atlasindx, massindx
;;       atlas = struct_addtags(atlas[atlasindx],struct_trimtags(mass[massindx],except='GALAXY'))
;;    endif
;;
;;    kcorrmassfile = 'atlas_stellar_mass_kcorrect.fits.gz'
;;    if file_test(analysis_path+kcorrmassfile,/regular)then begin
;;       splog, 'Reading and appending '+analysis_path+kcorrmassfile+'.'
;;       kcorrmass = mrdfits(analysis_path+kcorrmassfile,1,/silent)
;;       match, strtrim(atlas.galaxy,2), strtrim(kcorrmass.galaxy,2), atlasindx, kcorrmassindx
;;       atlas = struct_addtags(atlas[atlasindx],struct_trimtags(kcorrmass[kcorrmassindx],except='GALAXY'))
;;    endif

;;; finalize the UBVR photometry; apply an aperture correction to the
;;; synthesized B magnitudes, based on the RC3 B-band magnitude (do this
;;; because our errors are *much* smaller); use the synthesized U-B,
;;; B-V, and B-R colors to get the other bandpasses; correct the
;;; photometry for emission lines but do not let the correction be less
;;; than zero; use synthesized B-band photometry where no RC3
;;; measurements are available; we could also correct for internal
;;; extinction here, or to face-on inclination
;;
;;    if file_test(datapath+synthfile,/regular) then begin
;;
;;; just tabulate the apparent magnitudes and then recompute the
;;; absolute magnitudes and luminosities as in WRITE_ANCILLARY_DATA        
;;       
;;       ucor = (atlas.synth_u-atlas.synth_u_obs) > 0.0 ; <-- NOTE!
;;       bcor = (atlas.synth_b-atlas.synth_b_obs) > 0.0 ; <-- NOTE!
;;       vcor = (atlas.synth_v-atlas.synth_v_obs) > 0.0 ; <-- NOTE!
;;       rcor = (atlas.synth_r-atlas.synth_r_obs) > 0.0 ; <-- NOTE!
;;
;;;      ustats = im_stats(ucor,/verbose,/baremin)
;;;      bstats = im_stats(bcor,/verbose,/baremin,/no_head)
;;;      vstats = im_stats(vcor,/verbose,/baremin,/no_head)
;;;      rstats = im_stats(rcor,/verbose,/baremin,/no_head)
;;
;;; B
;;       
;;       rc3 = where((atlas.rc3_b gt -900.0) and (atlas.synth_b lt -900.0),nrc3)
;;       if (nrc3 ne 0L) then begin
;;          atlas[rc3].b_obs     = atlas[rc3].rc3_b
;;          atlas[rc3].b_obs_err = atlas[rc3].rc3_b_err
;;       endif
;;
;;       synth = where((atlas.rc3_b lt -900.0) and (atlas.synth_b gt -900.0),nsynth)
;;       if (nsynth ne 0L) then begin
;;          atlas[synth].b_obs     = atlas[synth].synth_b_obs
;;          atlas[synth].b         = atlas[synth].synth_b
;;          atlas[synth].b_obs_err = sqrt(atlas[synth].synth_b_obs_err^2 + magfloor_err^2)
;;          atlas[synth].b_err     = sqrt(atlas[synth].synth_b_err^2 + magfloor_err^2)
;;       endif
;;
;;       good = where((atlas.rc3_b gt -900.0) and (atlas.synth_b gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;
;;          bapercor = atlas[good].rc3_b-atlas[good].synth_b_obs
;;          
;;          atlas[good].b_obs     = atlas[good].synth_b_obs + bapercor
;;          atlas[good].b_obs_err = sqrt(atlas[good].synth_b_obs_err^2 + magfloor_err^2)
;;
;;          atlas[good].b         = atlas[good].b_obs + bcor[good] ; correct for emission lines (make fainter)
;;          atlas[good].b_err     = atlas[good].b_obs_err
;;
;;       endif
;;;      niceprint, atlas.galaxy, atlas.b_obs, atlas.b_obs_err, atlas.rc3_b, atlas.rc3_b_err
;;
;;; U       
;;       
;;       good = where((atlas.b gt -900.0) and (atlas.synth_u gt -900.0) and (atlas.synth_b gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;
;;          ub = atlas[good].synth_u-atlas[good].synth_b
;;          ub_err = sqrt(atlas[good].synth_u_err^2+atlas[good].synth_b_err^2)
;;       
;;          atlas[good].u         = atlas[good].b + ub
;;          atlas[good].u_obs     = atlas[good].b + ub - ucor[good] ; add the emission lines back in (make brighter)
;;          atlas[good].u_err     = sqrt(atlas[good].b_err^2 + ub_err^2)
;;          atlas[good].u_obs_err = sqrt(atlas[good].b_err^2 + ub_err^2)
;;
;;;         w = where(atlas.u_obs gt -900.0 and atlas.rc3_u gt -900.0)
;;;         plot, atlas[w].u_obs, atlas[w].rc3_u, ps=4, xsty=3, ysty=3, xr=[8,16], yr=[8,16]
;;;         oplot, !x.crange, !y.crange, thick=2
;;;         cc = get_kbrd(1)
;;
;;       endif
;;
;;; V
;;       
;;       good = where((atlas.b gt -900.0) and (atlas.synth_v gt -900.0) and (atlas.synth_b gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;
;;          bv = atlas[good].synth_b-atlas[good].synth_v
;;          bv_err = sqrt(atlas[good].synth_b_err^2+atlas[good].synth_v_err^2)
;;       
;;          atlas[good].v         = atlas[good].b - bv
;;          atlas[good].v_obs     = atlas[good].b - bv - vcor[good] ; add the emission lines back in (make brighter)
;;          atlas[good].v_err     = sqrt(atlas[good].b_err^2 + bv_err^2)
;;          atlas[good].v_obs_err = sqrt(atlas[good].b_err^2 + bv_err^2)
;;
;;;         w = where(atlas.v_obs gt -900.0 and atlas.rc3_v gt -900.0)
;;;         plot, atlas[w].v_obs, atlas[w].rc3_v, ps=4, xsty=3, ysty=3, xr=[8,16], yr=[8,16]
;;;         oplot, !x.crange, !y.crange, thick=2
;;;         cc = get_kbrd(1)
;;
;;       endif
;;
;;; R
;;       
;;       good = where((atlas.b gt -900.0) and (atlas.synth_r gt -900.0) and (atlas.synth_b gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;
;;          br = atlas[good].synth_b-atlas[good].synth_r
;;          br_err = sqrt(atlas[good].synth_b_err^2+atlas[good].synth_r_err^2)
;;       
;;          atlas[good].r         = atlas[good].b - br
;;          atlas[good].r_obs     = atlas[good].b - br - rcor[good] ; add the emission lines back in (make brighter)
;;          atlas[good].r_err     = sqrt(atlas[good].b_err^2 + br_err^2)
;;          atlas[good].r_obs_err = sqrt(atlas[good].b_err^2 + br_err^2)
;;
;;       endif
;;
;;    endif
;;
;;;   niceprint, atlas.galaxy, atlas.u, atlas.u_obs, atlas.b, atlas.b_obs, atlas.v, atlas.v_obs, atlas.r, atlas.r_obs
;;;   niceprint, atlas.galaxy, atlas.u_err, atlas.u_obs_err, atlas.b_err, atlas.b_obs_err, atlas.v_err, atlas.v_obs_err, atlas.r_err, atlas.r_obs_err
;;    
;;; compute luminosities and absolute magnitudes
;;    
;;    splog, 'Computing absolute magnitudes and luminosities.'
;;    band = ['U','B','V','R','U','B','V','R']
;;
;;    tags = [['u','b','v','r'],['u','b','v','r']+'_obs']
;;    tags_err = tags+'_err'
;;
;;    abstags = [['m_u','m_b','m_v','m_r'],['m_u','m_b','m_v','m_r']+'_obs']
;;    abstags_err = abstags+'_err'
;;    
;;    lumtags = [['u','b','v','r']+'_lum',['u','b','v','r']+'_lum_obs']
;;    lumtags_err = lumtags+'_err'
;;
;;    for iband = 0L, n_elements(tags)-1L do begin
;;
;;       true = tag_exist(atlas,tags[iband],index=tagsindx)
;;       true = tag_exist(atlas,tags_err[iband],index=tagsindx_err)
;;       true = tag_exist(atlas,abstags[iband],index=abstagsindx)
;;       true = tag_exist(atlas,abstags_err[iband],index=abstagsindx_err)
;;       true = tag_exist(atlas,lumtags[iband],index=lumtagsindx)
;;       true = tag_exist(atlas,lumtags_err[iband],index=lumtagsindx_err)
;;
;;       good = where((atlas.distance gt -900.0) and (atlas.(tagsindx) gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;
;;          mags = im_absolute_magnitudes(band[iband],atlas[good].(tagsindx),$
;;            atlas[good].distance,mag_err=atlas[good].(tagsindx_err),$
;;            distance_err=atlas[good].distance_err)
;;
;;          atlas[good].(abstagsindx)     = mags.absmag
;;          atlas[good].(abstagsindx_err) = mags.absmag_err
;;          
;;          atlas[good].(lumtagsindx)     = mags.lum
;;          atlas[good].(lumtagsindx_err) = mags.lum_err
;;          
;;       endif
;;       
;;    endfor

;;; write out    
;;    
;;    select_lines = ['OII_3727','H_BETA','OIII_5007','H_ALPHA']
;;    disttag = 'DISTANCE'
;;    disterrtag = 'DISTANCE_ERR'
;;
;;    snrcut_linedust = 3.0
;;    snrcut_abundance = 3.0
;;
;;    syserr = 0.0 ; think *carefully* about using a non-zero value!!
;;;   syserr = 4.0 ; systematic error [%]
;;    
;;; parse everything 
;;
;;    line = parse_ispeclinefit(datapath=datapath,outpath=outpath,prepend=atlas,$
;;      root=root,trimtags='GALAXY',select_lines=select_lines,syserr=syserr,$
;;      disttag=disttag,disterrtag=disterrtag,photerrtag=photerrtag,$
;;      outfile=outfile,electrondensity=1,snrcut_linedust=snrcut_linedust,$
;;      snrcut_abundance=snrcut_abundance,/match,/kauffmann,/write,$
;;      /odonnell,/nopropagate,_extra=extra)

;; ##################################################
;; OLD CODE WHICH RELIES ON THE RC3 PHOTOMETRY
;; ##################################################
;;
;;    if file_test(datapath+synthfile,/regular) then begin
;;       
;;; U
;;
;;       cor = (atlas.synth_u-atlas.synth_u_obs) > 0.0 ; <-- NOTE!
;;       good = where(atlas.rc3_u gt -900.0,ngood,comp=synth,ncomp=nsynth)
;;       if (ngood ne 0L) then begin
;;
;;          atlas[good].u_obs         = atlas[good].rc3_u
;;          atlas[good].u_obs_err     = atlas[good].rc3_u_err
;;          atlas[good].M_u_obs       = atlas[good].rc3_M_u
;;          atlas[good].M_u_obs_err   = atlas[good].rc3_M_u_err
;;          atlas[good].u_lum_obs     = atlas[good].rc3_u_lum
;;          atlas[good].u_lum_obs_err = atlas[good].rc3_u_lum_err
;;
;;          atlas[good].u         = atlas[good].rc3_u + cor[good]
;;          atlas[good].u_err     = atlas[good].rc3_u_err
;;          atlas[good].M_u       = atlas[good].rc3_M_u + cor[good]
;;          atlas[good].M_u_err   = atlas[good].rc3_M_u_err
;;          atlas[good].u_lum     = atlas[good].rc3_u_lum - 0.4*cor[good]
;;          atlas[good].u_lum_err = atlas[good].rc3_u_lum_err
;;       endif
;;       if (nsynth ne 0L) then begin
;;
;;          atlas[synth].u_obs         = atlas[synth].synth_u_obs
;;          atlas[synth].u_obs_err     = atlas[synth].synth_u_obs_err
;;          atlas[synth].M_u_obs       = atlas[synth].synth_M_u + cor[synth]
;;          atlas[synth].M_u_obs_err   = atlas[synth].synth_M_u_err
;;          atlas[synth].u_lum_obs     = atlas[synth].synth_u_lum - 0.4*cor[synth]
;;          atlas[synth].u_lum_obs_err = atlas[synth].synth_u_lum_err
;;
;;          atlas[synth].u         = atlas[synth].synth_u
;;          atlas[synth].u_err     = atlas[synth].synth_u_err
;;          atlas[synth].M_u       = atlas[synth].synth_M_u
;;          atlas[synth].M_u_err   = atlas[synth].synth_M_u_err
;;          atlas[synth].u_lum     = atlas[synth].synth_u_lum
;;          atlas[synth].u_lum_err = atlas[synth].synth_u_lum_err
;;       endif
;;
;;; B
;;       
;;       cor = (atlas.synth_b-atlas.synth_b_obs) > 0.0 ; <-- NOTE!
;;       good = where(atlas.rc3_b gt -900.0,ngood,comp=synth,ncomp=nsynth)
;;       if (ngood ne 0L) then begin
;;          atlas[good].b_obs         = atlas[good].rc3_b
;;          atlas[good].b_obs_err     = atlas[good].rc3_b_err
;;          atlas[good].M_b_obs       = atlas[good].rc3_M_b
;;          atlas[good].M_b_obs_err   = atlas[good].rc3_M_b_err
;;          atlas[good].b_lum_obs     = atlas[good].rc3_b_lum
;;          atlas[good].b_lum_obs_err = atlas[good].rc3_b_lum_err
;;
;;          atlas[good].b         = atlas[good].rc3_b + cor[good]
;;          atlas[good].b_err     = atlas[good].rc3_b_err
;;          atlas[good].M_b       = atlas[good].rc3_M_b + cor[good]
;;          atlas[good].M_b_err   = atlas[good].rc3_M_b_err
;;          atlas[good].b_lum     = atlas[good].rc3_b_lum - 0.4*cor[good]
;;          atlas[good].b_lum_err = atlas[good].rc3_b_lum_err
;;       endif
;;       if (nsynth ne 0L) then begin
;;          atlas[synth].b_obs         = atlas[synth].synth_b_obs
;;          atlas[synth].b_obs_err     = atlas[synth].synth_b_obs_err
;;          atlas[synth].M_b_obs       = atlas[synth].synth_M_b + cor[synth]
;;          atlas[synth].M_b_obs_err   = atlas[synth].synth_M_b_err
;;          atlas[synth].b_lum_obs     = atlas[synth].synth_b_lum - 0.4*cor[synth]
;;          atlas[synth].b_lum_obs_err = atlas[synth].synth_b_lum_err
;;
;;          atlas[synth].b         = atlas[synth].synth_b
;;          atlas[synth].b_err     = atlas[synth].synth_b_err
;;          atlas[synth].M_b       = atlas[synth].synth_M_b
;;          atlas[synth].M_b_err   = atlas[synth].synth_M_b_err
;;          atlas[synth].b_lum     = atlas[synth].synth_b_lum
;;          atlas[synth].b_lum_err = atlas[synth].synth_b_lum_err
;;       endif
;;
;;; V    
;;       
;;       cor = (atlas.synth_v-atlas.synth_v_obs) > 0.0 ; <-- NOTE!
;;       good = where(atlas.rc3_v gt -900.0,ngood,comp=synth,ncomp=nsynth)
;;       if (ngood ne 0L) then begin
;;          atlas[good].v_obs         = atlas[good].rc3_v
;;          atlas[good].v_obs_err     = atlas[good].rc3_v_err
;;          atlas[good].M_v_obs       = atlas[good].rc3_M_v
;;          atlas[good].M_v_obs_err   = atlas[good].rc3_M_v_err
;;          atlas[good].v_lum_obs     = atlas[good].rc3_v_lum
;;          atlas[good].v_lum_obs_err = atlas[good].rc3_v_lum_err
;;
;;          atlas[good].v         = atlas[good].rc3_v + cor[good]
;;          atlas[good].v_err     = atlas[good].rc3_v_err
;;          atlas[good].M_v       = atlas[good].rc3_M_v + cor[good]
;;          atlas[good].M_v_err   = atlas[good].rc3_M_v_err
;;          atlas[good].v_lum     = atlas[good].rc3_v_lum - 0.4*cor[good]
;;          atlas[good].v_lum_err = atlas[good].rc3_v_lum_err
;;       endif
;;       if (nsynth ne 0L) then begin
;;          atlas[synth].v_obs         = atlas[synth].synth_v_obs
;;          atlas[synth].v_obs_err     = atlas[synth].synth_v_obs_err
;;          atlas[synth].M_v_obs       = atlas[synth].synth_M_v + cor[synth]
;;          atlas[synth].M_v_obs_err   = atlas[synth].synth_M_v_err
;;          atlas[synth].v_lum_obs     = atlas[synth].synth_v_lum - 0.4*cor[synth]
;;          atlas[synth].v_lum_obs_err = atlas[synth].synth_v_lum_err
;;
;;          atlas[synth].v         = atlas[synth].synth_v
;;          atlas[synth].v_err     = atlas[synth].synth_v_err
;;          atlas[synth].M_v       = atlas[synth].synth_M_v
;;          atlas[synth].M_v_err   = atlas[synth].synth_M_v_err
;;          atlas[synth].v_lum     = atlas[synth].synth_v_lum
;;          atlas[synth].v_lum_err = atlas[synth].synth_v_lum_err
;;       endif
;;
;;    endif
       
