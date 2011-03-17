;+
; NAME:
;       ATLAS_STELLAR_MASS_PHOTO
;
; PURPOSE:
;       Estimate stellar masses for the spectral atlas using the Bell
;       et al. color-based methodology.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;       Use the emission-line corrected magnitudes with which to
;       compute stellar masses.  
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 May 17, U of A
;-

pro atlas_stellar_mass_photo, atlas, atlasinfo, mass, write=write

    outpath = atlas_path(/analysis)
    outfile = 'atlas_stellar_mass_photo.fits'
    
    if (n_elements(atlas) eq 0L) then atlas = read_integrated()
    if (n_elements(atlasinfo) eq 0L) then atlasinfo = atlas_read_info()
    ngalaxy = n_elements(atlas)
    natlasinfo = n_elements(atlasinfo)

    photo = {$
      galaxy:        '', $

      B:         -999.0, $
      B_err:     -999.0, $
      B_lum:     -999.0, $
      B_lum_err: -999.0, $
      V:         -999.0, $
      V_err:     -999.0, $
      V_lum:     -999.0, $
      V_lum_err: -999.0, $
      R:         -999.0, $
      R_err:     -999.0, $
      R_lum:     -999.0, $
      R_lum_err: -999.0, $
      H:         -999.0, $
      H_err:     -999.0, $
      H_lum:     -999.0, $
      H_lum_err: -999.0, $
      K:         -999.0, $
      K_err:     -999.0, $
      K_lum:     -999.0, $
      K_lum_err: -999.0}

    bigphoto = replicate(photo,natlasinfo)
    photo = replicate(photo,ngalaxy)

    bigphoto.galaxy = atlasinfo.galaxy
    photo.galaxy    = atlas.galaxy

    photo.b         = atlas.b
    photo.b_err     = atlas.b_err
    photo.b_lum     = atlas.b_lum
    photo.b_lum_err = atlas.b_lum_err
    photo.v         = atlas.v
    photo.v_err     = atlas.v_err
    photo.v_lum     = atlas.v_lum
    photo.v_lum_err = atlas.v_lum_err
    photo.r         = atlas.r
    photo.r_err     = atlas.r_err
    photo.r_lum     = atlas.r_lum
    photo.r_lum_err = atlas.r_lum_err

    photo.h         = atlas.twomass_h
    photo.h_err     = atlas.twomass_h_err
    photo.h_lum     = atlas.twomass_h_lum
    photo.h_lum_err = atlas.twomass_h_lum_err
    photo.k         = atlas.twomass_ks ; <-- NEED a color term!
    photo.k_err     = atlas.twomass_ks_err
    photo.k_lum     = atlas.twomass_ks_lum
    photo.k_lum_err = atlas.twomass_ks_lum_err

; copy "photo" into "bigphoto"

;   goodindx = cmset_op(atlasinfo.atlas_id,'AND',atlas.atlas_id,/index)
    doit = match_string(photo.galaxy,bigphoto.galaxy,/exact,findex=index)
;   niceprint, photo.galaxy, bigphoto[index].galaxy

    bigphoto[index] = photo

; compute stellar masses
    
    mass = compute_stellar_mass(bigphoto)
    mass = struct_addtags(replicate({galaxy: ''},natlasinfo),mass)
    mass.galaxy = atlasinfo.galaxy

    if keyword_set(write) then begin

       splog, 'Writing '+outpath+outfile+'.'
       mwrfits, mass, outpath+outfile, /create
       spawn, ['gzip -f '+outpath+outfile], /sh

    endif
    
return
end
    
