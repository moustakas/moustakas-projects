;+
; NAME:
;       NFGS_STELLAR_MASS_PHOTO
;
; PURPOSE:
;       Estimate stellar masses for the NFGS using the color-based
;       methodology. 
;
; COMMENTS:
;       Use the emission-line corrected magnitudes.  
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Aug 04, U of A
;-

pro nfgs_stellar_mass_photo, nfgs, nfgsinfo, mass, write=write

    outpath = nfgs_path(/analysis)
    outfile = 'nfgs_stellar_mass_photo.fits'
    
    if (n_elements(nfgs) eq 0L) then nfgs = read_nfgs()
    if (n_elements(nfgsinfo) eq 0L) then nfgsinfo = nfgs_read_info()
    ngalaxy = n_elements(nfgs)
    nnfgsinfo = n_elements(nfgsinfo)

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

    bigphoto = replicate(photo,nnfgsinfo)
    photo = replicate(photo,ngalaxy)

    bigphoto.galaxy = nfgsinfo.galaxy
    photo.galaxy    = nfgs.galaxy

    photo.b         = nfgs.b
    photo.b_err     = nfgs.b_err
    photo.b_lum     = nfgs.b_lum
    photo.b_lum_err = nfgs.b_lum_err
    photo.v         = nfgs.v
    photo.v_err     = nfgs.v_err
    photo.v_lum     = nfgs.v_lum
    photo.v_lum_err = nfgs.v_lum_err
    photo.r         = nfgs.r
    photo.r_err     = nfgs.r_err
    photo.r_lum     = nfgs.r_lum
    photo.r_lum_err = nfgs.r_lum_err

    photo.h         = nfgs.twomass_h
    photo.h_err     = nfgs.twomass_h_err
    photo.h_lum     = nfgs.twomass_h_lum
    photo.h_lum_err = nfgs.twomass_h_lum_err
    photo.k         = nfgs.twomass_ks ; <-- NEED a color term!
    photo.k_err     = nfgs.twomass_ks_err
    photo.k_lum     = nfgs.twomass_ks_lum
    photo.k_lum_err = nfgs.twomass_ks_lum_err

; copy "photo" into "bigphoto"

;   goodindx = cmset_op(nfgsinfo.nfgs_id,'AND',nfgs.nfgs_id,/index)
    doit = match_string(photo.galaxy,bigphoto.galaxy,/exact,findex=index)
;   niceprint, photo.galaxy, bigphoto[index].galaxy

; compute stellar masses
    
    mass = compute_stellar_mass(bigphoto)
    mass = struct_addtags(replicate({galaxy: ''},nnfgsinfo),mass)
    mass.galaxy = nfgsinfo.galaxy

    if keyword_set(write) then begin

       splog, 'Writing '+outpath+outfile+'.'
       mwrfits, mass, outpath+outfile, /create
       spawn, ['gzip -f '+outpath+outfile], /sh

    endif
    
return
end
