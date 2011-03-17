pro gto_starbursts_oh, atlas, atlasnodust
; jm05aug03uofa

    if (n_elements(atlas) eq 0L) then atlas = read_integrated(linefitnodust=atlasnodust)

    indx = index_gto_starbursts(atlas,gto,index_gto=index_gto)
    ngto = n_elements(gto)

    atlas1 = atlas[indx]
    atlasnodust1 = atlasnodust[indx]

    oh = im_abundance(atlasnodust1,snrcut=1.0)

    atlas_oh = {$
      galaxy:              '', $
      log12oh:         -999.0, $
      log12oh_err:     -999.0, $
      log12oh_lit:     -999.0, $
      log12oh_lit_err: -999.0, $
      log12oh_lit_ref:  '...'}
;     h_alpha:     fltarr(2), $
;     h_beta:      fltarr(2), $
;     oii_3727:    fltarr(2), $
;     oiii_5007:   fltarr(2), $
;     nii_6584:    fltarr(2)}
    atlas_oh = replicate(atlas_oh,ngto)

;   atlas_oh.galaxy = gto.nedgalaxy
    atlas_oh.galaxy = gto.galaxy
    atlas_oh[index_gto].log12oh         = oh.zstrong_12oh_oiiinii_niiha
    atlas_oh[index_gto].log12oh_err     = oh.zstrong_12oh_oiiinii_niiha_err
    atlas_oh[index_gto].log12oh_lit     = atlasnodust1.lit_log12oh
    atlas_oh[index_gto].log12oh_lit_err = atlasnodust1.lit_log12oh_err
    atlas_oh[index_gto].log12oh_lit_ref = strtrim(atlasnodust1.lit_log12oh_ref,2)
;   atlas_oh[index_gto].h_alpha     = atlasnodust1.h_alpha
;   atlas_oh[index_gto].h_beta      = atlasnodust1.h_beta
;   atlas_oh[index_gto].oiii_5007   = atlasnodust1.oiii_5007
;   atlas_oh[index_gto].nii_6584    = atlasnodust1.nii_6584

    replace = where(atlas_oh[index_gto].log12oh eq atlasnodust1.zstrong_12oh_niiha_pettini)
    atlas_oh[index_gto[replace]].log12oh = atlasnodust1[replace].zstrong_12oh_niiha_moustakas
    atlas_oh[index_gto[replace]].log12oh_err = atlasnodust1[replace].zstrong_12oh_niiha_moustakas_err
    
    struct_print, atlas_oh[index_gto]
    
    openw, lun, 'gto_starbursts_oh.dat', /get_lun
    struct_print, atlas_oh, lun=lun
    free_lun, lun

stop

return
end
