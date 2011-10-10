;+
; NAME:
;   BUILD_VAGC_GARCHING
;
; PURPOSE:
;   Build the various VAGC-matched Garching catalogs (emission-line
;   strengths, info, etc.)
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Mar 07, UCSD
;
; Copyright (C) 2010, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro build_vagc_garching, sample=sample, letter=letter, poststr=poststr

;   common com_vagc_garching
    
    if (n_elements(sample) eq 0) then sample = 'dr72'
    if (n_elements(letter) eq 0) then letter = 'bsafe'
    if (n_elements(poststr) eq 0) then poststr = '35'
    suffix = sample+letter+poststr

    mpapath = sdss_path(/mpa_dr7)
    mpacatpath = getenv('VAGC_REDUX')+'/garching_dr7/'
    vagcpath = getenv('LSS_REDUX')+'/'+sample+'/'+letter+'/'+poststr+'/'

    mpacat_outfile = vagcpath+'mpacat.'+suffix+'.fits'
    ispec_outfile = vagcpath+'ispecline.'+suffix+'.fits'
    massoh_outfile = vagcpath+'mpamassoh.'+suffix+'.fits'
    sfr_outfile = vagcpath+'totsfr.'+suffix+'.fits'

; read the VAGC/postlss catalog
    post = mrdfits(vagcpath+'post_catalog.'+suffix+'.fits.gz',1)
    ngal = n_elements(post)

; parse the info structure; matching to POST effectively removes
; duplicates (see MATCH_VAGC_GARCHING)
    mpacat = mrdfits(mpacatpath+'vagc_garching_dr7_catalog.fits.gz',1)
    match, post.object_position, mpacat.object_position, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
;   post = 0

;   out_mpacat = im_empty_structure(mpacat[0],empty_value=-999,ncopies=ngal)
;   out_mpacat[m1] = mpacat[m2]
;   im_mwrfits, temporary(out_mpacat), mpacat_outfile, /clobber

; parse the SFR measurements
    sfr = mrdfits(mpapath+'gal_totsfr_dr7_v5_2.fits.gz',1)
    out_sfr = im_empty_structure(sfr[0],empty_value=-999,ncopies=ngal)
    out_sfr[m1] = sfr[m2]
    im_mwrfits, temporary(out_sfr), sfr_outfile, /clobber

stop    
    
; parse the ispec emission-line measurements
    ispec = mrdfits(mpapath+'ispecline_dr7_v5_2.fits.gz',1)
    out_ispec = im_empty_structure(ispec[0],empty_value=-999,ncopies=ngal)
    out_ispec[m1] = ispec[m2]
    im_mwrfits, temporary(out_ispec), ispec_outfile, /clobber

; stellar mass and metallicity    
    massoh = mrdfits(mpapath+'mpamassoh_dr7_v5_2.fits.gz',1)
    out_massoh = im_empty_structure(massoh[0],empty_value=-999,ncopies=ngal)
    out_massoh[m1] = massoh[m2]
    im_mwrfits, temporary(out_massoh), massoh_outfile, /clobber

return
end
    
