pro build_deep2_unwise, query=query, parse=parse
; jm14may21siena - merge the individual unWISE pcat files and match
; them to the extended photometric catalog

    dr = 'dr4'
    
    catpath = deep2_path(/catalogs)
    wisepath = catpath+'unwise/'
    ff = file_search(wisepath+'*.fits.gz',count=nff)
    wise = mrdfits(ff[0],1)
    for ii = 1, nff-1 do wise = [wise,mrdfits(ff[ii],1)]

    phot = mrdfits(catpath+'zcat_ext.uniq.fits.gz',1)
    ngal = n_elements(phot)

    match, phot.objno, wise.objno, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]

; there are 58 objects missing from this list, all because they have
; R>24.1, but only 21 have Q>=3
    missing = lindgen(n_elements(phot))
    remove, m1, missing

    out = im_empty_structure(wise,ncopies=ngal)
    out[missing] = im_struct_assign(phot[missing],out[missing])
    out[m1] = wise[m2]

    im_mwrfits, out, catpath+'deep2.'+dr+'.unwise.fits', /clobber

return
end    
