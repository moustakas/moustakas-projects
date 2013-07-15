function read_deep2_zcat, dr=dr, all=all, photo=photo
; jm06aug28uofa
; jm13jun17siena - updated to DR4
; jm13jul14siena - by default read the sample of objects with (Q>=3)
;   redshifts *and* with good 1D spectra (as determined by
;   DEEP2_CHECK_SPEC1D_DR4), unless /ALL is set in which read
;   even the low-quality redshift objects

    catpath = deep2_path(/catalogs)
    if (n_elements(dr) eq 0) then dr = 'dr4'

    if keyword_set(all) then suffix = '' else suffix = '.Q34'
    zcatfile = 'zcat.'+dr+'.goodspec1d'+suffix+'.fits.gz'

    if (file_test(catpath+zcatfile) eq 0L) then begin
       splog, 'Redshift catalog '+catpath+zcatfile+' not found.'
       return, -1
    endif

    splog, 'Reading '+catpath+zcatfile
    zcat = mrdfits(catpath+zcatfile,1,/silent)

; also optionally read in the line-matched photometry file
    if arg_present(photo) then begin
       photofile = 'photo.'+dr+'.goodspec1d'+suffix+'.fits.gz'
       splog, 'Reading '+catpath+photofile
       photo = mrdfits(catpath+photofile,1,/silent)
    endif
    
return, zcat
end
