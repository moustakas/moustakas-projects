function read_deep2_zcat, dr=dr, all=all
; jm06aug28uofa
; jm13jun17siena - updated to DR4
; jm13jul14siena - read all the unique good redshifts with good 1D
;   spectra by default (as determined by DEEP2_CHECK_SPEC1D_DR4),
;   otherwise read the unique redshift catalog packaged by Cooper on
;   the DR4 website 

    zcatpath = deep2_path(/catalogs)
    if (n_elements(dr) eq 0L) then dr = 'dr4'

    version = '' ; dr4
;   version = 'v1_0' ; this is a DEEP2 internal version number for dr3

    suffix = version+'.uniq.good' ; default
    if keyword_set(all) then suffix = version

    zcatfile = 'zcat.'+dr+suffix+'.fits.gz' ; dr3 
;   zcatfile = 'zcat.deep2.'+dr+'.'+suffix+'.fits.gz' ; dr4

    if (file_test(zcatpath+zcatfile) eq 0L) then begin
       splog, 'Redshift catalog '+zcatpath+zcatfile+' not found.'
       return, -1
    endif

    splog, 'Reading '+zcatpath+zcatfile
    zcat = mrdfits(zcatpath+zcatfile,1,/silent)

; keep just the good ones    
    if keyword_set(all) eq 0 then zcat = zcat[where((zcat.z gt 0.0) and $
      (zcat.zquality ge 3L))]

return, zcat
end
