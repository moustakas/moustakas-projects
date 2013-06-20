function read_deep2_zcat, dr=dr, all=all, good=good
; jm06aug28uofa
; jm13jun17siena - updated to DR4

    zcatpath = deep2_path(/catalogs)
    if (n_elements(dr) eq 0L) then dr = 'dr4'

    version = '' ; dr4
;   version = 'v1_0' ; this is a DEEP2 internal version number for dr3
    if keyword_set(all) then suffix = version else suffix = version+'uniq'
    if keyword_set(good) then suffix = '.uniq.good'

    zcatfile = 'zcat.'+dr+suffix+'.fits.gz' ; dr3 
;   zcatfile = 'zcat.deep2.'+dr+'.'+suffix+'.fits.gz' ; dr4

    if (file_test(zcatpath+zcatfile) eq 0L) then begin
       splog, 'Redshift catalog '+zcatpath+zcatfile+' not found.'
       return, -1L
    endif

    splog, 'Reading '+zcatpath+zcatfile
    zcat = mrdfits(zcatpath+zcatfile,1,/silent)
    if keyword_set(good) then zcat = zcat[where((zcat.z gt 0.0) and $
      (zcat.zquality ge 3L))]

return, zcat
end
