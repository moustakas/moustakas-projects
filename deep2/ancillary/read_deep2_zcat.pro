function read_deep2_zcat, dr=dr, all=all, good=good
; jm06aug28uofa

    zcatpath = deep2_path(/analysis)
    if (n_elements(dr) eq 0L) then dr = 'dr3'

    version = 'v1_0' ; this is a DEEP2 internal version number
    if keyword_set(all) then suffix = '.'+version else suffix = '.'+version+'.uniq'
    if keyword_set(good) then suffix = '.uniq.good'

    zcatfile = 'zcat.'+dr+suffix+'.fits.gz'

    if (file_test(zcatpath+zcatfile) eq 0L) then begin
       splog, 'Redshift catalog '+zcatpath+zcatfile+' not found.'
       return, -1L
    endif

    zcat = mrdfits(zcatpath+zcatfile,1,/silent)
    if keyword_set(good) then zcat = zcat[where((zcat.z gt 0.0) and $
      (zcat.zquality ge 3L))]

return, zcat
end
