pro write_03bell, data, write=write
; jm04nov29uofa - read Table A1 from Bell 2003 on the Radio-FIR
;                 correlation and write out as a convenient binary
;                 FITS table

    root = '03bell'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = 'tableA1.dat'

    data = im_read_fmr(path+file,prefix='BELL')
    ngalaxy = n_elements(data)

; read the NED data and prepend everything

    ned = mrdfits(path+root+'_ned.fits.gz',1,/silent)
    bell = struct_addtags(ned,struct_trimtags(data,except='GALAXY'))
    
    if keyword_set(write) then begin

       splog, 'Writing '+path+root+'.fits.'
       mwrfits, bell, path+root+'.fits', /create
       spawn, ['gzip -f ']+path+root+'.fits', /sh
       
    endif

return
end    
