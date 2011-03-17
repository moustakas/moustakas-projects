pro write_96tully, data
; jm05mar02uofa

    root = '96tully'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    table2 = mrdfits(path+root+'_table2.fits.gz',1,/silent)
    table4 = mrdfits(path+root+'_table4.fits.gz',1,/silent)
    table5 = mrdfits(path+root+'_table5.fits.gz',1,/silent)
    ngalaxy = n_elements(table2)

    moredata = replicate({galaxy: '', alt_galaxy: '', ra: '', dec: '', B: -999.0, B_err: -999.0},ngalaxy)
    moredata.galaxy = strcompress(table2.name,/remove)
    moredata.alt_galaxy = 'PGC'+strcompress(table2.pgc,/remove)
    moredata.ra = repstr(strtrim(table2._raj2000,2),' ',':')
    moredata.dec = repstr(strtrim(table2._dej2000,2),' ',':')
    
    data = struct_addtags(moredata,struct_trimtags(table5,except=['RECNO','PGC','NAME']))

; retrieve the (total) B magnitude from the raw photometric catalog

    match = where(strmatch(table4.filter,'*B*') eq 1B)
    data.B = table4[match].Btot
    data.B_err = 0.02
    
; write out    

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, data, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh
    
return
end
    
