pro write_94deJong, data
; jm05mar01uofa

    root = '94deJong'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    table1 = mrdfits(path+root+'_table1.fits.gz',1,/silent)
    table6 = mrdfits(path+root+'_table6.fits.gz',1,/silent)
    ngalaxy = n_elements(table1)

    moredata = replicate({galaxy: '', alt_galaxy: '', ra: '', dec: ''},ngalaxy)
    moredata.galaxy = strcompress(table1.ugc,/remove)
    moredata.alt_galaxy = strcompress(table1.name,/remove)
    moredata.ra = repstr(strtrim(table1._raj2000,2),' ',':')
    moredata.dec = repstr(strtrim(table1._dej2000,2),' ',':')
    
    data = struct_addtags(moredata,struct_trimtags(table6,except=['UGC']))

; write out    

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, data, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh
    
return
end    
