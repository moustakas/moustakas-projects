pro write_03gildepaz, data
; jm05mar01uofa

    root = '03gildepaz'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    table1 = im_read_fmr(path+root+'_table1.dat')
    table5 = im_read_fmr(path+root+'_table5.dat')
    ngalaxy = n_elements(table5)
    
    ra = string(table1.rah,format='(I2.2)')+':'+string(table1.ram,format='(I2.2)')+':'+$
      repstr(string(table1.ras,format='(F4.1)'),' ','0')
    dec = table1.de_+string(table1.ded,format='(I2.2)')+':'+string(table1.dem,format='(I2.2)')+':'+$
      repstr(string(table1.des,format='(F4.1)'),' ','0')

    moredata = replicate({galaxy: '', ra: '', dec: ''},ngalaxy)
    moredata.galaxy = strtrim(table1.name,2)
    moredata.ra = ra
    moredata.dec = dec
    
    data = struct_addtags(moredata,struct_trimtags(table5,except=['NAME']))

; "uncorrect" the B magnitudes for Galactic extinction

;   data.bmag = data.bmag + table1.AB
    
; write out    

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, data, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh
    
return
end    
