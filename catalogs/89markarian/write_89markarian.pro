pro write_89markarian, data
; jm05jul25uofa

    root = '89markarian'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    data = mrdfits(path+'table7.fits.gz',1,/silent)
    ngalaxy = n_elements(data)

    moredata = replicate({galaxy: '', ra: '', dec: ''},ngalaxy)
    moredata.galaxy = 'MRK'+string(data.mrk,format='(I4.4)')
    moredata.ra = repstr(strtrim(data._raj2000,2),' ',':')
    moredata.dec = repstr(strtrim(data._dej2000,2),' ',':')
    
    data = struct_addtags(moredata,struct_trimtags(data,$
      except=['_RAJ2000','_DEJ2000','MRK','RA1950','DE1950',$
      'RECNO']))

; write out    

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, data, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh
    
return
end    
