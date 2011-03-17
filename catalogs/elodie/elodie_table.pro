function elodie_table
; jm03may13uofa

    elodie_path = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='elodie')

    readcol, elodie_path+'doc/table1', id, star, type, rv, format='L,A,A,F', /silent
    nstar = n_elements(id)

    info = mrdfits('elodie_stars.fits.gz',1,/silent) ; downloaded from VIZIER
    
    table = {$
      elodie_id: 0L, $
      ra:        '', $
      dec:       '', $
      star:      '', $
      sp_type:   '', $
      rv:    -999.0  $
      }
    table = replicate(table,nstar)

    table.elodie_id = id
    table.star = star
    table.sp_type = type
    table.rv = rv

    match = where(info.num eq table.elodie_id,nmatch)
    table.ra = info[match]._raj2000
    table.dec = info[match]._dej2000
    
return, table
end    
