function read_zcat, file, type=type
    readcol, file, ra, dec, z, zerr, r, flag, type, $
      format='D,D,F,F,F,I,A', /silent, comment='#'
    ngal = n_elements(ra)
    zcat = replicate({ra: 0D, dec: 0D, z: 0.0, $
      zerr: 0.0, flag: 0},ngal)
    zcat.ra = ra
    zcat.dec = dec
    zcat.z = z
    zcat.zerr = zerr
;   zcat.r = r
    zcat.flag = flag
;   zcat.type = type
return, zcat
end

function read_clash_redshifts
; jm13dec08siena - read the current compilation of redshifts for *all*
; the CLASH fields; currently no cuts on quality/flag
    path = getenv('IM_DATA_DIR')+'/clash/'
    zcat = read_zcat(path+'clash_redshifts_25apr2013.cat',type=type)
    keep = where(strtrim(type,2) eq 'Gal')
    zcat = zcat[keep]
return, zcat
end
