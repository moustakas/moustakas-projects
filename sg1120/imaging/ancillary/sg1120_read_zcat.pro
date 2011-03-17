function sg1120_read_zcat, suffix=suffix, silent=silent, write_regions=write_regions
; jm07aug22nyu - read the current redshift catalog
; jm08jun09nyu - use the updated redshift catalog
; jm09jun15nyu - updated to read the latest redshift catalog

    zcatpath = sg1120_path(/zcat)
    suffix = 'refR2.0chi2'
    catname_zcat = zcatpath+'sg1120_09jun10_'+suffix+'.zcat'

    if (file_test(catname_zcat,/regular) eq 0L) then begin
       splog, 'Redshift catalog '+catname_zcat+' not found'
       return, -1
    endif
    if (keyword_set(silent) eq 0) then splog, 'Reading '+catname_zcat
    readcol, catname_zcat, runid, rcat_id, z, q, rcat_ra, rcat_dec, $
      format='A,L,F,I,D,D', comment='#', /silent

    ngal = n_elements(runid)
    zcat = {$
      runid:   '', $
      id:      0L, $
      z:       0.0,$
      q:       0,  $
      ra:      0.0D,$
      dec:     0.0D}
    zcat = replicate(zcat,ngal)

    zcat.runid = runid
    zcat.id = rcat_id
    zcat.z = z
    zcat.q = q
    zcat.ra = rcat_ra
    zcat.dec = rcat_dec
    
;   zcatpath = sg1120_path(/mcatalogs)
;   catname_zcat = zcatpath+'sg1120_08jun08.zcat.sex'
;
;   if (file_test(catname_zcat,/regular) eq 0L) then begin
;      splog, 'Redshift catalog '+catname_zcat+' not found'
;      return, -1
;   endif
;   if (not keyword_set(silent)) then splog, 'Reading '+catname_zcat
;   zcat1 = rsex(catname_zcat)
;
;   zcat = im_struct_trimtags(zcat1,select=tag_names(zcat1[0]),$
;     newtags='ZCAT_'+tag_names(zcat1[0]))
;
;   if keyword_set(write_regions) then begin
;      regionfile = repstr(catname_zcat,'.sex','.radec.reg')
;      splog, 'Writing '+regionfile
;      good = where(zcat.zcat_ra gt -900.0)
;      write_ds9_regionfile, zcat[good].zcat_ra, zcat[good].zcat_dec, $
;        filename=regionfile, color='blue', symbol='box'
;   endif
    
return, zcat
end
