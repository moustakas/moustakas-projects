pro write_jacoby_fits, jacoby, write=write
; jm02apr5uofa
; jm03may12uofa - add radial velocity data
; jm08jul24nyu - added E(B-V), U_B_obs, and B_V_obs to output
;                structure 
    
    jacoby_path = filepath('',root_dir=getenv('CATALOGS_DIR'),$
      subdirectory='84jacoby')
    atlas_path = jacoby_path+'atlas/'
    
    flist = file_search(atlas_path+'jhc*.fits',count=fcount)

    flux = mrdfits(flist[0],0,h,/silent)
    flux = float(flux) + abs(min(flux)) & flux = flux/max(flux)
    nflux = n_elements(flux)

    w0 = 3510.0D
    wpc = 1.4D
    wave = w0+dindgen(nflux)*wpc
    
;   wave = sxpar(h,'W0')+lindgen(nflux)*sxpar(h,'WPC')
    
    jacoby = {$
      jhc_id:                 0L, $
      star:                   '', $
      sp_type:                '', $
      sp_type_lit:            '', $
      ebv:                -999.0, $
      U_B:                -999.0, $ ; dereddened
      B_V:                -999.0, $ ; dereddened
      U_B_obs:            -999.0, $ ; observed
      B_V_obs:            -999.0, $ ; observed
      rv:                 -999.0, $
      rv_err:             -999.0, $
      flux:        dblarr(nflux), $
      wave:        dblarr(nflux), $
      header:      strarr(n_elements(h)), $
      comment:     ''}
    jacoby = replicate(jacoby,fcount)

    jacoby[0].jhc_id = 1L
    jacoby[0].star = strcompress(sxpar(h,'OBJECT'),/remove)
    jacoby[0].sp_type = strcompress(strmid(h[14],11,10),/remove)
    jacoby[0].ebv = float(strmid(h[16],9,21))
    jacoby[0].U_B = float(strmid(h[19],9,21))
    jacoby[0].B_V = float(strmid(h[20],9,21))
    jacoby[0].U_B_obs = float(strmid(h[17],9,21))
    jacoby[0].B_V_obs = float(strmid(h[18],9,21))
    jacoby[0].flux = flux
    jacoby[0].wave = wave
    jacoby[0].header = h

; read the FITS files    

    for i = 1L, fcount-1L do begin

       jacoby[i].jhc_id = i+1L
       
       flux = mrdfits(flist[i],0,h,/silent)
       flux = float(flux) + abs(min(flux)) & flux = flux/max(flux)

;      w0 = double(sxpar(h,'W0'))
;      wpc = double(sxpar(h,'WPC'))
       wave = w0+dindgen(nflux)*wpc
       
       star = strcompress(sxpar(h,'OBJECT'),/remove)

; fix star names to be SIMBAD-compatible       

       if strmatch(star,'BD*') eq 1B then if strmatch(star,'BD-*') eq 0B then star = repstr(star,'BD','BD+')
       if strmatch(star,'HD17971N',/fold) eq 1B then star = 'HD17971'
       if strmatch(star,'YALE1755',/fold) eq 1B then star = 'BD+051668'
       if strmatch(star,'O1015',/fold) eq 1B then star = 'BD+56521'
       
       jacoby[i].star = star
       jacoby[i].sp_type = strcompress(strmid(h[14],11,10),/remove)
       jacoby[i].ebv = float(strmid(h[16],9,21))
       jacoby[i].U_B = float(strmid(h[19],9,21))
       jacoby[i].B_V = float(strmid(h[20],9,21))
       jacoby[i].U_B_obs = float(strmid(h[17],9,21))
       jacoby[i].B_V_obs = float(strmid(h[18],9,21))
       jacoby[i].flux = flux
       jacoby[i].wave = wave
       jacoby[i].header = h

    endfor

; read the radial velocity data

    readcol, atlas_path+'rv.dat', name, sp_type_lit, rv, rv_err, $
      format='A,X,A,F,F,X', /silent, comment='#'
    doit = match_string(name,jacoby.star,/exact,index=indx)
    jacoby[indx].rv = rv
    jacoby[indx].rv_err = rv_err
    jacoby[indx].sp_type_lit = sp_type_lit    

; read additional simbad info

    readcol, atlas_path+'simbad_info.dat', name, comment, format='A,A', /silent, comment='#'
    comment = repstr(comment,'_',' ')
    doit = match_string(name,jacoby.star,/exact,index=indx)
    jacoby[indx].comment = comment
    
;; generate a SIMBAD query file
;    
;    exclude = ['HZ948','HZ227','TRA14','O2311','LSIVP24','42LSI','HD249240','HD249384']
;    doit = match_string(exclude,jacoby.star,index=indx,/exact)
;    good = lindgen(fcount)
;    remove, indx, good
;
;    openw, lun, 'simbad_query.txt', /get_lun
;    struct_print, struct_trimtags(jacoby[good],select=['STAR']), /no_head, lun=lun
;    free_lun, lun

; the following objects have problems with SIMBAD:  O1015, HZ948,
; HZ227, TRA14, YALE1755, O2311, LSIVP24, 42LSI
    
; read the Jacoby atlas information binary fits file from Vizier and
; output a summary text table

;   info = mrdfits('jacoby_atlas_info.fits',1,hinfo)
;
;   good = where(strcompress(info.x_raj2000,/remove) ne '')
;   info = info[good]
;   srtra = sort(info.x_raj2000)
;   info = info[srtra]

;   print_struct, info, ['NAME','X_RAJ2000','X_DEJ2000','SP','LUM','E_B_V_'], $
;     file=path+'jacoby_atlas_info.txt'

    if keyword_set(write) then begin
    
       outfile = jacoby_path+'jacoby_atlas.fits'
       splog, 'Writing '+outfile+'.'
       mwrfits, jacoby, outfile, /create

    endif
       
return
end
