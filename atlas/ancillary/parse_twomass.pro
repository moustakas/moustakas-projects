pro parse_twomass, write=write
; jm03mar7uofa
; 15" search radius with 2MASS
    
    path = atlas_path(/analysis2d)

; read the input to twomass    

    nedpath = atlas_path(/ned)
    outpath = atlas_path(/analysis2d)

    data = mrdfits(nedpath+'atlas2d_basic_data.fits',1)
    ngalaxy = n_elements(data)
 
    nirdata = create_struct('galaxy', '', 'ra', 0.0, 'dec', 0.0, 'ra_2mass', 0.0, $
      'dec_2mass', 0.0, 'dx', -999.0, 'J', -999.0, 'H', -999.0, 'K', -999.0, $
      'J_err', -999.0, 'H_err', -999.0, 'K_err', -999.0)
    nirdata = replicate(nirdata,ngalaxy)

    nirdata.galaxy = data.galaxy
    nirdata.ra = 15.0*im_hms2dec(data.ra)
    nirdata.dec = im_hms2dec(data.dec)
    
; read the data file
    
    readcol, path+'twomass.out', centr_u, dist_x, pang_x, galaxy_u, ra_u, dec_u, $
      ra, dec, j_m, h_m, k_m, j_msig, h_msig, k_msig, format='I,F,F,A,F,F,F,F,A,A,A,A,A,A', $
      /silent, skipline=20

    unique = uniq(centr_u)
    nunique = n_elements(unique)
    
    splog, 'Matched '+strn(nunique)+' unique galaxies.'

    allindx = centr_u[uniq(centr_u)]

    for i = 0L, nunique-1L do begin
    
       ii = where(centr_u eq allindx[i],count)

       if (count eq 0L) then begin
          splog, 'Problem here.'
          stop
       endif

       if (count eq 1L) then begin

          ii = ii[0]
          doit = match_string(galaxy_u[ii],data.galaxy,index=indx)

;         if strtrim(j_m[ii]) eq 'null' then stop          
;         niceprint, galaxy_u[ii], j_m[ii]
          
          nirdata[indx].ra_2mass = ra[ii]
          nirdata[indx].dec_2mass = dec[ii]
          nirdata[indx].dx = dist_x[ii]

          if (float(j_m[ii]) gt -99.0) and (strtrim(j_m[ii]) ne 'null') then nirdata[indx].J = float(j_m[ii])
          if (float(h_m[ii]) gt -99.0) and (strtrim(h_m[ii]) ne 'null') then nirdata[indx].H = float(h_m[ii])
          if (float(k_m[ii]) gt -99.0) and (strtrim(k_m[ii]) ne 'null') then nirdata[indx].K = float(k_m[ii])

          if (float(j_m[ii]) gt -99.0) and (strtrim(j_m[ii]) ne 'null') then nirdata[indx].J_err = float(j_msig[ii])
          if (float(h_m[ii]) gt -99.0) and (strtrim(h_m[ii]) ne 'null') then nirdata[indx].H_err = float(h_msig[ii])
          if (float(k_m[ii]) gt -99.0) and (strtrim(k_m[ii]) ne 'null') then nirdata[indx].K_err = float(k_msig[ii])
          
       endif
       
       if (count gt 1L) then begin

          splog, 'Found multiple matches for galaxy #'+strn(galaxy_u[ii[0]])+'.'
          
       endif

    endfor

    if keyword_set(write) then begin

       splog, 'Writing '+outpath+'atlas2d_2mass.fits.'
       mwrfits, nirdata, outpath+'atlas2d_2mass.fits', /create

    endif

return
end
