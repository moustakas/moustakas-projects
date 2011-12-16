function addphot, in, all
    filt = clash_filterlist(short=short,zpt=zpt)
    out = in
    nall = n_elements(all)
    for ii = 0, n_elements(short)-1 do begin
       gd = where(finite(all.(tag_indx(all,short[ii]+'_fluxerr'))))
       totflux = total(all[gd].(tag_indx(all[gd],short[ii]+'_flux')))
       totferr = sqrt(total(all[gd].(tag_indx(all[gd],short[ii]+'_fluxerr'))^2))
       if (totflux le 0) then begin
          out.(tag_indx(out,short[ii]+'_mag')) = 99.0 
          out.(tag_indx(out,short[ii]+'_magerr')) = -2.5*alog10(totferr)+zpt[ii]
       endif else begin
          out.(tag_indx(out,short[ii]+'_mag')) = -2.5*alog10(totflux)+zpt[ii]
          out.(tag_indx(out,short[ii]+'_magerr')) = 2.5*totferr/totflux/alog(10)
       endelse
    endfor
return, out
end

function read_z9arc
; jm11dec14ucsd - read the z=9 photometry

    isedpath = clash_path(/ised)
    datapath = clash_path(/z9arc)

    filt = clash_filterlist(short=short,weff=weff,zpt=zpt,nice=nice,/useirac)
    nfilt = n_elements(filt)

    readcol, datapath+'test.cat', ap, f160w, f160werr, f140w, f140werr, f125w, f125werr, f110w, f110werr, $
      f105w, f105werr, f850lp, f850lperr, f814w, f814werr, f775w, f775werr, f625w, f625werr, f606w, f606werr, $
      f555w, f555werr, f475w, f475werr, f435w, f435werr, f390w, f390werr, f336w, f336werr, f275w, f275werr, $
      f225w, f225werr, ch1, ch1err, ch2, ch2err, skip=1, /silent
    ngal = n_elements(ap)
    
    mag = [[f225w],[f275w],[f336w],[f390w],[f435w],[f475w],[f606w],[f625w],$
      [f775w],[f814w],[f850lp],[f105w],[f110w],[f125w],[f140w],[f160w],[ch1],[ch2]]
    magerr = [[f225werr],[f275werr],[f336werr],[f390werr],[f435werr],[f475werr],[f606werr],[f625werr],$
      [f775werr],[f814werr],[f850lperr],[f105werr],[f110werr],[f125werr],[f140werr],[f160werr],[ch1err],[ch2err]]
    
    cat = {ap: 0.0}
    for ii = 0, nfilt-1 do cat = create_struct(cat,$
      short[ii]+'_mag',-99.0,short[ii]+'_magerr',-99.0)
    cat = replicate(cat,ngal)
    cat.ap = ap
       
    for ii = 0, nfilt-1 do begin
       magindx = tag_indx(cat,short[ii]+'_mag')
       magerrindx = tag_indx(cat,short[ii]+'_magerr')
       cat.(magindx) = mag[*,ii]
       cat.(magerrindx) = magerr[*,ii]
    endfor

return, cat
end
