function ages_zptoffsets_do_kcorrect, in_redshift, in_maggies, $
  in_ivarmaggies, $
  filterlist=filterlist
; jm11mar25ucsd - simple wrapper to compute K-corrections 

    h100 = 0.7
    vname = 'default' ; 'default.nolines'
    
    ngal = n_elements(in_redshift)
    nfilt = n_elements(filterlist)
    kcorr = {$
      k_zobj:                          -999.0, $
      k_maggies:                fltarr(nfilt), $
      k_ivarmaggies:            fltarr(nfilt), $
      k_bestmaggies:            fltarr(nfilt), $
      k_mass:                          -999.0, $
      k_coeffs:                     fltarr(5), $
      k_chi2:                          -999.0, $

      k_ugriz_absmag_01:         fltarr(5)-999.0,$
      k_ugriz_absmag_ivar_01:    fltarr(5)-999.0,$
      k_ugriz_kcorrect_01:       fltarr(5)-999.0,$
      k_ugriz_absmag_05:         fltarr(5)-999.0,$
      k_ugriz_absmag_ivar_05:    fltarr(5)-999.0,$
      k_ugriz_kcorrect_05:       fltarr(5)-999.0}
    kcorr = replicate(kcorr,ngal)

    kcorr.k_zobj = in_redshift
    kcorr.k_maggies = in_maggies
    kcorr.k_ivarmaggies = in_ivarmaggies
    
; compute k-corrections    
    splog, 'Computing ugriz K-corrections'
    ugriz_kcorrect_01 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,sdss_filterlist(),band_shift=0.1,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=ugriz_absmag_01,$
      ivarabsmag=ugriz_absmag_ivar_01,$;clineflux=cflux,uvflux=uvflux,$
      /silent,vname=vname,h100=h100) ; AB, band_shift=0.1

    kcorr.k_bestmaggies = bestmaggies
    kcorr.k_mass = alog10(mass) ; h=0.7, Chabrier
    kcorr.k_coeffs = coeffs
    kcorr.k_chi2 = chi2
    kcorr.k_ugriz_absmag_01         = ugriz_absmag_01
    kcorr.k_ugriz_absmag_ivar_01    = ugriz_absmag_ivar_01
    kcorr.k_ugriz_kcorrect_01       = ugriz_kcorrect_01

return, kcorr
end

pro ages_apply_zptoffsets, zpt, maggies, ivarmaggies, filters=filters, $
  absolute=absolute
; jm11mar25ucsd - apply the derived zeropoint offsets; maggies and
; ivarmaggies modified on output

    if (n_elements(filters) eq 0) then filters = strtrim(zpt[0].filterlist,2)
    dim = size(maggies,/dim)
    nfilt = dim[0]
    ngal = dim[1]

    if (nfilt ne n_elements(filters)) then message, 'Filter mismatch!'
    match, filters, strtrim(zpt[0].filterlist,2), m1, m2
    zptoffset = fltarr(nfilt)

    if keyword_set(absolute) then $
      zptoffset[m1] = total(zpt.zptoffset[m2],2) else $
        zptoffset[m1] = total(zpt.relative_zptoffset[m2],2)

    bigzptoffset = rebin(reform(zptoffset,nfilt,1),nfilt,ngal)
    factor = 10^(-0.4*bigzptoffset)
    maggies = maggies*factor
    ivarmaggies = ivarmaggies/factor^2

return
end

function init_zptout, filters, niter=niter

    nfilt = n_elements(filters)
    zptout = {filterlist: strarr(nfilt), iter: 0, converged: 0, $
      chi2cut: 0.0, nobj: lonarr(nfilt), used: intarr(nfilt)+1, $
      zptoffset: fltarr(nfilt), relative_zptoffset: fltarr(nfilt), $
      zptoffset_cumu: fltarr(nfilt), relative_zptoffset_cumu: fltarr(nfilt)}
    zptout = replicate(zptout,niter)

    zptout.filterlist = filters
    zptout.iter = lindgen(niter)

return, zptout
end

pro ages_derive_zptoffsets, clobber=clobber, debug=debug
; jm11mar25ucsd - iteratively derive zeropoint offsets

    zptpath = ages_path(/mycatalogs)+'zptoffsets/'
    niter = 20 ; 10
    minphot = 8 ; 3
    minerr = 0.03

    outfile = zptpath+'ages_zptoffsets.fits'
    logfile = repstr(outfile,'.fits','.log')
    if file_test(outfile+'.gz') and (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+outfile+' exists; use /CLOBBER!'
       return
    endif
    outfile_data = repstr(outfile,'.fits','_data.fits')

; read the sample and apply some cuts
    ages = read_ages(/phot)
    bootes_to_maggies, ages, maggies, ivarmaggies, filterlist=filters, $
      /nozpoffset, /nominerror
    maggies = maggies[0:7,*]
    ivarmaggies = ivarmaggies[0:7,*]
    filters = filters[0:7]
    select = (where(filters eq 'ndwfs_I.par'))[0]
    nfilt = n_elements(filters)

    keep = where((ages.z gt 0.05) and (ages.z lt 0.7) and $
      (total(maggies gt 0,1) ge minphot) and $
      (ages.i_tot gt 16.5) and (ages.i_tot lt 19.5),nobj)
    splog, 'Number of objects = '+strtrim(nobj,2)

    maggies = maggies[*,keep]
    ivarmaggies = ivarmaggies[*,keep]
    zobj = ages[keep].z

; apply a minimum photometric error across all the bands       
    mag = maggies2mag(maggies,ivar=ivarmaggies,magerr=magerr)
    good = where(magerr gt 0.0,ngood)
    if (ngood ne 0L) then begin
       magerr[good] = sqrt(magerr[good]^2 + minerr^2)
       junk = mag2maggies(mag,magerr=magerr,ivarmaggies=ivarmaggies)
    endif
    
    zptout = init_zptout(filters,niter=niter)
    mask = intarr(nfilt,nobj)+1 ; all are good

; iterate       
    splog, file=logfile
    t0 = systime(1)
    for iter = 0, niter-1 do begin

; apply the zero-point offsets and then fit
       good = where(total(mask,1) gt 0,ngood)
       in_zobj = zobj[good]
       in_maggies = maggies[*,good]
       in_ivarmaggies = ivarmaggies[*,good]

       isused = rebin(reform(zptout[0].used,nfilt,1),nfilt,ngood)
       zptout[iter].nobj = total(in_maggies gt 0,2)
       
       ages_apply_zptoffsets, zptout, in_maggies, in_ivarmaggies
       
       result = ages_zptoffsets_do_kcorrect(in_zobj,in_maggies,$
         isused*in_ivarmaggies,filterlist=filters)
       if (iter eq 0) then mwrfits, result, outfile_data, /create

; derive the zeropoint tweak for each band using the best model fits
       zptout[iter].chi2cut = weighted_quantile(result.k_chi2,quant=0.95)
       mask[*,good] = mask[*,good] and rebin(reform(result.k_chi2 lt $
         zptout[iter].chi2cut,1,ngood),nfilt,ngood)

       norm = total(mask[*,good]*in_ivarmaggies*result.k_bestmaggies^2,2,/double)
       dscale = total(mask[*,good]*in_ivarmaggies*in_maggies*result.k_bestmaggies,2,/double)/$
         (norm+(norm eq 0))*(norm ne 0)
       chi2 = total(mask[*,good]*in_ivarmaggies*(in_maggies-rebin(reform(dscale,nfilt,1),nfilt,nobj)*$
         result.k_bestmaggies)^2,/double,2)/total(mask[*,good],2)
       dzpt = +2.5*alog10(dscale)

       small = where(abs(dzpt) le 0.005,nsmall)
       if (nsmall ne 0) then dzpt[small] = 0.0
       zptout[iter].zptoffset = dzpt
       
; adjust the zeropoint corrections relative to the selection band 
       zptout[iter].relative_zptoffset = zptout[iter].zptoffset - $
         zptout[iter].zptoffset[select[0]]

; cumulative correction
       zptout[iter].zptoffset_cumu = total(zptout.zptoffset,2)
       zptout[iter].relative_zptoffset_cumu = total(zptout.relative_zptoffset,2)
       
; write out info       
       if (iter eq 0) then begin
          splog, format='("Time for the first iteration = '+$
            '",G0," minutes")', (systime(1)-t0)/60.0
          print
       endif
       
       splog, '##################################################'
       splog, 'Iteration '+string(iter,format='(I2.2)')
       splog, 'Chi^2 cut: '+strtrim(string(zptout[iter].chi2cut,format='(F12.2)'),2)
       for jj = 0, nfilt-1 do splog, filters[jj], $
         zptout[iter].zptoffset[jj], zptout[iter].relative_zptoffset[jj], $
         zptout[iter].zptoffset_cumu[jj], zptout[iter].relative_zptoffset_cumu[jj], $
         zptout[iter].used[jj], zptout[iter].nobj[jj], $
         format='(A30,F12.4,F12.4,F12.4,F12.4,2x,I0,2x,I0)'

       if keyword_set(debug) then begin
          for jj = 0, nfilt-1 do begin
             omag =-2.5*alog10(in_maggies[jj,*])
             bmag = -2.5*alog10(result.k_bestmaggies[jj])
             plot, omag, bmag-omag, ps=3, title=filters[jj], yr=[-0.5,0.5]
             oplot, !x.crange, dzpt[jj]*[1,1]
             cc = get_kbrd(1)
          endfor
       endif

       if (total(abs(zptout[iter].zptoffset) gt 0.0) eq 0.0) then begin
          zptout[iter].converged = 1
          splog, 'Zeropoints have converged'
          zptout = zptout[0:iter]

; fit one last time
          in_zobj = zobj[good]
          in_maggies = maggies[*,good]
          in_ivarmaggies = ivarmaggies[*,good]
          ages_apply_zptoffsets, zptout, in_maggies, in_ivarmaggies
          result = ages_zptoffsets_do_kcorrect(in_zobj,in_maggies,$
            isused*in_ivarmaggies,filterlist=filters)
          mwrfits, result, outfile_data
          spawn, 'gzip -f '+outfile_data, /sh
          break                 ; leave
       endif 
       
; write out
       mwrfits, zptout, outfile, /create
       spawn, 'gzip -f '+outfile, /sh
       splog, /close
    endfor
       
return
end 
