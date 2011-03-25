pro mz_remove_old_zptoffsets, maggies, ivarmaggies, old_zptoffset=old_zptoffset
; jm10oct19ucsd - remove any zeropoint offsets applied to the
; optical/near-IR bands in MZ_TO_MAGGIES, if any

    dim = size(maggies,/dim)
    nfilt = dim[0] & ngal = dim[1]
    junk = mz_filterlist(zpoffset=old_zptoffset)
    
    bigold_zptoffset = rebin(reform(old_zptoffset,nfilt,1),nfilt,ngal)
    factor = 10^(+0.4*bigold_zptoffset)
    maggies = maggies*factor
    ivarmaggies = ivarmaggies/factor^2
    
return
end

pro mz_apply_zptoffsets, zptout, maggies, ivarmaggies
; jm10oct19ucsd - apply the derived zeropoint offsets; maggies and
; ivarmaggies modified on output 
;
    dim = size(maggies,/dim)
    nfilt = dim[0] & ngal = dim[1]
    zptoffset = total(zptout.zptoffset,2) ; cumulative correction
    bigzptoffset = rebin(reform(zptoffset,nfilt,1),nfilt,ngal)
    factor = 10^(-0.4*bigzptoffset)
    maggies = maggies*factor
    ivarmaggies = ivarmaggies/factor^2
return
end

function init_zptout, nfilt, niter=niter
    zptout = {iter: 0, nobj: lonarr(nfilt), zptoffset: fltarr(nfilt), $
      zptoffset_orig: fltarr(nfilt)}
    zptout = replicate(zptout,niter)
    zptout.iter = lindgen(niter)
return, zptout
end

pro derive_mz_zptoffsets, clobber=clobber
; jm10oct19ucsd - iteratively fit the AGES galaxies with either
; isedfit or K-correct to derive zeropoint offsets, if any 

    zptpath = ages_path(/projects)+'mz/zptoffsets/'
    niter = 15

;   paramfile = isedpath+field[ii]+'_isedfit.par'
;   params = read_isedfit_paramfile(paramfile,sfhgrid='02')

    outfile = zptpath+'ages_zptoffsets.fits'
    logfile = repstr(outfile,'.fits','.log')
    if file_test(outfile+'.gz') and (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+outfile+' exists; use /CLOBBER!'
       return
    endif

    sample = read_mz_sample(/parent)
;   sample = sample[500:800]
    
    maggies = sample.k_maggies
    ivarmaggies = sample.k_ivarmaggies
    zobj = sample.z

    filters = mz_filterlist()
    nfilt = n_elements(filters)
    ifaint = mz_ifaint(select_filter=select_filter)
    select = where(filters eq select_filter)
    
    zptout = init_zptout(nfilt,niter=niter)

; remove old zeropoint offsets
    mz_remove_old_zptoffsets, maggies, ivarmaggies

;   splog, 'Testing!!'
;   test = mrdfits(outfile+'.gz',1)
;   mz_apply_zptoffsets, test, maggies, ivarmaggies

; iterate       
    splog, file=logfile
    t0 = systime(1)
    for iter = 0, niter-1 do begin
; apply the zero-point offsets and then fit
       in_maggies = maggies
       in_ivarmaggies = ivarmaggies
       mz_apply_zptoffsets, zptout, in_maggies, in_ivarmaggies

       result = mz_kcorrect(zobj,in_maggies,in_ivarmaggies,filterlist=filters,/just_ugriz)
       
       for jj = 0, nfilt-1 do begin
          good = where(in_maggies[jj,*] gt 0.0,nobj)
          zptout[iter].nobj[jj] = nobj
          if (nobj ne 0L) then begin
             dmag = reform(+2.5*alog10(in_maggies[jj,good]/result[good].k_bestmaggies[jj]))
             dmagmed = median(dmag)
             if (abs(dmagmed) lt 0.003) then dmagmed = 0.0
             zptout[iter].zptoffset_orig[jj] = dmagmed
          endif
       endfor

; adjust the zeropoint corrections so that they are relative to the
; selection band
       zptout[iter].zptoffset = zptout[iter].zptoffset_orig - $
         zptout[iter].zptoffset_orig[select[0]]
       
; write out info       
       if (iter eq 0) then begin
          splog, format='("Time for the first iteration = '+$
            '",G0," minutes")', (systime(1)-t0)/60.0
          print
       endif
       
       splog, '##################################################'
       splog, 'Iteration '+string(iter,format='(I2.2)')
       cumu_zptoffset_orig = total(zptout.zptoffset_orig,2) ; cumulative correction
       cumu_zptoffset = total(zptout.zptoffset,2)
       for jj = 0, nfilt-1 do splog, filters[jj], $
         zptout[iter].zptoffset_orig[jj], cumu_zptoffset_orig[jj], $
         zptout[iter].zptoffset[jj], cumu_zptoffset[jj], $
         format='(A25,F12.4,F12.4,F12.4,F12.4)'

       if (total(abs(zptout[iter].zptoffset) gt 0.0) eq 0.0) then begin
          splog, 'Zeropoints have converged'
          break
       endif
    endfor
       
; write out
    im_mwrfits, zptout, outfile, /clobber
    splog, /close

return
end 
