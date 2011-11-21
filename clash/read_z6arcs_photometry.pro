function read_z6arcs_photometry, adi=adi
; jm11nov18ucsd - read the *final* photometry of the arcs

    datapath = clash_path(/macs0329_z6arcs)

    filt = clash_filterlist(/useirac,weff=weff)
    nfilt = n_elements(filt)
    
    irac = rsex(datapath+'macs0329_z6arcs_irac.cat')
    adi = rsex(datapath+'macs0329_z6arcs.cat')
    struct_print, adi
    nobj = n_elements(adi)

    hstindx = lindgen(16)
    iracindx = lindgen(2)+16
    
    cat = replicate({weff: weff, maggies: fltarr(nfilt), ivarmaggies: fltarr(nfilt), $
      limit: intarr(nfilt), mag: fltarr(nfilt), magerr_up: fltarr(nfilt), $
      magerr_lo: fltarr(nfilt), magerr: fltarr(nfilt)},nobj)
    file = datapath+'z6arc1.'+['1','2','3','4']+'HSTphotometry.cat'
    for jj = 0, nobj-1 do begin
       phot = rsex(file[jj])
       cat[jj].weff = phot.lam
       
       good = where(finite(phot.fluxerr))
       cat[jj].maggies[hstindx[good]] = phot[good].flux*10D^(-0.4*phot[good].zpext)
       cat[jj].ivarmaggies[hstindx[good]] = 1.0/(phot[good].fluxerr*10D^(-0.4*phot[good].zpext))^2

; for the plot
       sig = where(phot.signif gt 2.0,comp=insig)
       cat[jj].mag[sig] = phot[sig].mag
       cat[jj].magerr_up[sig] = phot[sig].magerrhi
       cat[jj].magerr_lo[sig] = phot[sig].magerrlo
       cat[jj].magerr[sig] = phot[sig].magerr
       
       cat[jj].limit[insig] = 1
       cat[jj].mag[insig] = phot[insig].mag2sig
    endfor

; add IRAC
    cat.maggies[iracindx[0]] = 10D^(-0.4*irac.ch1_ab)
    cat.ivarmaggies[iracindx[0]] = 1.0/(0.4*alog(10)*(cat.maggies[iracindx[0]]*irac.ch1_aberr))^2
    cat.mag[iracindx[0]] = irac.ch1_ab
    cat.magerr_up[iracindx[0]] = irac.ch1_aberr
    cat.magerr_lo[iracindx[0]] = irac.ch1_aberr
    
    cat.maggies[iracindx[1]] = 10D^(-0.4*irac.ch2_ab)
    cat.ivarmaggies[iracindx[1]] = 1.0/(0.4*alog(10)*(cat.maggies[iracindx[1]]*irac.ch2_aberr))^2
    cat.mag[iracindx[1]] = irac.ch2_ab
    cat.magerr_up[iracindx[1]] = irac.ch2_aberr
    cat.magerr_lo[iracindx[1]] = irac.ch2_aberr

; 2-sigma limit for ch2 on arc 1.4
    cat[3].ivarmaggies[iracindx[1]] = 1.0/(10D^(-0.4*irac[3].ch2_ab))^2.0
    cat[3].maggies[iracindx[1]] = 0.0
    cat[3].limit[iracindx[1]] = 1
    cat[3].mag[iracindx[1]] = irac[3].ch2_ab
    cat[3].magerr_up[iracindx[1]] = 0.0
    cat[3].magerr_lo[iracindx[1]] = 0.0
    
return, cat
end
