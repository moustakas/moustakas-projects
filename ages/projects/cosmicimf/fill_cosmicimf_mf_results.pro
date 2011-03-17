pro fill_cosmicimf_mf_results, indx, ngal, mf_fit, mf_data, $
  fit, binmass, phi, phierr, fullbin, number, nolss=nolss
; jm10mar23ucsd - see COSMICIMF_MF
    mf_fit[indx] = im_struct_assign(fit,mf_fit[indx],/nozero)
; convert the fractional LSS error to Phi* units
    if (keyword_set(nolss) eq 0) then begin
       lss = ages_sigmalss()
       mf_fit[indx].phistar_lss_err = mf_fit[indx].phistar*lss[indx]
       mf_fit[indx].rho_lss_err = mf_fit[indx].rho*lss[indx]
    endif
; fill the data structure    
    nbins = n_elements(phi)
    mf_data[indx].ngal = ngal
    mf_data[indx].nbins = nbins
    mf_data[indx].fullbin[0:nbins-1] = fullbin
    mf_data[indx].number[0:nbins-1] = number
    mf_data[indx].mass[0:nbins-1] = binmass
    mf_data[indx].phi[0:nbins-1] = phi
    mf_data[indx].phierr[0:nbins-1] = phierr
return
end
