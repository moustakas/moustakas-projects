pro fill_sfrm_results, indx, ngal, mf_fit, mf_data, $
  fit, binmass, phi, phierr, fullbin, number
; jm10feb16ucsd - see SFR_MF
    mf_fit[indx] = im_struct_assign(fit,mf_fit[indx],/nozero)
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
