pro fill_lf24_results, indx, ngal, lf24_fit, lf24_data, $
  fit, binlum, phi, phierr, fullbin, number
; jm10mar25ucsd - see COSMICILF24_LF24
    lf24_fit[indx] = im_struct_assign(fit,lf24_fit[indx],/nozero)
    nbins = n_elements(phi)
    lf24_data[indx].ngal = ngal
    lf24_data[indx].nbins = nbins
    lf24_data[indx].fullbin[0:nbins-1] = fullbin
    lf24_data[indx].number[0:nbins-1] = number
    lf24_data[indx].l24[0:nbins-1] = binlum
    lf24_data[indx].phi[0:nbins-1] = phi
    lf24_data[indx].phierr[0:nbins-1] = phierr
return
end
