function ediscs_filterlist
; jm09aug17ucsd - written
    filterlist = [['FORS_B','FORS_V','FORS2_R','FORS_I']+$
      '_ccd',['SOFI_J','SOFI_Ks']]+'_atm.par'
return, filterlist
end
