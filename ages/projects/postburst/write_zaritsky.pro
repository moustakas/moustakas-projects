pro write_zaritsky
; jm04sep14uofa

    ages = read_ages(/silent)
    tags = ['ZOSU_PASSFIELD','CATALOG_NUMBER','OII_3727','OII_3727_EW','BABS_H_DELTA_EW',$
      'BABS_H_GAMMA_EW','BABS_H_BETA_EW','BABS_H_ALPHA_EW','CONTINUUM_CHI2']
    newtags = ['PLATE','CATALOG_NUMBER','F_OII','EW_OII','EW_HD','EW_HG','EW_HB','EW_HA','CHI2']
    line = im_struct_trimtags(ages,select=tags,newtags=newtags,format=['I3','I5',replicate('E12.5',2),$
      replicate('F12.5',11)])

    outpath = ages_path(/zaritsky)
    openw, lun, outpath+'ages_table.dat', /get_lun
    struct_print, line, lun=lun
    free_lun, lun

return
end
    
