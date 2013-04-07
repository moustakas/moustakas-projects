pro cww
; jm01apr9uofa
; convert the SED text files from HyperZ to SIRTFz format.  HyperZ
; format has the SEDs in angstroms and in f_lambda; convert to microns
; and f_nu (arbitrary units)
    
    path = '/home/ioannis/sirtf/sirtfz/seds/cww/'

    fnames = ['CWW_E_ext.sed','CWW_Im_ext.sed','CWW_Sbc_ext.sed','CWW_Scd_ext.sed']
    outnames = ['CWW_E.sed','CWW_Im.sed','CWW_Sbc.sed','CWW_Scd.sed']

    for j = 0L, n_elements(fnames)-1L do begin

       readfast, path+fnames[j], sed, nlines=nlines

       lum = sed[0,*] * sed[0,*] * sed[1,*] / 2.99793D18 ; f_nu [arbitrary units]
       lambda = sed[0,*]*1D-4                            ; wavelength [mu]

       openw, lun, path+outnames[j], /get_lun
       for k = 0L, nlines-1L do printf, lun, lambda[k], lum[k], format='(2x,E9.3,2x,E10.4)'
       free_lun, lun

    endfor
    
return
end
