pro cww_sb
; jm01sep19uofa
; convert the SED text files from BPZ to SIRTFz format.  BPZ format
; has the SEDs in angstroms and in f_lambda; convert to microns and
; f_nu (arbitrary units)
    
    path = '/home/ioannis/sirtf/sirtfz/seds/cww/'

    fnames = ['SB2_kin.sed','SB3_kin.sed']
    outnames = ['SB2.sed','SB3.sed']

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
