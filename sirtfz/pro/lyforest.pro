pro lyforest
; jm01aug8uofa
; generate a table of differential extinction corrections due to
; intervening intergalactic gas (Madau 1995, ApJ, 441, 18)

    path = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='lib')
    fname = 'lyforest.dat'
    
    fwave = (findgen(10/0.00001)+1.0)*0.00001 ; wavelength vector (micron) [1 Angstrom - 10 micron]
    nwave = size(wave,/n_elements)
    
    zarray = findgen(10.05/0.05)*0.05 ; redshift vector [0,10]
    nz = size(zarray,/n_elements)

; define some constants from the paper
    
    lambda = [1216.0,1026.0,973.0,950.0]*1E-4 ; Ly-alpha through Ly-delta (micron)
    
    aconst = [3.6E-3, $ ; Ly-alpha
              1.7E-3, $ ; Ly-beta
              1.2E-3, $ ; Ly-gamma 
              9.3E-4]   ; Ly-delta

    openw, lun, path+fname, /get_lun
    printf, lun, '# Lyman Forest Depletion '+strmid(systime(),4,20)
    printf, lun, '# z, D_A, D_B (Madau 1995, ApJ, 441, 18) '

    for i = 0L, nz-1L do begin  ; loop on redshift

       z = zarray[i]
       zemit = (1.0+z)

; constants
       
       dlya = 120.0E * 1E-4 * zemit ; micron
       dlyb = 95.0E * 1E-4 * zemit  ; micron

; calculate d_A and d_B

       get_element, fwave, 1050.0*1E-4*zemit, wmina
       get_element, fwave, 1170.0*1E-4*zemit, wmaxa

       get_element, fwave, 920.0*1E-4*zemit, wminb
       get_element, fwave, 1015.0*1E-4*zemit, wmaxb

; 1 - <D_A>

       d_a = (1.0/dlya) * int_tabulated(fwave[wmina:wmaxa],$
                                        exp(-aconst[0]*(fwave[wmina:wmaxa]/lambda[0])^(3.46)),/double)
; 1 - <D_B>
       
       d_b = (1.0/dlyb) * (int_tabulated(fwave[wminb:wmaxb],$
                                         exp(-(aconst[1]*(fwave[wminb:wmaxb]/lambda[1])^(3.46) + $
                                               aconst[2]*(fwave[wminb:wmaxb]/lambda[2])^(3.46) + $
                                               aconst[3]*(fwave[wminb:wmaxb]/lambda[3])^(3.46))),/double))

       d_a = d_a > 0.0
       d_b = d_b > 0.0

       printf, lun, z, d_a, d_b, format='(3F15.7)'
       
    endfor
    free_lun, lun

return
end
