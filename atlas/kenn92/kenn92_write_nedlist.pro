pro kenn92_write_nedlist, write=write
; jm03may27uofa
; jm05aug02uofa - updated    

    datapath = kenn92_path(/analysis)
    nedpath = datapath+'ned/'
    rob = mrdfits(datapath+'rob2.fits',0,h,/silent)

    ngalaxy = 55L
    galaxy = strarr(ngalaxy)
    galaxy[0] = 'NGC3379'

    for i = 1L, ngalaxy-1L do galaxy[i] = sxpar(h,'APID'+string(i+1,format='(I0)'))

; modify KENN92_NED_INPUT.TXT to GALAXY_LIST.TXT
    
    if keyword_set(write) then begin
    
       splog, 'Writing '+nedpath+'kenn92_ned_input.txt.'
       openw, lun, nedpath+'kenn92_ned_input.txt', /get_lun
       niceprintf, lun, galaxy
       free_lun, lun

    endif else niceprint, galaxy

return
end    
