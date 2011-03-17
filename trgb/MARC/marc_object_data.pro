; jm00june27ucb

pro marc_object_data

	light = 2.99E5	; km/s

	hst = sread('/deep1/ioannis/trgb/hst_object.dat')
	keck = sread('/deep1/ioannis/trgb/keck_object.dat')
        
        openw, lun1, '/deep1/ioannis/trgb/hst_object.txt', /get_lun
        openw, lun2, '/deep1/ioannis/trgb/keck_object.txt', /get_lun
        openw, lun3, '/deep1/ioannis/trgb/data_object.txt', /get_lun

        for k = 0L, n_elements(hst)-1L do begin
            printf, lun1, $
              hst[k].object, $
              hst[k].longitude, hst[k].latitude, $
;             hst[k].v_helio, $
            format = '(A10,1x,3F12.5)'
;             format = '(3F12.5)'
            printf, lun3, $
              hst[k].object, $
              hst[k].longitude, hst[k].latitude, $
;             hst[k].v_helio, $
            format = '(A10,1x,3F12.5)'
 ;            format = '(3F12.5)'
        endfor

        free_lun, lun1

        for k = 0L, n_elements(keck)-1L do begin
            printf, lun2, $
              keck[k].object, $
              keck[k].longitude, keck[k].latitude, $
;             keck[k].v_helio, $
            format = '(A10,1x,3F12.5)'
;             format = '(3F12.5)'
            printf, lun3, $
              keck[k].object, $
              keck[k].longitude, keck[k].latitude, $
;             keck[k].v_helio, $
            format = '(A10,1x,3F12.5)'
;             format = '(3F12.5)'
        endfor

        free_lun, lun2
        free_lun, lun3        

return
end
