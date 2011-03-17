pro sirtfz_shutdown, sedcube=sedcube, bandcube=bandcube, all=all
; destroy the widgets, kill the windows, and free the pointers.
; allows you to destroy pointers in structures piecemeal
    
    common sirtf_simulations, sirtf

    if keyword_set(sirtf) then begin ; check for definition
    
       if keyword_set(all) then begin

          sedcube = 1L
          bandcube = 1L

          if (xregistered('sirtfz')) then widget_control, sirtf.widgets.base_id, /destroy
          if ptr_valid(sirtf.redshift.zarray) then ptr_free, sirtf.redshift.zarray
          if ptr_valid(sirtf.filters) then ptr_free, sirtf.filters
          if obj_valid(progress) then obj_destroy, progress
          
       endif
       
       if keyword_set(sedcube) then begin
          
          nseds = size(sirtf.sedcube,/n_elements)
          if long(total(ptr_valid(sirtf.sedcube.lambda)))/nseds then ptr_free, sirtf.sedcube[*].lambda
          if long(total(ptr_valid(sirtf.sedcube.mlum_nu)))/nseds then ptr_free, sirtf.sedcube[*].mlum_nu
          if long(total(ptr_valid(sirtf.sedcube.mlum_mu)))/nseds then ptr_free, sirtf.sedcube[*].mlum_mu
          
       endif
       
       if keyword_set(bandcube) then begin
          
          nbands = size(sirtf.bandcube,/n_elements)
          if long(total(ptr_valid(sirtf.bandcube.rband)))/nbands then ptr_free, sirtf.bandcube[*].rband
          if long(total(ptr_valid(sirtf.bandcube.wband)))/nbands then ptr_free, sirtf.bandcube[*].wband
          
       endif

    endif
       
return
end
