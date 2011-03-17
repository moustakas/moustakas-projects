function make_zarray, zmin, zmax, dz

    if zmin eq float(0) then k = zmin+dz else k = zmin
    zarray = (findgen((zmax-zmin)/dz)*dz)+k
    
return, zarray
end

pro redshift_event, event

    common sirtf_simulations

    widget_control, event.id, get_uvalue = event_name

    case event_name of
       'zmin': begin
          if (*event.value gt sirtf.redshift.zmax) then $
            temp = dialog_message('The minimum redshift cannot exceed the maximum redshift!',$
                                  /error,dialog_parent=event.top) else begin
             if ptr_valid(sirtf.redshift.zarray) then ptr_free, sirtf.redshift.zarray    ; free the existing pointer
             new_array = make_zarray(*event.value,sirtf.redshift.zmax,sirtf.redshift.dz) ; make the new array
             sirtf.redshift.zmin = *event.value
             sirtf.redshift.zarray = ptr_new(new_array)
          endelse
       end
       'zmax': begin
          if (*event.value lt sirtf.redshift.zmin) then $
            temp = dialog_message('The maximum redshift cannot be less than the minimum redshift!',$
                                  /error,dialog_parent=event.top) else begin
             if ptr_valid(sirtf.redshift.zarray) then ptr_free, sirtf.redshift.zarray     
             new_array = make_zarray(sirtf.redshift.zmin,*event.value,sirtf.redshift.dz)
             sirtf.redshift.zmax = *event.value
             sirtf.redshift.zarray = ptr_new(new_array)
          endelse
       end
       'deltaz': begin
          zdiff = sirtf.redshift.zmax-sirtf.redshift.zmin
          if (*event.value gt zdiff) then $
            temp = dialog_message('This is an invalid value for the redshift binsize.',$
                                  /error,dialog_parent=event.top) else begin
             if ptr_valid(sirtf.redshift.zarray) then ptr_free, sirtf.redshift.zarray     
             new_array = make_zarray(sirtf.redshift.zmin,sirtf.redshift.zmax,*event.value)
             sirtf.redshift.dz = *event.value
             sirtf.redshift.zarray = ptr_new(new_array)
          endelse
       end
       'help': help = dialog_message(['Create a redshift array by selecting a starting',$
                                      'redshift (zmin) and an ending redshift (zmax).',$
                                      'The binsize is given by delta z.  You can use',$
                                      'the TAB button to toggle between the fields.'],$
                                     /info,dialog_parent=event.top)
       'done': widget_control, event.top, /destroy
    endcase

    sirtf.redshift.nz = n_elements(*sirtf.redshift.zarray)

return
end

pro redshift_control
; jm01jan12uofa

    common sirtf_simulations

    zmin = sirtf.redshift.zmin
    zmax = sirtf.redshift.zmax
    dz = sirtf.redshift.dz

    base = widget_base(title='Model Redshift Grid',/column,/align_bottom,$
                       group_leader=sirtf.widget_ids.base_id,/base_align_center)
    field1 = fsc_inputfield(base,title='Minimum z',/doublevalue,decimal=6,/positive,$
                            /cr_only,value=zmin,uvalue='zmin',event_pro='redshift_event',$
                            labelsize=60,xsize=7,/frame)
    field2 = fsc_inputfield(base,title='Maximum z',/doublevalue,decimal=6,/positive,$
                            /cr_only,value=zmax,uvalue='zmax',event_pro='redshift_event',$
                            labelsize=60,xsize=7,/frame)
    field3 = fsc_inputfield(base,title='Resolution',/doublevalue,decimal=6,/positive,$
                            /cr_only,value=dz,uvalue='deltaz',event_pro='redshift_event',$
                            labelsize=60,xsize=7,/frame)
    help_button = widget_button(base,value='Help',uvalue='help')
    done_button = widget_button(base,value='Done',uvalue='done')

    widget_control, base, /realize;, set_uvalue = {zmin: zmin, zmax: zmax, dz: dz}, /no_copy
    xmanager, 'redshift_control', base, event_handler='redshift_event', /no_block

return
end
