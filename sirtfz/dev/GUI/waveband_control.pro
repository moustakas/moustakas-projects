pro waveband_event, event

    common sirtf_simulations

    event_name = tag_names(event,/structure_name)

    case event_name of
       'WIDGET_LIST': begin ; list option

          widget_control, event.id, get_uvalue = filters
          indices = widget_info(event.id,/list_select)  ; which filters have been selected?
          nindices = n_elements(indices)                ; how many were selected?

          if nindices eq 1L then begin  ; one selected
             if (indices eq -1L) then $ ; none selected
               temp = dialog_message(['You must select at least one filter!'],$
                                     /error,dialog_parent=event.top) else begin
                if ptr_valid(sirtf.filters) then ptr_free, sirtf.filters ; free the previous filters
                sirtf.filters = ptr_new(filters[indices])                ; assign the new filters
             endelse
          endif else begin              ; multiple selected
             if ptr_valid(sirtf.filters) then ptr_free, sirtf.filters
             sirtf.filters = ptr_new(filters[indices])
          endelse

       end
       'WIDGET_BUTTON': begin   ; button option
          widget_control, event.id, get_uvalue=option
          case option of
             'help': help = dialog_message(['Select one or more filters with the',$
                                            'mouse button, using the SHIFT and CNTL',$
                                            'keys to make multiple selections.'],$
                                           /info,dialog_parent=event.top)
             'done': widget_control, event.top, /destroy
          endcase
       end
    endcase
    
return
end

pro waveband_control
; jm01jan8uofa
; droplist widget that allows the user to select a particular bandpass

    common sirtf_simulations

    filters = sirtf.bandcube.bandnames ; available filters
    
    base = widget_base(title='Filter Selection',/column,$
                       /align_bottom,group_leader=sirtf.widget_ids.base_id,$
                       /base_align_center)
    list = widget_list(base,uvalue=filters,event_pro='waveband_event', $
                       /multiple,value=filters,ysize=8.5)
    help_button = widget_button(base,value='Help',uvalue='help')
    done_button = widget_button(base,value='Done',uvalue='done')

    widget_control, base, /realize
    xmanager, 'waveband_control', base, event_handler = 'waveband_event'
    
return
end

