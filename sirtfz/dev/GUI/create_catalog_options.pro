pro options_event, event

    common sirtf_simulations

    event_name = tag_names(event,/structure_name)

stop

    case event_name of
       'WIDGET_LIST': begin ; list option

          widget_control, event.id, get_uvalue = info

;         case 

          print, select
          stop

       end
       'WIDGET_BUTTON': begin   ; button option
          widget_control, event.id, get_uvalue=option
          case option of
             'help': help = dialog_message(['Help!'],/info,dialog_parent=event.top)
             'done': widget_control, event.top, /destroy
          endcase
       end
    endcase
    
return
end

pro create_catalog_options, filters, dnds, dndz
; jm01jan28uofa

    common sirtf_simulations

    options = ['Number Counts','Redshift Distribution','Simulated Image']
    
    mainbase = widget_base(title='Create Catalog Options',/row,$
                           /align_bottom,group_leader=sirtf.widget_ids.base_id,$
                           /base_align_center)

    list1 = widget_list(mainbase,event_pro='options_event',value=options,ysize=3.5)
    list2 = widget_list(mainbase,uvalue=filters,event_pro='options_event', $
                        value=filters,ysize=3.5)

    buttonbase = widget_base(mainbase,/column)
    help_button = widget_button(buttonbase,value='Help',uvalue='help')
    done_button = widget_button(buttonbase,value='Done',uvalue='done')

    info = {options: options, filters: filters, index: 0L}

    widget_control, mainbase, set_uvalue=info, /realize
    xmanager, 'create_catalog_options', mainbase, event_handler = 'options_event'
    
return
end
