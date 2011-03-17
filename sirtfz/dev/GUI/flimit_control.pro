pro flimit_event, event

    common sirtf_simulations
    
    widget_control, event.id, get_uvalue = event_name

    case event_name of
       'flux': begin ; need to update the magnitude field
          event.objref -> getproperty, name=filter
          indx = filter_match(filter)
          sirtf.bandcube[indx].flimit = *event.value
       end          
       'mag' : begin ; need to update the flux field
          event.objref -> getproperty, name=filter
          indx = filter_match(filter)
          flux = 10D^(-0.4D*(*event.value-16.40D)) ; convert to flux
          sirtf.bandcube[indx].flimit = flux
       end
       'help': help = dialog_message(['Change the default limiting flux', $
                                      'for a given filter (units are mJy).'],$
                                     /info,dialog_parent=event.top)
       'done': widget_control, event.top, /destroy


       else: begin          
          match = filter_match(event_name)
          if match[0] eq -1L then $
            err = dialog_message(['No known filters match the selected filter.'],$
                                 /error,dialog_parent=event.top) else $
            sirtf.bandcube[match].flimit = *event.value
       end
       
       else: temp = dialog_message('Not a valid menu option!',/error, $
                                      dialog_parent=event.top)
    endcase

return
end

pro flimit_control

    common sirtf_simulations

    filters = sirtf.bandcube.bandnames ; available filters

    flimit = string(sirtf.bandcube.flimit,format='(G12.4)')
    fmag = -2.5D*alog10(float(flimit))+16.40D ; mJy

    print, sirtf.bandcube.flimit
    
    biggerbase = widget_base(title='Limiting Flux',/column,group_leader=sirtf.widget_ids.base_id)
    bigbase = widget_base(biggerbase,/row,/base_align_center)
    col1 = widget_base(bigbase,/column)
    col2 = widget_base(bigbase,/column)
    col3 = widget_base(bigbase,/column)

    label = widget_label(col1,/align_center,value='Filter',xsize=3.0,ysize=1.,units=2)
    label = widget_label(col2,/align_center,value='Flux (mJy)',xsize=3.,units=2,ysize=1.)
    label = widget_label(col3,/align_center,value='AB mag',xsize=2.8,units=2,ysize=1.)

    for k = 0L, n_elements(filters)-1L do $
      label = widget_label(col1,value=filters[k],xsize=3.,/frame,/align_center,$
                           ysize=1.255,units=2)

    for k = 0L, n_elements(filters)-1L do begin
       field = fsc_inputfield(col2,title='',value=flimit[k],/cr_only,xsize=5.5,/frame,$
                              labelsize=0,uvalue='flux',event_pro='flimit_event',$
                              name=filters[k])
       field->resize, 85
    endfor

    for k = 0L, n_elements(filters)-1L do begin
       field = fsc_inputfield(col3,title='',value=fmag[k],/cr_only,xsize=4.5,/frame,$
                              labelsize=0,/float,decimal=3,uvalue='mag',event_pro='flimit_event',$
                              name=filters[k])
       field->resize, 60
    endfor
    
    buttonbase = widget_base(biggerbase,/row,/base_align_right)
    help_button = widget_button(buttonbase,value='Help',uvalue='help',/align_right)
    done_button = widget_button(buttonbase,value='Done',uvalue='done',/align_right)

    widget_control, biggerbase, /realize
    xmanager, 'flimit_control', biggerbase, event_handler='flimit_event', /no_block

return
end


