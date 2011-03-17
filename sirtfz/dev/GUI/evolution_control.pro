pro evolution_event, event

    common sirtf_simulations
    
    widget_control, event.id, get_uvalue = event_name

    case event_name of
       'alpha': sirtf.evolution.alpha = *event.value
       'beta': sirtf.evolution.beta = *event.value
       'help': help = dialog_message(['These parameters allow you to control',$
                                      'the amount of luminosity and density',$
                                      'evolution to apply to the luminosity',$
                                      'function.'],/info,dialog_parent=event.top)
       'done': widget_control, event.top, /destroy
    endcase

return
end

pro evolution_control

    common sirtf_simulations

    base = widget_base(title='Evolution Parameters',/column,$;/align_bottom,$
                       group_leader=sirtf.widget_ids.base_id,/base_align_center)
    evol1 = fsc_inputfield(base,title='Luminosity Evolution',/doublevalue,event_pro='evolution_event', $
                           uvalue='alpha',value=sirtf.evolution.alpha,/cr_only,/positive,decimal=2,$
                           /frame,labelsize=120,xsize=5)
    evol2 = fsc_inputfield(base,title='Density Evolution',/doublevalue,event_pro='evolution_event', $
                           uvalue='beta',value=sirtf.evolution.beta,/cr_only,/positive,decimal=2,$
                           /frame,labelsize=120,xsize=5)
    help_button = widget_button(base,value='Help',uvalue='help')
    done_button = widget_button(base,value='Done',uvalue='done')

    widget_control, base, /realize
    xmanager, 'evolution_control', base, event_handler='evolution_event', /no_block

return
end


