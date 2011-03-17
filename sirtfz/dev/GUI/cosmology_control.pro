pro cosmology_event, event

    common sirtf_simulations
    common cosmology
    
    widget_control, event.id, get_uvalue = event_name

    case event_name of
       'omega_0': cosmo_params.omega_0 = *event.value
       'omega_lambda': cosmo_params.omega_lambda = *event.value
       'h_100': cosmo_params.h_100 = *event.value
       'help': help = dialog_message(['Select the cosmological parameters.'],/info,$
                                     dialog_parent=event.top)
       'done': widget_control, event.top, /destroy
    endcase

return
end

pro cosmology_control
; widget sequence that allows the user to change the cosmological parameters
        
    common sirtf_simulations
    common cosmology
    
    base = widget_base(title='Cosmology Parameters',/column,$
                       group_leader=sirtf.widget_ids.base_id,/base_align_center)
    cosmo1 = fsc_inputfield(base,title='Matter Density',/doublevalue,event_pro='cosmology_event', $
                            uvalue='omega_0',value=cosmo_params.omega_0,/cr_only,$
                            /positive,decimal=2,/frame,labelsize=90,xsize=5)
    cosmo2 = fsc_inputfield(base,title='Vacuum Density',/doublevalue,event_pro='cosmology_event', $
                            uvalue='omega_lambda',value=cosmo_params.omega_lambda,$
                            /cr_only,/positive,decimal=2,/frame,labelsize=90,xsize=5)
    cosmo3 = fsc_inputfield(base,title='Hubble Constant',/doublevalue,event_pro='cosmology_event', $
                            uvalue='h_100',value=cosmo_params.h_100,/cr_only,/positive,$
                            decimal=2,labelsize=90,/frame,xsize=5)
    help_button = widget_button(base,value='Help',uvalue='help')
    done_button = widget_button(base,value='Done',uvalue='done')
    
    widget_control, base, /realize;, set_uvalue = {omega_0: omega_0, omega_lambda: omega_lambda, h_100: h_100}
    xmanager, 'cosmology_control', base, event_handler='cosmology_event', /no_block

return
end


