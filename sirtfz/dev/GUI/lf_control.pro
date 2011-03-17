pro lf_event, event

    common sirtf_simulations

    event_name = tag_names(event,/structure_name)

    case event_name of
       'WIDGET_BUTTON': begin
          widget_control, event.id, get_uvalue=option
          case option of
             'help': help = dialog_message(['Help!  Put in the reference to Guiderdoni.'],$
                                           /info,dialog_parent=event.top)
             'done': widget_control, event.top, /destroy
          endcase
       end
    endcase

return
end

pro lf_control
; jm01jan13uofa
; luminosity function control:  allow the user to examine the input
; luminosity function

    common sirtf_simulations

    base = widget_base(title='Luminosity Function',/column,/align_bottom,$
                       group_leader=sirtf.widget_ids.base_id,/base_align_right)
    draw_window = widget_draw(base,uvalue='draw_window',scr_xsize=620L,scr_ysize=540L)
    buttonbase = widget_base(base,/row,group_leader=sirtf.widget_ids.base_id)
    help_button = widget_button(buttonbase,value='Help',uvalue='help',units=1,xsize=1.5)
    done_button = widget_button(buttonbase,value='Done',uvalue='done',units=1,xsize=1.5)

    widget_control, base, /realize
    xmanager, 'lf_control', base, event_handler='lf_event', /no_block

; plot the luminosity function

    lum_lf = sirtf.lf.lum_lf
    phi_lf = sirtf.lf.phi_lf

    plotfaves

    widget_control, draw_window, get_value=window_index
    wset, window_index
    plotsym, 0, 1, /fill
    plot, lum_lf, phi_lf, ps=8, xsty=3, ysty=3, xrange=[6.5,14], yrange=[-9,0], $
      title='IRAS 60 '+textoidl('\mu')+'m Luminosity Function', $
      xtitle='log '+textoidl('\nu')+'L'+textoidl('_{\nu}')+' (L'+textoidl('_{'+sunsymbol()+'}')+')', $
      ytitle='log '+textoidl('\Phi')+' ('+textoidl('\nu')+'L'+textoidl('_{\nu}')+') Mpc'+$
      textoidl('^{-3}')+' mag'+textoidl('^{-1}')

    plotfaves, /restore
    
return
end
