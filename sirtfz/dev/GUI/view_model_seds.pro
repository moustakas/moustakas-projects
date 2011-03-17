pro sed_window_event, event

    common sirtf_simulations

    widget_control, event.top, get_uvalue=sed_data, /no_copy
    widget_control, event.id, get_uvalue=event_name

    wset, sed_data.windx

    case event_name of
       'overplot filters': begin

          scale = (sed_data.fmax-sed_data.fmin)/10.0
    
          indx = filter_match(*sirtf.filters) ; currently selected filter bandpasses
          nmatch = n_elements(indx)

          for k = 0L, nmatch-1L do $
            oplot, *(sirtf.bandcube[indx[k]].wband), *(sirtf.bandcube[indx[k]].rband)+sed_data.fmin, $
            line=2, thick=2.5

          stop
       end
       'overplot data': 
       'write ps':
       'help':
    endcase    

    widget_control, event.top, set_uvalue=sed_data, /no_copy

return
end

pro sed_window, wave_z, flux_z, ytitle=ytitle, title=title, info

    common sirtf_simulations

    widget_control, info.child_id, /realize
    xmanager, 'sed_window', info.child_id, event_handler='sed_window_event', /no_block
          
    widget_control, info.draw_id, get_value=window_index
    wset, window_index
    plot, wave_z, flux_z, xsty=3, ysty=3, /xlog, /ylog, xtit=textoidl('\lambda')+' ('+$
      textoidl('\mu')+'m)', ytitle=ytitle, title=title, xrange=[0.1,10000]

; store some SED data in an internal structure and recover it in the
; event handler    

    widget_control, info.draw_id, get_value = window_index

    sed_data = {windx: window_index, fmax: max(flux_z), fmin: min(flux_z)}
    widget_control, info.child_id, set_uvalue = sed_data, /no_copy
    
return
end

pro view_seds_event, event

    common sirtf_simulations
    
    widget_control, event.top, get_uvalue = info, /no_copy
    event_name = tag_names(event,/structure_name)

    case event_name of

       'FSC_FIELD': info.sedz = *event.value ; new redshift

       'FSC_DROPLIST_EVENT': begin ; droplist options

          indx = event.index
          if info.sedz gt 0.0 then begin
             redshift_sed, info.sedz, *sirtf.sedcube[indx].lambda, *sirtf.sedcube[indx].mlum, $ ; redshift the SED
               wave_z, flux_z, /jansky
             ytitle = textoidl('f_{\nu}')+' (mJy)'
          endif else begin
             wave_z = *sirtf.sedcube[indx].lambda
             flux_z = *sirtf.sedcube[indx].mlum
             ytitle = textoidl('L_{\nu}')+' (W Hz'+textoidl('^{-1}')+')'
          endelse

          sed_window, wave_z, flux_z, ytitle=ytitle, title=sirtf.sedcube[indx].galaxy, info ; plotit

       end

       'WIDGET_BUTTON': begin   ; button option
          widget_control, event.id, get_uvalue=option
          if option eq 'done' then widget_control, event.top, /destroy
       end
       else: temp = dialog_message('There was a problem here!',/error,dialog_parent=event.top)
    endcase

return
end

pro view_model_seds
; jm01jan10uofa
; examine the Devriendt model SEDs

    common sirtf_simulations
    common cosmology

    gals = string(sirtf.sedcube.galaxy,format='(A9)')    ; model galaxy names
    gtypes = sirtf.sedcube.gtype   ; model galaxy types
    lbol = string(sirtf.sedcube.lum60_sun,format='(F5.2)') ; 60mu luminosity (solar units)

    list = gals+' ('+strcompress(lbol,/remove)+')'

; parent base: contains the list and redshift of the SEDs

    base = widget_base(title='SED Models',/row,/align_bottom,group_leader=sirtf.widget_ids.base_id,$
                       /base_align_right)
    sed_list = fsc_droplist(base,event_pro='view_seds_event',index=0,uvalue='sed_list',$
                            value=[list[0],list[1],list[2],list[3],list[4],list[5],list[6],$
                                   list[7],list[8],list[9],list[10],list[11],list[12],list[13],$
                                   list[14],list[15]],labelsize=90)
    field = fsc_inputfield(base,title='Redshift:',event_pro='view_seds_event',labelsize=60, $
                           uvalue='redshift',value=0.0,/cr_only,/positive,decimal=6,$
                           xsize=5,/doublevalue)
    button = widget_button(base,value='Done',uvalue='done',units=1,xsize=1.0)

; draw window base in which the SEDs are plotted (realize these
; widgets only once an SED has been selected)
    
    dbase = widget_base(title='Spectral Energy Distribution Model',/column,$
                        /base_align_center,group_leader=base)
    draw = widget_draw(dbase,scr_xsize=620L,scr_ysize=540L)

    button_base = widget_base(dbase,/row,/base_align_bottom)
    button = widget_button(button_base,value='Overplot Filters',uvalue='overplot filters',$
                           units=1,xsize=1.5)
    button = widget_button(button_base,value='Overplot Data',uvalue='overplot data',$
                           units=1,xsize=1.5)
    button = widget_button(button_base,value='Write PS',uvalue='write ps',units=1,xsize=1.5)
    button = widget_button(button_base,value='Help',uvalue='help',units=1,xsize=1.5)
    
; fill the information structure

    info = {base_id: base, child_id: dbase, draw_id: draw, sedz: float(0.0)}

; realize the parent base
    
    widget_control, base, /realize, set_uvalue = info, /no_copy
    xmanager, 'view_model_seds', base, event_handler = 'view_seds_event', /no_block
    
return
end

