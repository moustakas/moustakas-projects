pro sirtfz_event, event
; event handler for the main widget menu

    common sirtf_simulations, sirtf
    common cosmology, cosmo_params
    
    widget_control, event.id, get_uvalue = event_name

    case event_name of

       'Quit': begin
          sirtfz_shutdown, /all
          widget_control, event.top, /destroy
       end
       'CWW (HyperZ)' : temp = dialog_message('This feature has not been implemented yet!',$
                                              /error,dialog_parent=event.id)

;         sirtfz_shutdown, /sedcube 
;         sirtf.templates = 'CWW'
;         sirtf.sedcube = read_sedcube(templates='CWW')
       'CWW (Benitez)': temp = dialog_message('These SED templates are not available yet!',$
                                              /error,dialog_parent=event.id)
       'Devriendt': temp = dialog_message('These SED templates are not available yet!',$
                                              /error,dialog_parent=event.id)
       'GISSEL': temp = dialog_message('These SED templates are not available yet!',$
                                              /error,dialog_parent=event.id)
       'New': temp = dialog_message('This feature has not been implemented yet!',$
                                              /error,dialog_parent=event.id)
       'Select': temp = dialog_message('This feature has not been implemented yet!',$
                                              /error,dialog_parent=event.id)
       'Extract Redshifts': begin
          if strcompress(sirtf.catalog,/remove) eq '' then begin
             temp = dialog_message('Please select a photometric catalog to analyze.',$
                                   /error,dialog_parent=event.id) 
          endif else begin
             fileid = fsc_fileselect(sirtf.baseid,directoryname=filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='results'),$
                                     filename=strmid(sirtf.catalog,0,strpos(sirtf.catalog,'.'))+'.dat',$
                                     /frame,objectref=objectref)
             filename = objectref->getfilename()
             stop
             galz, outname=filename
          endelse
       end
       'Analyze': temp = dialog_message('This feature has not been implemented yet!',$
                                              /error,dialog_parent=event.id)
       else: print, 'Not a valid menu option!'

    endcase


return
end
    
pro sirtfz
;+
; NAME:
;	SIRTFZ
;
; PURPOSE:
;	A general user-interactive widget program to study photometric
;	redshift methods in the infrared.  And other stuff.
;
; COMMON BLOCKS:
;	sirtf_simulations
;	cosmology
;    
; COMMENTS:
;
;
; MODIFICATION HISTORY:
;	John Moustakas, 2001 January
;-

; Notes:  because "sirtf_menu_desc" is an array of structures i can't
; specify independent event handling programs, which would be
; extremely useful for the a "selecting catalog" menu option.

    common sirtf_simulations, sirtf
    common cosmology, cosmo_params

    sirtfz_initialize           ; initialize the common blocks and read in the models

    base = widget_base(title='SIRTFz',/row,/base_align_bottom,uname='sirtf_base',$
                       app_mbar=sirtf_menu)
    sirtf.baseid = base
    
    tmp_struct = {cw_pdmenu_s, flags: 0, name: ''}
    sirtf_menu_desc = [ {cw_pdmenu_s, 1, 'File'},                  $
                        {cw_pdmenu_s, 2, 'Quit'},     $
                        {cw_pdmenu_s, 1, 'SEDs'},        $
                        {cw_pdmenu_s, 3, 'Select'},   $
                        {cw_pdmenu_s, 0, 'CWW (HyperZ)'},$ ; observe each SED through the selected filters
                        {cw_pdmenu_s, 0, 'CWW (Benitez)'},        $
                        {cw_pdmenu_s, 0, 'Devriendt'},  $ ; select an existing photometric catalog
                        {cw_pdmenu_s, 2, 'GISSEL'},       $
                        {cw_pdmenu_s, 1, 'Catalogs'},       $
                        {cw_pdmenu_s, 0, 'New'},  $
                        {cw_pdmenu_s, 2, 'Select'}, $
                        {cw_pdmenu_s, 1, 'SIRTFz'},          $
                        {cw_pdmenu_s, 0, 'Extract Redshifts'}, $
                        {cw_pdmenu_s, 2, 'Analyze'} ]

    sirtf_menu = cw_pdmenu(sirtf_menu,sirtf_menu_desc,/mbar,/help, $ ; main pull-down menu
                           /return_name,uname='sirtf_menu',$
                           font='-adobe-times-bold-r-normal--12-120-75-75-p-67-iso8859-1')
    
; center the top-level base (since widget_info can't figure out the
; size of a menubar widget, we need to guess the geometry)

    device, get_screen_size=xysize

;   geometry = widget_info(base,/geometry)
;   widget_control, base, xoffset=xysize[0]/2L-geometry.scr_xsize/2L, $
;     yoffset=xysize[1]/2L-geometry.scr_ysize/2L    

    widget_control, base, xoff=xysize[0]/2L-170, yoff=xysize[1]/2L-250
    widget_control, base, /realize
    widget_control, sirtf_menu, event_pro = 'sirtfz_event'
    xmanager, 'sirtf_menu', base, event_handler='sirtfz_event', $
      /no_block, group_leader=base;, cleanup='sirtfz_shutdown'

return
end




