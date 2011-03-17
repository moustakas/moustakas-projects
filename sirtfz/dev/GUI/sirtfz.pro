pro sirtfz_shutdown, windowid
; destroy the widgets, kill the windows, and free the pointers

    common sirtf_simulations, sirtf

    if (xregistered('sirtfz')) then widget_control, sirtf.widgets.base_id, /destroy

    clear_data_catalog
    clear_model_catalog
    ptr_free, sirtf.bandcube[*].rband
    ptr_free, sirtf.bandcube[*].wband
    ptr_free, sirtf.sedcube[*].lambda
    ptr_free, sirtf.sedcube[*].mlum
    
    if ptr_valid(sirtf.redshift.zarray) then ptr_free, sirtf.redshift.zarray
    if ptr_valid(sirtf.filters) then ptr_free, sirtf.filters

    if obj_valid(progress) then obj_destroy, progress

;   wkill
    
    plotfaves, /restore

return
end

pro sirtfz_event, event
; event handler for the main widget menu

    common sirtf_simulations, sirtf
    common cosmology, cosmo_params
    
    widget_control, event.id, get_uvalue = event_name

    case event_name of

       'Compile': compile_sirtfz
       'Quit': begin
          sirtfz_shutdown
          widget_control, event.top, /destroy
       end
; ----------------------------------------------------------------------
       'Parameters'            : cosmology_control
       'Evolution'             : evolution_control
; ----------------------------------------------------------------------
       'Filters'               : waveband_control
       'Limiting Flux'         : flimit_control
       'Redshift Grid'         : redshift_control
       'Luminosity Function'   : lf_control
; ----------------------------------------------------------------------
       'Plot SEDs'             : view_model_seds
       'Observe SEDs'          : observe_models
       'Color-Color'           : cz_degeneracy
       'Color-Redshift'        : cz_degeneracy
; ----------------------------------------------------------------------
       'Display Catalog'       : temp = dialog_message('This feature has not been implemented yet!',$
                                                       /error,dialog_parent=event.id)
; ----------------------------------------------------------------------
       'Create Catalog'        : create_catalog
       'Display Image'         : temp = dialog_message('This feature has not been implemented yet!',$
                                                       /error,dialog_parent=event.id)
       'Calibrate Priors'      : temp = dialog_message('This feature has not been implemented yet!',$
                                                       /error,dialog_parent=event.id)
       'Extract Redshifts' : zmatch

       else: begin ; was a photometric catalog selected?

          if (event_name eq 'HDF North') or (event_name eq 'HDF South') or $
            (event_name eq 'NTT Deep Field') or (event_name eq 'NOAO Deep Field') then begin

             clear_data_catalog
             temp = read_data_catalog(event_name)

          endif else print, 'Not a valid menu option!'

       end
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

    sirtfz_initialize ; initialize the common blocks and read in the models

    base = widget_base(title='SIRTF',/row,/base_align_bottom,uname='sirtf_base',$
                       app_mbar=sirtf_menu)
    sirtf.widget_ids.base_id = base
    
    tmp_struct = {cw_pdmenu_s, flags: 0, name: ''}
    sirtf_menu_desc = [ {cw_pdmenu_s, 1, 'File'},                  $
                        {cw_pdmenu_s, 0, 'Compile'},  $ ; compile all procedures
                        {cw_pdmenu_s, 2, 'Quit'},     $
                        {cw_pdmenu_s, 1, 'Cosmology'}, $ 
                        {cw_pdmenu_s, 0, 'Parameters'},$ ; cosmological parameters
                        {cw_pdmenu_s, 2, 'Evolution'}, $ ; luminosity/density evolution
                        {cw_pdmenu_s, 1, 'Input'},              $
                        {cw_pdmenu_s, 0, 'Filters'},            $ ; observing filters
                        {cw_pdmenu_s, 0, 'Limiting Flux'},      $ ; limiting flux in each band
                        {cw_pdmenu_s, 0, 'Redshift Grid'},      $ ; redshift min, max and binsize
                        {cw_pdmenu_s, 2, 'Luminosity Function'},$ ; display the luminosity function
                        {cw_pdmenu_s, 1, 'SEDs'},        $
                        {cw_pdmenu_s, 0, 'Plot SEDs'},   $ ; display a model SED
                        {cw_pdmenu_s, 0, 'Observe SEDs'},$ ; observe each SED through the selected filters
                        {cw_pdmenu_s, 0, 'Color-Color'}, $
                        {cw_pdmenu_s, 2, 'Color-Redshift'},$   ; generate color-redshift degeneracy plots
                        {cw_pdmenu_s, 1, 'Database'},        $
                        {cw_pdmenu_s, 1, 'Select Catalog'},  $ ; select an existing photometric catalog
                        {cw_pdmenu_s, 0, 'HDF North'},       $
                        {cw_pdmenu_s, 0, 'HDF South'},       $
                        {cw_pdmenu_s, 0, 'NTT Deep Field'},  $
                        {cw_pdmenu_s, 2, 'NOAO Deep Field'}, $
                        {cw_pdmenu_s, 2, 'Display Catalog'}, $ ; display an image of the catalog on the sky
                        {cw_pdmenu_s, 1, 'Simulations'},          $
                        {cw_pdmenu_s, 0, 'Create Catalog'},       $
                        {cw_pdmenu_s, 0, 'Display Image'},        $
                        {cw_pdmenu_s, 1, 'Photometric Redshifts'},$
                        {cw_pdmenu_s, 0, 'Calibrate Priors'},    $
                        {cw_pdmenu_s, 2, 'Extract Redshifts'} ]

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
      cleanup='sirtfz_shutdown', /no_block, group_leader=base

return
end




