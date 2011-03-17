pro compile_sirtfz
; jm01jan26uofa
; this routine is more a debugging tool than anything.  since the main
; widget has been set up as a non-blocking widget, modifications can
; be made to subroutines, then compiled by pressing this menu option.
; that way the user does not need to exit from the main program, and
; all the data don't have to be re-read.  the procedure and function
; names need to be updated as the application is developed.

    routines = ['zmatch','clear_data_catalog','clear_model_catalog',$
                'cosmology_control','create_catalog','evolution_control',$
                'flimit_control','lf_control','likelihood',$
                'observe_each_sed','observe_models','redshift_control',$
                'select_catalog','view_model_seds','waveband_control',$
                'band_flux','define_data_catalog','flux_error','read_bandcube',$
                'read_sedcube','read_lf','sirtf_datapath','filter_match',$
                'create_catalog_options','cz_degeneracy','read_data_catalog',$
                'cmsave']

    resolve_routine, routines, /either

return
end
