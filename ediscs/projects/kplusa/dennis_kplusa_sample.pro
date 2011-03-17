pro dennis_kplusa_sample
; jm08jun20nyu - write out some measurements for Dennis in his attempt
;                to reproduce Bianca's qualitative results

    aa = read_ediscs(/ancillary)
    ss = read_ediscs(/specfit)

    sel = struct_addtags(struct_trimtags(aa,select=['cluster',$
      'cluster_z','cluster_sigma','memberflag']),$
      struct_trimtags(ss,select=['galaxy',$
      'z_obj','oii_3727_ew','lick_hd_a_raw',$
      'lick_hg_a_raw']))
    struct_print, sel, file='ediscs_kplusa.dat'
      
return
end
    
