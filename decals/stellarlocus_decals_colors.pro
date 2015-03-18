function stellarlocus_decals_colors, maggies, filterlist=filterlist
; jm09jul28ucsd - support routine
    case filterlist of
       'grz': begin
          colors = {$
            gr: reform(-2.5*alog10(maggies[0,*]/maggies[1,*])), $
            rz: reform(-2.5*alog10(maggies[1,*]/maggies[2,*]))}
       end
    endcase

return, colors
end
