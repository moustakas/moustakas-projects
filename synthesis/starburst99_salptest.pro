pro salptest

;   tail -2 salptest?.ewidth1
    

    t1 = im_read_starburst99('salptest1.spectrum1',ew='salptest1.ewidth1',qua='salptest1.quanta1')
    t2 = im_read_starburst99('salptest2.spectrum1',ew='salptest2.ewidth1',qua='salptest2.quanta1')
    t3 = im_read_starburst99('salptest3.spectrum1',ew='salptest3.ewidth1',qua='salptest3.quanta1')
    t4 = im_read_starburst99('salptest4.spectrum1',ew='salptest4.ewidth1',qua='salptest4.quanta1')
    t5 = im_read_starburst99('salptest5.spectrum1',ew='salptest5.ewidth1',qua='salptest5.quanta1')

    
    djs_plot, t1.time/1E6, t1.ha, psym=-6, xsty=3, ysty=3, yr=[40.5,42]
    djs_oplot, t2.time/1E6, t2.ha, psym=-6, color='red'
    djs_oplot, t3.time/1E6, t3.ha, psym=-6, color='blue'
    djs_oplot, t4.time/1E6, t4.ha, psym=-6, color='green'
    djs_oplot, t5.time/1E6, t5.ha, psym=-6, color='orange'

stop    
    
return
end
    
