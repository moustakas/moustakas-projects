pro phys130_lab_f12_grades, alldata, test=test, sendit=sendit, final=final
; jm12oct05siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/130-F12/grades/'
    
    date = '12dec04' ; update this
    semester = 'Fall 2012'
    class = 'Physics 130 - General Physics I - Lab'
    
; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'phys130_lab_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = 'Lab'
    weight = 1.0

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, /lab, $
      test=test, sendit=sendit, alldata=alldata, final=final
    struct_print, alldata

return
end
