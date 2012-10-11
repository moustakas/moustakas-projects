pro phys390_f12_grades
; jm12oct05siena - parse the grades for this class 

    date = '12oct05' ; update this
    semester = 'Fall 2012'
    class = 'Physics 390: Introductory Astrophysics I'
    
; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = '~/Dropbox/teaching/12fall/phys390/grades/phys390_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)
    nstudent = n_elements(data)

; as a special case for this class alone, assign everybody 100% on
; their participation grade
    assign = [assign,'Participation']
    part = replicate({participation: 100.0, participation_details: 'APOD Discussion', $
      participation_date: 'various dates', participation_points: 100.0},nstudent)
    data = struct_addtags(data,part)
    
; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Participation','Homework','Midterm',$
      'Research Paper','Research Presentation','Final']
    weight = [0.05,0.4,0.1,0.15,0.15,0.15]

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester

stop    

return
end
