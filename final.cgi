#!/usr/bin/env python3

import jinja2
import cgi
import re
import mysql.connector



# This line tells the template loader where to search for template files
templateLoader = jinja2.FileSystemLoader( searchpath="./templates" )

# This creates your environment and loads a specific template
env = jinja2.Environment(loader=templateLoader)
template = env.get_template('final.html')

form = cgi.FieldStorage()
oligo_seq = str(form.getfirst('keyword'))
#seq_list=list()
#reverse_list=list()
#Length_list=list()
#GC_list=list()
#tm_list=list()


def reverse_complement(oligo_seq):
    oligo_join=(','.join(oligo_seq))
    oligo_seq_split=oligo_join.split(',')
    oligo_seq_split.reverse()
    n=len(oligo_seq)

    for x in range(0, n):
        if oligo_seq_split[x]=="C":
            oligo_seq_split[x]="G"
        elif oligo_seq_split[x]=="G":
            oligo_seq_split[x]="C"
        elif oligo_seq_split[x]=="A":
            oligo_seq_split[x]="T"
        elif oligo_seq_split[x]=="T":
            oligo_seq_split[x]="A"
    str = ""
    oligo_seq_rev = str.join( oligo_seq_split )
    return oligo_seq_rev

def GC_content(oligo_seq):
    G = oligo_seq.count("G")
    C = oligo_seq.count("C")
    GC_content=((G+C)/len(oligo_seq))*100
    return GC_content

def Tm(oligo_seq):
    G = oligo_seq.count("G")
    C = oligo_seq.count("C")
    A = oligo_seq.count("A")
    T = oligo_seq.count("T")
    if len(oligo_seq)> 14:
        Tm= 64.9 +41*(G+C-16.4)/(A+T+G+C)
    else:
        Tm= (A+T)*2 + (G+C)*4
    return Tm


if re.search(r"[^ATGC]", oligo_seq):
    seq=oligo_seq
    reverse="An invalid character was found!"
    Length="------------"
    GC="------------"
    tm="------------"
else:
    seq=oligo_seq
    reverse=reverse_complement(oligo_seq)
    Length=len(oligo_seq)
    GC=('%.1f' % GC_content(oligo_seq))
    tm=('%.1f' % Tm(oligo_seq))


conn = mysql.connector.connect(user='hzhan141', password='huiminzhang',
                                   host='localhost', database='hzhan141_chado')
curs = conn.cursor()
    #curs.execute('''CREATE TABLE history
    #       (Sequence TEXT NOT NULL,
     #      Reverse TEXT NOT NULL, 
    #       Length INTEGER NOT NULL, 
     #      GC INTEGER NOT NULL, 
     #      Melting_temperature INTEGER NOT NULL)''')
qry = "INSERT INTO history(Sequence, Reverse, Length, GC, Melting_temperature) VALUES (%s, %s, %s, %s, %s)"
curs.execute(qry, (seq, reverse, Length, GC, tm))
conn.commit()
curs.close()



print("Content-type: text/html\n\n")
print(template.render(seq=seq, reverse=reverse,Length=Length, GC=GC, tm=tm))



#print("Content-type: text/html\n\n")
#print(template.render(seq=seq, reverse=reverse,Length=Length, GC=GC, tm=tm ))


#try :
    #conn = sqlite3.connect("/var/www/html/hzhan141/final/history.db")
    #curs = conn.cursor()
    #qry = "INSERT INTO history (Sequence, Reverse, Length, GC, Melting_temperature) VALUES (?, ?, ?, ?, ?)"
    #curs.execute(qry, (seq, reverse, Length, GC, tm))  
    #conn.commit()
   # conn.close()
#except:
   # e = sys.exc_info()[0]
   # print("Content-Type: text/plain\n")
    #print(e)



#print("Content-type: text/html\n\n")
#print(template.render(seq=seq, reverse=reverse,Length=Length, GC=GC, tm=tm ))


#conn = sqlite3.connect('history.db')
#try:
#    conn.execute('''CREATE TABLE history
 #          (Sequence TEXT NOT NULL,
  #         Reverse TEXT NOT NULL, 
   #        Length INTEGER NOT NULL, 
    #       GC INTEGER NOT NULL, 
     #      Melting_temperature INTEGER NOT NULL)''')
    #print ("Table created successfully")
#except:
 #   pass

#conn.close()






#def main():
 #   conn = sqlite3.connect('history.db')
 #   curs = conn.cursor()
   
 #   qry ="INSERT INTO history (Sequence, Reverse, Length, GC, Melting_temperature) VALUES (?, ?, ?, ?, ?)"
 #   curs.execute(qry, (seq, reverse, Length, GC, tm))

#    qry_2="select * from history"
#    curs.execute(qry_2)
  #  curs.commit()
   # curs.close()

#if __name__ == '__main__':
 #   main()

#print("Content-type: text/html\n\n")

