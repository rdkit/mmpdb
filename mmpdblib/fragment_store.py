import psycopg2
import subprocess

db_home="/home/oriol/dev/mmpdb/new_db/"
user="postgres"
password="postgres"
db="mmpdb"
host="localhost"


class FragmentStore:
    def __init__(self,filename):
        self.text_file = open(db_home + filename, "w")
        self.text_file.write("constant_smiles|constant_symmetry_class|num_cuts|id|variable_symmetry_class|variable_smiles|attachment_order|enumeration_label\n")
    def insert(self, constant_smiles, constant_symmetry_class, num_cuts, id, variable_symmetry_class, variable_smiles,
               attachment_order, enumeration_label):
        self.text_file.write(constant_smiles + "|" + constant_symmetry_class + "|" + str(num_cuts) + "|" + str(id) + "|" + variable_symmetry_class + "|" + variable_smiles + "|" + attachment_order + "|" + enumeration_label+"\n")

    def close(self):
        self.text_file.close()

class ConstantStore:
    def __init__(self):
        self.text_file = open(db_home +"/salida_cte.txt", "w")
        self.text_file.write("constant_smiles|constant_with_H_smiles\n")
    def insert(self, constant_smiles, constant_with_H_smiles):
        self.text_file.write(constant_smiles + "|" + constant_with_H_smiles+"\n")

    def close(self):
        self.text_file.close()

class MainStore:
    def __init__(self):
        self.text_file = open(db_home +"salida_global.txt", "w")
        self.text_file.write("normalized_smiles|id\n")

    def insert(self, normalized_smiles, id):
        self.text_file.write(normalized_smiles + "|" + str(id) + "\n")

    def close(self):
        self.text_file.close()

def pg_load():
    myinput = open(db_home+'script.sql')
    p = subprocess.Popen(["psql","-U", user,"-h",host,"-d", db], stdin=myinput)
    p.wait()

def pg_load_sc():
    myinput = open(db_home + 'import_h_cte.sql')
    p = subprocess.Popen(["psql", "-U", user, "-h", host, "-d", db], stdin=myinput)
    p.wait()




def get_ct(conn,smiles):
    cur = conn.cursor()
    cur.execute("""SELECT constant_with_h_smiles from constant_unique where constant_smiles='%s' """ % smiles)
    rows = cur.fetchall()
    if (cur.rowcount>0):
        #print("Row")
        for row in rows:
            #print(row)
            r=row[0]
    else:
        r = None
    return(r)

def get_id(conn,smiles):
    cur = conn.cursor()
    cur.execute("""SELECT id from main where normalized_smiles='%s' """ % smiles)
    rows = cur.fetchall()
    lv=[]
    if (cur.rowcount>0):
        print("Id Row")
        for row in rows:
            print(row)
            lv.append(row[0])
    else:
        lv = []
    return(lv)

def addsc():
    conn = psycopg2.connect("dbname='" + db + "' user='" + user + "' host='" + host + "' password='" + password + "'")
    fstore = FragmentStore("salida_single_cut.txt")
    cur = conn.cursor()
    cur.execute("""SELECT constant_smiles,constant_symmetry_class,num_cuts,variable_part from index_agg where num_cuts='1';""")
    rows = cur.fetchall()
    i=0
    for row in rows:
        i=i+1
        if (i%100==0): print("Row: " + str(i)+"\n")
        smiles=get_ct(conn,row[0])
        list_id=get_id(conn,smiles)
        for other_id in list_id:
            constant_smiles=row[0]
            constant_symmetry_class = row[1]
            num_cuts = row[2]
            js= row[3]
            for v in js:
                fstore.insert(constant_smiles, constant_symmetry_class, num_cuts, v['id'], "1", "[*:1][H]", "0", "N-CTE")

    conn.close()
    fstore.close()

