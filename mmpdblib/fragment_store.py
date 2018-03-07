import psycopg2
import subprocess

db_home="/home/oriol/dev/mmpdb/new_db/"
filename="compound_test"
user="postgres"
password="postgres"
db="mmpdb"
host="localhost"
filename_index =db_home + filename + ".txt"
filename_cte = db_home + filename + "_cte.txt"
filename_main = db_home + filename + "_main.txt"
filename_sc = db_home +filename + "_single_cut.txt"
filename_idrecord = db_home +filename + "_idrecord.txt"

class FragmentStore:
    def __init__(self,filename):
        self.text_file = open(filename, "w")
        self.text_file.write("constant_smiles|constant_symmetry_class|num_cuts|id|variable_symmetry_class|variable_smiles|attachment_order|enumeration_label\n")
    def insert(self, constant_smiles, constant_symmetry_class, num_cuts, id, variable_symmetry_class, variable_smiles,
               attachment_order, enumeration_label):
        self.text_file.write(constant_smiles + "|" + constant_symmetry_class + "|" + str(num_cuts) + "|" + str(id) + "|" + variable_symmetry_class + "|" + variable_smiles + "|" + attachment_order + "|" + enumeration_label+"\n")

    def close(self):
        self.text_file.close()

class ConstantStore:
    def __init__(self):
        self.text_file = open(filename_cte, "w")
        self.text_file.write("constant_smiles|constant_with_H_smiles\n")
    def insert(self, constant_smiles, constant_with_H_smiles):
        self.text_file.write(constant_smiles + "|" + constant_with_H_smiles+"\n")

    def close(self):
        self.text_file.close()

class MainStore:
    def __init__(self):
        self.text_file = open(filename_main, "w")
        self.text_file.write("normalized_smiles|id\n")

    def insert(self, normalized_smiles, id):
        self.text_file.write(normalized_smiles + "|" + str(id) + "\n")

    def close(self):
        self.text_file.close()

class IdRecordStore:
    def __init__(self):
        self.text_file = open(filename_idrecord, "w")
        self.text_file.write("id|input_smiles|num_normalized_heavies|normalized_smiles\n")

    def insert(self,id, input_smiles,num_normalized_heavies,normalized_smiles):
        self.text_file.write(id+"|"+input_smiles+"|"+str(num_normalized_heavies)+"|"+normalized_smiles+"\n")

    def close(self):
        self.text_file.close()

def pg_create_tables():
    print("Creating tables\n")
    myinput = open(db_home+'create_tables.sql')
    p = subprocess.Popen(["psql","-U", user,"-h",host,"-d", db], stdin=myinput)
    p.wait()

def pg_load():
    print("Loading CSV\n")
    conn = psycopg2.connect("dbname='" + db + "' user='" + user + "' host='" + host + "' password='" + password + "'")
    print("C:"+str(conn))
    #import INDEX
    cur=conn.cursor()
    print(filename_index+ "\n")
    file_object= open(filename_index)
    SQL_STATEMENT = """
        COPY index(constant_smiles,constant_symmetry_class,num_cuts,id,variable_symmetry_class,variable_smiles,attachment_order,enumeration_label)
        FROM STDIN WITH
            CSV
            HEADER
            DELIMITER AS '|'
        """
    r=cur.copy_expert(sql=SQL_STATEMENT, file=file_object)
    print("R: "+str(r))
    print(SQL_STATEMENT)
    conn.commit()
    cur.close()
    #import constant
    cur=conn.cursor()
    file_object= open(filename_cte)
    SQL_STATEMENT = """
        COPY constant(constant_smiles,constant_with_H_smiles)
        FROM STDIN WITH
            CSV
            HEADER
            DELIMITER AS '|'
        """
    r=cur.copy_expert(sql=SQL_STATEMENT, file=file_object)
    conn.commit()
    cur.close()
    #import main
    cur=conn.cursor()
    file_object= open(filename_main)
    SQL_STATEMENT = """
        COPY main(normalized_smiles,id)
        FROM STDIN WITH
            CSV
            HEADER
            DELIMITER AS '|'
        """
    r=cur.copy_expert(sql=SQL_STATEMENT, file=file_object)
    conn.commit()
    cur.close()
    #import idrecord
    cur=conn.cursor()
    file_object= open(filename_idrecord)
    SQL_STATEMENT = """
        COPY idrecord(id,input_smiles,num_normalized_heavies,normalized_smiles)
        FROM STDIN WITH
            CSV
            HEADER
            DELIMITER AS '|'
        """
    r=cur.copy_expert(sql=SQL_STATEMENT, file=file_object)
    conn.commit()
    cur.close()
    conn.close()

def pg_transform():
    print("Transforming\n")
    #tranform
    myinput = open(db_home + 'scripts.sql')
    p = subprocess.Popen(["psql", "-U", user, "-h", host, "-d", db], stdin=myinput)
    p.wait()

def pg_reagg():
    print("Reagg\n")
    #reagg
    myinput = open(db_home + 'script_reagg.sql')
    p = subprocess.Popen(["psql", "-U", user, "-h", host, "-d", db], stdin=myinput)
    p.wait()

def pg_load_sc():
    print("Load SC\n")
    #myinput = open(db_home + 'import_h_cte.sql')
    #p = subprocess.Popen(["psql", "-U", user, "-h", host, "-d", db], stdin=myinput)
    #p.wait()
    conn = psycopg2.connect("dbname='" + db + "' user='" + user + "' host='" + host + "' password='" + password + "'")
    # import sc
    cur = conn.cursor()
    file_object = open(filename_sc)
    SQL_STATEMENT = """
           COPY index(constant_smiles,constant_symmetry_class,num_cuts,id,variable_symmetry_class,variable_smiles,attachment_order,enumeration_label)
           FROM STDIN WITH
               CSV
               HEADER
               DELIMITER AS '|'
           """
    r = cur.copy_expert(sql=SQL_STATEMENT, file=file_object)
    conn.commit()
    cur.close()
    conn.close()

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
        for row in rows:
            lv.append(row[0])
    else:
        lv = []
    return(lv)

def addsc():
    print("Add SC\n")
    ## Add the single cut hydrogen transformations

    # The algorithm is:
    #   - for each single cut constant, get its with-hydrogen version
    #   - if the with-hydrogen version matches an actual record
    #   - add the records using the [*:1][H] variable fragment
    #
    # for (constant_smiles, constant_symmetry_class, num_cuts), matches in index.items():
    #     print("Single cut:")
    #     print(constant_smiles,constant_symmetry_class,num_cuts)
    #     if num_cuts != 1:
    #         continue
    #     constant_with_H_smiles = constant_smiles_to_hydrogen_constant_smiles.get(constant_smiles, None)
    #     if constant_with_H_smiles is None:
    #         continue
    #     other_ids = normalized_smiles_to_ids.get(constant_with_H_smiles, [])
    #     for other_id in other_ids:
    #         # NOTE: this is hard-coded to "[*:1][H]", and must match the
    #         # same string used in fragment.py's _hydrogen_cut_smiles
    #         print("--------------")
    #         print(constant_smiles, constant_symmetry_class, num_cuts)
    #         print(" -> ")
    #         print(other_id, "1", "[*:1][H]", "0", "N")
    #         print("--------------")
    #         fstore.insert(constant_smiles,
    #                       constant_symmetry_class,
    #                       num_cuts,
    #                       other_id, "1", "[*:1][H]", "0", "N-CTE")
    #         matches.append((other_id, "1", "[*:1][H]", "0", "N"))
    conn = psycopg2.connect("dbname='" + db + "' user='" + user + "' host='" + host + "' password='" + password + "'")
    fstore = FragmentStore(filename_sc)
    cur = conn.cursor()
    cur.execute("""SELECT constant_smiles,constant_symmetry_class,num_cuts,variable_part from index_agg where num_cuts='1';""")
    rows = cur.fetchall()
    i=0
    for row in rows:
        i=i+1
        #if (i%1000==0): print("Row: " + str(i)+"\n")
        smiles=get_ct(conn,row[0])
        if (smiles is None):
            continue
        list_id=get_id(conn,smiles)
        for other_id in list_id:
            print("Id: " + other_id+"\n")
            constant_smiles=row[0]
            constant_symmetry_class = row[1]
            num_cuts = row[2]
            js= row[3]
            for v in js:
                print("------------------\n")
                print("DBRow: " + str(i))
                print(row)
                print("\t" + str(js) + "\n")
                print("\t" + str(v)+"\n")
                print("Iiiiiiiiiiiiiiiiiiiiiiiiiiiii\n")
                print("\t" +constant_smiles +"|"+ constant_symmetry_class+"|"+num_cuts+"|"+other_id+"|"+"1"+"|"+ "[*:1][H]"+"|"+"0"+"|"+"N-CTE" + "\n")
                print("------------------\n")
                fstore.insert(constant_smiles, constant_symmetry_class, num_cuts, other_id, "1", "[*:1][H]", "0", "N-CTE")

    conn.close()
    fstore.close()

def dump(index,id_to_record):
    text_file = open(db_home +"index_mem.txt", "w")
    text_file.write("constant_smiles|constant_symmetry_class|num_cuts)|idr|variable_symmetry_class|variable_smiles|attachment_order|enumeration_label\n")
    for (constant_smiles, constant_symmetry_class, num_cuts), matches in sorted(index.items()):
        #print(constant_smiles, constant_symmetry_class, num_cuts, constant_smiles, constant_symmetry_class, matches)
        #print(type(matches))
        #print(len(matches))
        for idr, variable_symmetry_class, variable_smiles, attachment_order, enumeration_label in matches:
            text_file.write(constant_smiles +"|"+ constant_symmetry_class +"|"+str(num_cuts)+"|"+ idr+"|"+ variable_symmetry_class+"|"+ variable_smiles+"|"+ attachment_order+"|"+ enumeration_label+"\n")
    text_file.close()

    text_file = open(db_home + "index_id_mem.txt", "w")
    text_file.write(
        "id|input_smiles|num_normalized_heavies,normalized_smiles\n")
    for idr,record in id_to_record.items():
        text_file.write(record.id+"|"+  record.input_smiles+"|"+  str(record.num_normalized_heavies)+"|"+  record.normalized_smiles+"\n")
    text_file.close()