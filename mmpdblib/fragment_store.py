import subprocess
import psycopg2
import psycopg2.extras

#db_home = "/home/oriol/dev/mmpdb/test_db/"
#filename = "compound_test"
#db = "mmpdb_test_new"

db_home = "/home/oriol/dev/mmpdb/compound_all_db/"
filename = "compound_all_new"
db = "mmpdb_all_new"

#db_home = "/home/oriol/dev/mmpdb/compound_100k/"
#filename = "compound_100k_new"
#db = "mmpdb_100k_new"

db_scripts="/home/oriol/dev/mmpdb/new_db_scripts/"
user = "postgres"
password = "postgres"
host = "localhost"
filename_index = db_home + filename + ".txt"
filename_cte = db_home + filename + "_cte.txt"
filename_main = db_home + filename + "_main.txt"
filename_sc = db_home + filename + "_single_cut.txt"
filename_idrecord = db_home + filename + "_idrecord.txt"
filename_rejected = db_home + filename + "_rejected.txt"


class FragmentStore:
    """
    """
    def __init__(self, filename):
        """

        Args:
            filename:
        """
        self.text_file = open(filename, "w")
        self.text_file.write(
            "constant_smiles|constant_symmetry_class|num_cuts|id|variable_symmetry_class|variable_smiles|attachment_order|enumeration_label\n")

    def insert(self, constant_smiles, constant_symmetry_class, num_cuts, id, variable_symmetry_class, variable_smiles,
               attachment_order, enumeration_label):
        """

        Args:
            constant_smiles:
            constant_symmetry_class:
            num_cuts:
            id:
            variable_symmetry_class:
            variable_smiles:
            attachment_order:
            enumeration_label:
        """
        self.text_file.write(constant_smiles + "|" + constant_symmetry_class + "|" + str(num_cuts) + "|" + str(
            id) + "|" + variable_symmetry_class + "|" + variable_smiles + "|" + attachment_order + "|" + enumeration_label + "\n")

    def close(self):
        """

        """
        self.text_file.close()


class ConstantStore:
    """
    """
    def __init__(self):
        """

        """
        self.text_file = open(filename_cte, "w")
        self.text_file.write("constant_smiles|constant_with_H_smiles\n")

    def insert(self, constant_smiles, constant_with_H_smiles):
        """

        Args:
            constant_smiles:
            constant_with_H_smiles:
        """
        self.text_file.write(constant_smiles + "|" + constant_with_H_smiles + "\n")

    def close(self):
        """

        """
        self.text_file.close()


class MainStore:
    """
    """
    def __init__(self):
        """

        """
        self.text_file = open(filename_main, "w")
        self.text_file.write("normalized_smiles|id\n")

    def insert(self, normalized_smiles, id):
        """

        Args:
            normalized_smiles:
            id:
        """
        self.text_file.write(normalized_smiles + "|" + str(id) + "\n")

    def close(self):
        """

        """
        self.text_file.close()


class IdRecordStore:
    """
    """
    def __init__(self):
        """

        """
        self.text_file = open(filename_idrecord, "w")
        self.text_file.write("id|input_smiles|num_normalized_heavies|normalized_smiles\n")

    def insert(self, id, input_smiles, num_normalized_heavies, normalized_smiles):
        """

        Args:
            id:
            input_smiles:
            num_normalized_heavies:
            normalized_smiles:
        """
        self.text_file.write(
            id + "|" + input_smiles + "|" + str(num_normalized_heavies) + "|" + normalized_smiles + "\n")

    def close(self):
        """

        """
        self.text_file.close()


class RejectedStore:
    """
    """
    def __init__(self):
        """

        """
        self.text_file = open(filename_rejected, "w")
        self.text_file.write("number|id|input_smiles\n")

    def insert(self, i, rid, input_smiles):
        """

        Args:
            i:
            rid:
            input_smiles:
        """
        self.text_file.write(i + "|" + rid + "|" + input_smiles + "\n")

    def close(self):
        """

        """
        self.text_file.close()


def pg_create_tables():
    """

    """
    print("Creating tables\n")
    myinput = open(db_scripts + 'create_tables.sql')
    p = subprocess.Popen(["psql", "-U", user, "-h", host, "-d", db], stdin=myinput)
    p.wait()


def pg_load():
    """

    """
    print("Loading CSVs\n")
    print("Loading index\n")
    conn = psycopg2.connect("dbname='" + db + "' user='" + user + "' host='" + host + "' password='" + password + "'")
    # import INDEX
    cur = conn.cursor()
    file_object = open(filename_index)
    SQL_STATEMENT = """
        COPY index(constant_smiles,constant_symmetry_class,num_cuts,id,variable_symmetry_class,variable_smiles,attachment_order,enumeration_label)
        FROM STDIN WITH
            CSV
            HEADER
            DELIMITER AS '|'
        """
    r = cur.copy_expert(sql=SQL_STATEMENT, file=file_object)
    print(SQL_STATEMENT)
    conn.commit()
    cur.close()
    # import constant
    print("Loading constant\n")
    cur = conn.cursor()
    file_object = open(filename_cte)
    SQL_STATEMENT = """
        COPY constant(constant_smiles,constant_with_H_smiles)
        FROM STDIN WITH
            CSV
            HEADER
            DELIMITER AS '|'
        """
    r = cur.copy_expert(sql=SQL_STATEMENT, file=file_object)
    conn.commit()
    cur.close()
    # import main
    print("Loading main\n")
    cur = conn.cursor()
    file_object = open(filename_main)
    SQL_STATEMENT = """
        COPY main(normalized_smiles,id)
        FROM STDIN WITH
            CSV
            HEADER
            DELIMITER AS '|'
        """
    r = cur.copy_expert(sql=SQL_STATEMENT, file=file_object)
    conn.commit()
    cur.close()
    # import idrecord
    print("Loading idrecord\n")
    cur = conn.cursor()
    file_object = open(filename_idrecord)
    SQL_STATEMENT = """
        COPY idrecord(id,input_smiles,num_normalized_heavies,normalized_smiles)
        FROM STDIN WITH
            CSV
            HEADER
            DELIMITER AS '|'
        """
    r = cur.copy_expert(sql=SQL_STATEMENT, file=file_object)
    conn.commit()
    cur.close()
    conn.close()


def pg_transform():
    """

    """
    print("Transforming\n")
    # tranform
    myinput = open(db_scripts + 'scripts.sql')
    p = subprocess.Popen(["psql", "-U", user, "-h", host, "-d", db], stdin=myinput)
    p.wait()


def pg_reagg():
    """

    """
    print("Reagg\n")
    # reagg
    myinput = open(db_scripts + 'script_reagg.sql')
    p = subprocess.Popen(["psql", "-U", user, "-h", host, "-d", db], stdin=myinput)
    p.wait()


def pg_load_sc():
    """

    """
    print("Load Single cut ct\n")
    # myinput = open(db_home + 'import_h_cte.sql')
    # p = subprocess.Popen(["psql", "-U", user, "-h", host, "-d", db], stdin=myinput)
    # p.wait()
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


def get_ct(conn, smiles):
    """

    Args:
        conn:
        smiles:

    Returns:

    """
    cur = conn.cursor()
    cur.execute("""SELECT constant_with_h_smiles from constant_unique where constant_smiles='%s' """ % smiles)
    rows = cur.fetchall()
    if (cur.rowcount > 0):
        # print("Row")
        for row in rows:
            # print(row)
            r = row[0]
    else:
        r = None
    return (r)


def get_id(conn, smiles):
    """

    Args:
        conn:
        smiles:

    Returns:

    """
    cur = conn.cursor()
    cur.execute("""SELECT id from main where normalized_smiles='%s' """ % smiles)
    rows = cur.fetchall()
    lv = []
    if (cur.rowcount > 0):
        for row in rows:
            lv.append(row[0])
    else:
        lv = []
    return (lv)


def addsc():
    """

    """
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
    cur = conn.cursor(name="index_agg_cursor")
    cur.itersize = 50000
    cur.execute(
        """SELECT constant_smiles,constant_symmetry_class,num_cuts,variable_part from index_agg where num_cuts='1';""")
    rows = cur.fetchall()
    i = 0
    for row in rows:
        i = i + 1
        if (i % 1000 == 0): print("Single cut Row: " + str(i))
        smiles = get_ct(conn, row[0])
        if (smiles is None):
            continue
        list_id = get_id(conn, smiles)
        for other_id in list_id:
            # print("Id: " + other_id+"\n")
            constant_smiles = row[0]
            constant_symmetry_class = row[1]
            num_cuts = row[2]
            js = row[3]
            for v in js:
                # print("------------------\n")
                # print("DBRow: " + str(i))
                # print(row)
                # print("\t" + str(js) + "\n")
                # print("\t" + str(v)+"\n")
                # print("Iiiiiiiiiiiiiiiiiiiiiiiiiiiii\n")
                # print("\t" +constant_smiles +"|"+ constant_symmetry_class+"|"+num_cuts+"|"+other_id+"|"+"1"+"|"+ "[*:1][H]"+"|"+"0"+"|"+"N-CTE" + "\n")
                # print("------------------\n")
                fstore.insert(constant_smiles, constant_symmetry_class, num_cuts, other_id, "1", "[*:1][H]", "0", "N")

    conn.close()
    fstore.close()

class FragmentIndexDB(object):
    """
    """
    def __init__(self):
        """

        """
        print("Open Id to record iteration")
        self.conn = psycopg2.connect("dbname='" + db + "' user='" + user + "' host='" + host + "' password='" + password + "'")
        self.cur = self.conn.cursor(cursor_factory=psycopg2.extras.NamedTupleCursor)
        self.cur.itersize = 100000
        self.cur.execute("""SELECT id,input_smiles,normalized_smiles,num_normalized_heavies from idrecord limit 1000;""")
        #self.rows = self.cur.fetchall()

    def close(self):
        """

        """
        print("Closing IX iteration")
        self.conn.close()

    def __len__(self):
        """

        """
        raise("FragmentIndexDB len error")

    def test_iter(self):
        """

        """
        for i in self.iter_constant_matches():
            print(i)

    def iter_constant_matches(self):
        """

        Returns:

        """
        return self.rows

    def get_input_record(self,id):
        """

        Args:
            id:

        Returns:

        """
        cur=self.conn.cursor( cursor_factory=psycopg2.extras.NamedTupleCursor)
        cur.itersize = 50000
        cur.execute("""SELECT id,input_smiles,normalized_smiles,num_normalized_heavies from idrecord where id=%s;""",(id,))
        rows = cur.fetchall()
        cur.close()
        return(rows[0])


class IndexIteration:
    """
    """
    def __init__(self):
        """

        """
        print("Open IX iteration")
        self.conn = psycopg2.connect("dbname='" + db + "' user='" + user + "' host='" + host + "' password='" + password + "'")
        self.cur = self.conn.cursor(name="index_agg_cursor", cursor_factory=psycopg2.extras.NamedTupleCursor)
        self.cur.itersize = 10000
        self.cur.execute("""SELECT num_cuts,constant_smiles,constant_symmetry_class,variable_part from index_agg limit 100000;""")
        #self.rows = self.cur.fetchall()


    def close(self):
        """

        """
        print("Closing IX iteration")
        self.conn.close()


def test_ix():
    """

    """
    it=IndexIteration()
    for num_cuts,constant_smiles,csc,vp in it.rows:
        print(constant_smiles)
        print(vp)
        for offset1, m1 in enumerate(vp):
            print("\t" + str(offset1))
            print("\t" + str(m1))
            print("\t"+ m1['id'])
    it.close()

def test_idrecord():
    """

    """
    it=IdToRecordIteration()
    for id,input_smiles,normalized_smiles,num_normalized_heavies in it.rows:
        print(id)
        print(input_smiles)
        print(normalized_smiles)
        print(num_normalized_heavies)
    it.close()


def test_db():
    """

    """
    conn = psycopg2.connect("dbname='mmpdb_all' user='" + user + "' host='" + host + "' password='" + password + "'")
    cur = conn.cursor(name="index_agg_cursor", cursor_factory=psycopg2.extras.NamedTupleCursor)
    cur.itersize = 100
    cur.execute("""SELECT * from index_agg where jsonb_array_length(variable_part)>1 limit 10 ;""")

    i = 0
    for r in cur:
        i = i + 1
        print(i)
        print(r.variable_part)
        print(type(r.variable_part))
        print(len(r.variable_part))
        for offset1,m1 in enumerate(r.variable_part):
            for m2 in r.variable_part[offset1 +1:]:
                print(offset1)
                print(m1['id'])
                print(m2['id'])
                print("--------------")

    conn.close()

def test_sframe():
        """

        Returns:

        """
        import  turicreate
        from turicreate import SFrame
        conn = psycopg2.connect(
            "dbname='mmpdb_all' user='" + user + "' host='" + host + "' password='" + password + "'")
        sf=SFrame.from_sql(conn,"SELECT * from index;")
        conn.close()
        return(sf)


def dump(index, id_to_record):
    """

    Args:
        index:
        id_to_record:
    """
    text_file = open(db_home + "index_mem.txt", "w")
    text_file.write(
        "constant_smiles|constant_symmetry_class|num_cuts)|idr|variable_symmetry_class|variable_smiles|attachment_order|enumeration_label\n")
    for (constant_smiles, constant_symmetry_class, num_cuts), matches in sorted(index.items()):
        # print(constant_smiles, constant_symmetry_class, num_cuts, constant_smiles, constant_symmetry_class, matches)
        # print(type(matches))
        # print(len(matches))
        for idr, variable_symmetry_class, variable_smiles, attachment_order, enumeration_label in matches:
            text_file.write(constant_smiles + "|" + constant_symmetry_class + "|" + str(
                num_cuts) + "|" + idr + "|" + variable_symmetry_class + "|" + variable_smiles + "|" + attachment_order + "|" + enumeration_label + "\n")
    text_file.close()

    text_file = open(db_home + "index_id_mem.txt", "w")
    text_file.write(
        "id|input_smiles|num_normalized_heavies,normalized_smiles\n")
    for idr, record in id_to_record.items():
        text_file.write(record.id + "|" + record.input_smiles + "|" + str(
            record.num_normalized_heavies) + "|" + record.normalized_smiles + "\n")
    text_file.close()
