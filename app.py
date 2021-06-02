import time

from flask import Flask, render_template, request
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd

from tkinter import Tk
from tkinter.filedialog import askopenfilename, askdirectory

import mysql.connector

import openpyxl

import os.path

import glob

app = Flask(__name__)


@app.route('/')
def home_page():
    # excelfilter()
    # hits()
    return render_template("homepage.html", data=get_sql_data())


@app.route('/search', methods=["POST", "GET"])
def search_page():
    if request.method == "POST":
        return render_template("searchpage.html", data=search_database(
            request.form["input"], request.form.get("select"))[0],
                               fragment=search_database(
                                   request.form["input"],
                                   request.form.get("select"))[1])
    else:
        return render_template("searchpage.html", data=get_sql_data()[0]
                               , fragment=get_sql_data()[1])


@app.route('/about')
def about_page():
    return render_template("aboutpage.html", data=get_sql_data())


def excelfilter():
    Tk().withdraw()
    filename = askopenfilename()

    excelread = pd.read_excel(filename, header=None,
                              names=["ForwardRead", "ForwardSequence",
                                     "ScoreMatrixForward",
                                     "ReverseRead", "ReverseSequence",
                                     "ScoreMatrixReverse"])

    exceldata = excelread.to_dict(orient="records")

    seq = []

    for i in exceldata:
        read_1 = i.get("ForwardSequence")
        read_2 = i.get("ReverseSequence")
        seq.append(read_1)
        seq.append(read_2)
        header_1 = i.get("ForwardRead")
        header_2 = i.get("ReverseRead")

        send_to_sql_fragment(header_1, read_1)
        send_to_sql_fragment(header_2, read_2)

    # count = 0
    # for x in seq:
    #    count += 1
    #    if count >= 156:
    #        blast(x, count)


def blast(x, count):
    print("Blasting...", count)
    result_handle = NCBIWWW.qblast("blastx", "nr", x,
                                   matrix_name="BLOSUM62",
                                   expect=1e-07,
                                   hitlist_size=5)

    time.sleep(5)

    file_name = "blast" + str(count) + ".xml"

    blast_results = "C:/Users/fozia/OneDrive - HAN/" \
                    "Tutor Course 4/Programma/blastresults"

    map = os.path.join(blast_results, file_name)

    with open(map, "w") as out_handle:
        out_handle.write(result_handle.read())
        out_handle.close()

    print("Blasted...........")


def hits():
    Tk().withdraw()
    path = askdirectory()

    for i in range(0, 200):
        result_path = path + r"\blast" + str(i) + ".xml"
        with open(result_path) as out_handle:
            print("****newfile****")
            blast_records = NCBIXML.parse(out_handle)
            blast_record = next(blast_records)
            if len(blast_record.alignments) > 0:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        print("***Alignment***")
                        ac = alignment.accession
                        evalue = hsp.expect
                        perc_iden = hsp.identities
                        score = hsp.score

                        title = alignment.title
                        x = title.rfind("| ")
                        name = title[x + 2:]

                        if "'" in name:
                            name = name.replace("'", "")
                    send_to_sql_hit(ac, evalue, score,
                                    perc_iden, name,
                                    protein_function(ac), taxonomy(ac),
                                    (i + 1))

            else:
                ac = "None"
                evalue = -1
                perc_iden = -1
                score = -1
                name = "None"

                send_to_sql_hit(ac, evalue, score, perc_iden,
                                name, 0, 0, (i + 1))

                family = "None"
                genus = "None"
                species = "None"
                send_to_sql_micro_organismen(family, genus,
                                             species)

                function = "None"
                send_to_sql_functie(function)


def send_to_sql_hit(ac, evalue, score, perc_iden, name, f_id, m_id,
                    fr_id):
    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11")

    cursor = conn.cursor()

    print(f_id, " ", m_id, " ", fr_id)

    if f_id == 0:
        cursor.execute(
            "INSERT INTO hit(accesie_code, e_value, score, identity, "
            "eiwit_naam, fragmenten_id) "
            "VALUES('" + ac + "', '" + str(evalue) + "', '" + str(
                score) + "', '" + str(
                perc_iden) + "', '" + name + "'," + str(fr_id) + ")")

        conn.commit()
        cursor.close()
        conn.close()
    else:
        cursor.execute(
            "INSERT INTO hit(accesie_code, e_value, score, identity, "
            "eiwit_naam, functie_id, micro_organisme_id, fragmenten_id)"
            "VALUES('" + ac + "', '" + str(evalue) + "', '" + str(
                score) + "', '" + str(
                perc_iden) + "', '" + name + "', " + str(
                f_id) + "," + str(m_id) + "," + str(fr_id) + ")")

        conn.commit()
        cursor.close()
        conn.close()


def send_to_sql_fragment(headerdata, sequencedata):
    print("Fragment Sending...")
    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11"
                                   )
    cursor = conn.cursor()

    cursor.execute("INSERT INTO fragmenten(header, sequentie) "
                   "VALUES('" + headerdata + "',"
                                             " '" + sequencedata + "')")

    conn.commit()
    cursor.close()
    conn.close()
    print("Fragment Send...")


def taxonomy(accession_codes):
    Entrez.email = "Your.Name.Here@example.org"

    search = Entrez.efetch(db="protein", id=accession_codes,
                           rettype="gb",
                           retmode="text")

    genbank_search = SeqIO.parse(search, "genbank")
    for genbank_entry in genbank_search:
        tax = genbank_entry.annotations["taxonomy"]

        try:
            family = tax[4]
        except IndexError:
            family = "None"
        try:
            genus = tax[5]
        except IndexError:
            genus = "None"
        try:
            species = tax[6]
        except IndexError:
            species = "None"

    search.close()
    return send_to_sql_micro_organismen(family, genus, species)


def send_to_sql_micro_organismen(family, genus, species):
    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11")

    cursor = conn.cursor()

    cursor.execute("SELECT * FROM micro_organisme")

    for e in cursor.fetchall():
        if e[0] == family and e[1] == genus and e[2] == species:
            print("It already exist: ")
            print(e[0] + " " + family)
            print(e[1] + " " + genus)
            print(e[2] + " " + species)
            cursor.close()
            conn.close()
            return e[3]

    cursor.execute(
        "INSERT INTO micro_organisme(familie, geslacht, soort)"
        "VALUES('" + family + "', '" + genus + "', '" + species + "')")

    conn.commit()

    print("New Send Organism found!: ")
    print(family)
    print(genus)
    print(species)
    cursor.execute("SELECT MAX(id) FROM micro_organisme")
    newid = cursor.fetchone()

    cursor.close()
    conn.close()

    return newid[0]


def protein_function(accession_codes):
    Entrez.email = "Your.Name.Here@example.org"

    search = Entrez.efetch(db="protein", id=accession_codes,
                           rettype="gb",
                           retmode="text")

    genbank_search = SeqIO.parse(search, "genbank")
    for genbank_entry in genbank_search:
        function = genbank_entry.description

        if "'" in function:
            function = function.replace("'", "")

        return send_to_sql_functie(function)

    search.close()


def send_to_sql_functie(function):
    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11")

    cursor = conn.cursor()

    cursor.execute("SELECT * FROM functie")

    for e in cursor.fetchall():
        if e[0] == function:
            print("It already exist: ")
            print(e[0] + " " + function)
            cursor.close()
            conn.close()
            return e[1]

    cursor.execute("INSERT INTO functie(eiwit_functie) "
                   "VALUES('" + function + "')")

    conn.commit()

    print("New Send Organism found!: ")
    print(function)
    cursor.execute("SELECT MAX(id) FROM functie")
    newid = cursor.fetchone()

    cursor.close()
    conn.close()
    return newid[0]


def get_sql_data():
    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11")

    cursor = conn.cursor()

    cursor.execute(
        "select * from hit"
        " left join fragmenten f on f.id = hit.fragmenten_id"
        " left join functie f2 on f2.id = hit.functie_id"
        " left join micro_organisme f3"
        " on f3.id = hit.micro_organisme_id")
    data = cursor.fetchall()

    cursor.execute("select * from fragmenten")
    fragments = cursor.fetchall()

    conn.commit()
    cursor.close()
    conn.close()
    return data, fragments


def search_database(data, selection):
    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11")
    cursor = conn.cursor()

    if selection == "Fragment":
        cursor.execute(
            " select * from hit"
            " left join fragmenten f on f.id = hit.fragmenten_id"
            " left join functie f2 on f2.id = hit.functie_id"
            " left join micro_organisme f3"
            " on f3.id = hit.micro_organisme_id"
            " where f.header like '%" + data + "%'")
        rows = cursor.fetchall()

        cursor.execute(
            " select distinct f.id, f.header from hit "
            " left join fragmenten f on f.id = hit.fragmenten_id"
            " where f.header like'%" + data + "%'")
        fragments = cursor.fetchall()
    elif selection == "Hit":
        cursor.execute(
            " select * from hit"
            " left join fragmenten f on f.id = hit.fragmenten_id"
            " left join functie f2 on f2.id = hit.functie_id"
            " left join micro_organisme f3"
            " on f3.id = hit.micro_organisme_id"
            " where eiwit_naam like '%" + data + "%'")
        rows = cursor.fetchall()

        cursor.execute(
            " select distinct f.id, f.header from hit "
            " left join fragmenten f on f.id = hit.fragmenten_id"
            " where eiwit_naam like'%" + data + "%'")
        fragments = cursor.fetchall()
    elif selection == "Function":
        cursor.execute(
            " select * from hit"
            " left join fragmenten f on f.id = hit.fragmenten_id"
            " left join functie f2 on f2.id = hit.functie_id"
            " left join micro_organisme f3"
            " on f3.id = hit.micro_organisme_id"
            " where f2.eiwit_functie like '%" + data + "%'")
        rows = cursor.fetchall()

        cursor.execute(
            " select distinct f.id, f.header from hit "
            " left join fragmenten f on f.id = hit.fragmenten_id"
            " left join functie f2 on f2.id = hit.functie_id"
            " where f2.eiwit_functie like'%" + data + "%'")
        fragments = cursor.fetchall()
    elif selection == "Micro-Organism":
        cursor.execute(
            " select * from hit"
            " left join fragmenten f on f.id = hit.fragmenten_id"
            " left join functie f2 on f2.id = hit.functie_id"
            " left join micro_organisme f3"
            " on f3.id = hit.micro_organisme_id"
            " where f3.familie like '%" + data + "%'")
        rows = cursor.fetchall()

        cursor.execute(
            " select distinct f.id, f.header from hit "
            " left join fragmenten f on f.id = hit.fragmenten_id"
            " left join micro_organisme f3"
            " on f3.id = hit.micro_organisme_id"
            " where f3.familie like'%" + data + "%'")
        fragments = cursor.fetchall()

    conn.commit()
    cursor.close()
    conn.close()
    return rows, fragments


if __name__ == '__main__':
    app.run()
