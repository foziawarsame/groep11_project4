"""

Made by: Hans Baker, Fozia Warsame and Joost Kamperman
Bin-1d
Made 10-06-2021

"""

import time
from textwrap import wrap

from flask import Flask, render_template, request
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd

# from tkinter import Tk
# from tkinter.filedialog import askopenfilename, askdirectory

import mysql.connector

import os.path

import matplotlib.pyplot as plt

app = Flask(__name__)


@app.route('/')
def home_page():
    """ Sends you to the home page

    :return: All the data from the database and HTML file for the
    home page
    """

    # These commended functions were used to turn the excel file into
    # blast results, turn them into different files and then read those
    # files into a database.

    # excelfilter()
    # hits()
    return render_template("homepage.html", data=get_sql_data())


@app.route('/search', methods=["POST", "GET"])
def search_page():
    """Sends you to the search page

    input: (potentially) the filled in input from the search page

    :return: All the data from the database unless you filled in the
    input from the search page. Then you will get the data that align
    with the given input. You will also get the HTML file for the
    search page
    """

    # This if statement is made to check when the page is asked for
    # if you did so with a POST command. If it's POST then it will
    # assume you are searching for something and use a different
    # function to find data that align with the given input. If with
    # something else it will then just show all the data.

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
    """Sends you to the about page

    :return: All the data from the database and HTML file for the
    about page
    """

    return render_template("aboutpage.html", data=get_sql_data())


@app.route('/graphs')
def graph_page():
    """Sends you to the graph page

    :return: HTML file for the graph page
    """

    # The two functions were made to turn the data from the database
    # two graphs.

    make_graph_eiwit()
    make_graph_familie()
    return render_template("grafiekpage.html")


def excelfilter():
    """Selects an Excel file that holds data to perform a Blast.
    Splits the data into a dictionary so if can then be blasted
    separately. A dictionary holding the Read, Sequence and Score of
    both the forward and reverse

    :return: header_1 - str - header of the forward sequence
             header_2 - str - header of the reverse sequence
             read_1 - str - forward sequence
             read_2 - str - reverse sequence
    """

    # Open a window to select the proper directory for the file.
    try:
        Tk().withdraw()
        filename = askopenfilename()

        # Reads and turn the Excel File into a dictionary.
        excelread = pd.read_excel(filename, header=None,
                                  names=["ForwardRead",
                                         "ForwardSequence",
                                         "ScoreMatrixForward",
                                         "ReverseRead",
                                         "ReverseSequence",
                                         "ScoreMatrixReverse"])

        # Turn excel data in a dictionary.
        exceldata = excelread.to_dict(orient="records")

        seq = []

        # This for loop send every header and sequence separately to
        # a function that will put them into the data and append the
        # sequences into a list for further use.

        for i in exceldata:
            read_1 = i.get("ForwardSequence")
            read_2 = i.get("ReverseSequence")
            seq.append(read_1)
            seq.append(read_2)
            header_1 = i.get("ForwardRead")
            header_2 = i.get("ReverseRead")

            send_to_sql_fragment(header_1, read_1)
            send_to_sql_fragment(header_2, read_2)

        # A for loop to send the both the sequence and count to the
        # blast function.

        count = 0
        for x in seq:
            count += 1
            blast(x, count)
    except AssertionError:
        print("No file was selected")


def blast(x, count):
    """Blast the given sequence and transfer the result into a xml file
    for further use.

    :param x - str - the sequence
    :param count - int - count of the read
    :return: a xml file containing the result of the blast
    """

    print("Blasting...", count)
    # Parameters to execute the BLAST
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
    """Transfer all xml files in the selected folder into a database.
    By sending the results to other functions

    :return: ac - str - the accessioncodes
            evalue - float - the e-value
            perc_iden - int - the percentage identities
            score - int - the BLAST score
            name - str - the protein name
            i - int - ids for in the database
    """
    try:
        Tk().withdraw()
        path = askdirectory()

        # Loops 200 times with the selected directory to take apart
        # all valuable data from the files and then transfer that
        # data to the database.
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
                                        protein_function(ac),
                                        taxonomy(ac),
                                        (i + 1))

                else:
                    # If the given file does not contain any hits it
                    # will then create an empty hit with empty values
                    # for the other tables.

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
    except FileNotFoundError:
        print("No folder was selected")


def send_to_sql_hit(ac, evalue, score, perc_iden, name, f_id, m_id,
                    fr_id):
    """Transfer the given hit data into the hit table

    :param ac - str - the accessioncodes
           evalue - float - the e-value
           perc_iden - int - the percentage identities
           score - int - the BLAST score
           name - str - the protein name
           f_id - int - function id for in the database
           m_id - int - micro-organism id for in the database
           fr_id - int - fragment id for the database
    :return: A hit send into the hit table of the database
    """

    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11")

    cursor = conn.cursor()

    print(f_id, " ", m_id, " ", fr_id)

    # When a hit is empty impossible values are given to indicate to
    # the program that this read has no hits. So an if statement is used
    # to keep the foreign keys empty so we also don't need to make
    # impossible data to every table.

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
    """Transfer the given header and sequence data into the fragment
    table

    :param headerdata - str - header of the sequence
           sequencedata - str - the DNA sequence
    :return: A read send into the fragment table of the database
    """

    print("Fragment Sending...")
    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11"
                                   )
    cursor = conn.cursor()
    # Sends the data to the database
    cursor.execute("INSERT INTO fragmenten(header, sequentie) "
                   "VALUES('" + headerdata + "',"
                                             " '" + sequencedata + "')")

    conn.commit()
    cursor.close()
    conn.close()
    print("Fragment Send...")


def taxonomy(accession_codes):
    """Searches for the taxonomy of the organism that codes for the
    found protein, using GenBank

    :param accession_codes - str - the accession code from the BLAST
    hit/protein
    :return: family - str - the family name of the found
    organism genus - str - the genus name of the found organism
    species - str - the species name of the found organism
    """

    Entrez.email = "Your.Name.Here@example.org"

    search = Entrez.efetch(db="protein", id=accession_codes,
                           rettype="gb",
                           retmode="text")

    genbank_search = SeqIO.parse(search, "genbank")
    for genbank_entry in genbank_search:
        tax = genbank_entry.annotations["taxonomy"]

        # Gets the right names from the taxonomy. If the name can't
        # be found None will be put in the database
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
    """Transfer the given family, genus and species data into
    the organism table

    :param family - str - the family name of the found organism
           genus - str - the genus name of the found organism
           species - str - the species name of the found organism
    :return: A organism send into the organism table of the database
    """

    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11")

    cursor = conn.cursor()

    cursor.execute("SELECT * FROM micro_organisme")

    # A for loop is made to check if the given organism already
    # exist in the database. If so then it will print that the organism
    # exist and only return the id of the organism. This prevents
    # excessive data.

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
    """Searches for the protein function of the found protein with
    GenBank

    :param accession_codes - str - the accession code from the BLAST
    hit/protein
    :return: function - str - the protein function
    """
    Entrez.email = "Your.Name.Here@example.org"

    search = Entrez.efetch(db="protein", id=accession_codes,
                           rettype="gb",
                           retmode="text")

    genbank_search = SeqIO.parse(search, "genbank")
    for genbank_entry in genbank_search:
        function = genbank_entry.description
        # Gets function from GenBank and replaces a symbol so it can
        # be put in the database
        if "'" in function:
            function = function.replace("'", "")

        return send_to_sql_functie(function)

    search.close()


def send_to_sql_functie(function):
    """Transfer the given function data into the function table

    :param function - str - the protein function
    :return: A function send into the function table of the database
    """

    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11")

    cursor = conn.cursor()

    cursor.execute("SELECT * FROM functie")

    # A for loop is made to check if the given function already
    # exist in the database. If so then it will print that the function
    # exist and only return the id of the function. This prevents
    # excessive data.

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
    """To fetch all data from every table of the database

    :return: data - dict - all hit results
            fragments - dict - the fragment data, reads and sequences
    """

    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11")

    cursor = conn.cursor()
    # Gets data form the database with a MySQL query
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
    """To fetch all data from the selected table of the database that
    align with the given search word

    :param data - str - input from the website
           selection - str - selection on filter from the website
    :return: rows - dict - selected data from database
            fragments - dict - fragments of the selected data
    """

    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11")
    cursor = conn.cursor()

    # The if statement is to check which table was selected on the
    # search page so it only grabs the hits which data align.

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


def make_graph_eiwit():
    """Turn the five most used functions into a graph for the graph page

    :return: A graph png for the graph page
    """

    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11")

    cursor = conn.cursor()

    cursor.execute(
        " select count(replace(eiwit_naam, substr(eiwit_naam, "
        " instr(eiwit_naam, '['), instr(eiwit_naam, ']')), '')), "
        " replace(eiwit_naam, "
        " substr(eiwit_naam, instr(eiwit_naam, '['), "
        " instr(eiwit_naam, ']')), '') "
        " from hit "
        " group by replace(eiwit_naam, substr(eiwit_naam, "
        " instr(eiwit_naam, '['), instr(eiwit_naam, ']')), '') "
        " order by count(replace(eiwit_naam, substr(eiwit_naam, "
        " instr(eiwit_naam, '['), instr(eiwit_naam, ']')), '')) desc "
        " limit 5 ")
    data_eiwit = cursor.fetchall()

    left = [1, 2, 3, 4, 5]

    height = []
    tick_label = []

    for e in data_eiwit:
        tick_label.append(e[1])
        height.append(e[0])

    tick_label = ['\n'.join(wrap(l, 12)) for l in tick_label]
    # Makes the graph
    plt.bar(left, height, tick_label=tick_label,
            width=0.6, color=['red', 'green'])

    plt.xlabel('eiwit naam')
    plt.ylabel('hoeveelheid hits met deze eiwit')
    plt.title('Meest gevonden eiwitten')
    plt.savefig("static/assets/test_image.png", bbox_inches="tight")
    conn.commit()
    cursor.close()
    conn.close()


def make_graph_familie():
    """Turn the five most used families into a graph for the graph page

    :return: A graph png for the graph page
    """

    conn = mysql.connector.connect(host="145.74.104.145",
                                   user="uxhap",
                                   password="Aa635152!",
                                   db="groepje_11")

    cursor = conn.cursor()
    cursor.execute(
        " select count(f3.familie), f3.familie "
        " from hit "
        " left join micro_organisme "
        " f3 on f3.id = hit.micro_organisme_id "
        " group by f3.familie "
        " order by count(f3.familie) desc "
        " limit 5")

    data_familie = cursor.fetchall()

    left = [1, 2, 3, 4, 5]

    height = []
    tick_label = []

    for e in data_familie:
        tick_label.append(e[1])
        height.append(e[0])

    tick_label = ['\n'.join(wrap(l, 11)) for l in tick_label]
    # Makes the graph
    plt.bar(left, height, tick_label=tick_label,
            width=0.6, color=['red', 'green'])

    plt.xlabel('familie naam')
    plt.ylabel('hoeveelheid hits met deze familie')
    plt.title('Meest gevonden familie')
    plt.savefig("static/assets/test_image2.png", bbox_inches="tight")
    conn.commit()
    cursor.close()
    conn.close()


if __name__ == '__main__':
    app.run()
