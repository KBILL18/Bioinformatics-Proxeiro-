"""
Ερώτημα 2 -- Βήμα 2 (BLAST Query)

Αναζήτηση ομόλογων αλληλουχιών της ανθρώπινης Cytochrome C, περιορισμένη
στην κλάση των θηλαστικών (Mammalia) μέσω του φίλτρου entrez_query.

ΣΗΜΑΝΤΙΚΟ: τρέξε το ΜΙΑ φορά. Είναι το πιο αργό βήμα ολόκληρης της
εργασίας -- υπολόγισε 5 έως 15 λεπτά. Το αποτέλεσμα αποθηκεύεται σε XML,
ώστε η ανάλυση (q2_analyze.py) να μπορεί να επαναλαμβάνεται ελεύθερα.
"""

import os
from Bio import SeqIO
from Bio.Blast import NCBIWWW

FASTA_IN = "human_cytc.fasta"
XML_OUT = "q2_blast.xml"

# Πόσες εμφανίσεις να ζητηθούν. Θέλουμε αρκετές ώστε, μετά την αφαίρεση
# των διπλοεγγραφών ανά είδος, να μείνουν 15-20 διαφορετικά θηλαστικά
# από διάφορες τάξεις.
POSES_EMFANISEIS = 250


def trekse_blast():
    if os.path.exists(XML_OUT):
        print(f"[!] Το {XML_OUT} υπάρχει ήδη -- παραλείπεται η αναζήτηση.")
        print("    Σβήσε το αν θέλεις να ξανατρέξει.")
        return

    if not os.path.exists(FASTA_IN):
        raise SystemExit(f"Δεν βρέθηκε το {FASTA_IN}. Τρέξε πρώτα το q2_fetch.py.")

    record = SeqIO.read(FASTA_IN, "fasta")
    print(f"[>] Ερώτημα: {record.id}  ({len(record.seq)} αμινοξέα)")
    print("    Αλγόριθμος : blastp")
    print("    Βάση       : nr")
    print("    Φίλτρο     : Mammalia[Organism]")
    print("    Περίμενε, αυτό παίρνει συνήθως 5-15 λεπτά...\n")

    apotelesma = NCBIWWW.qblast(
        program="blastp",
        database="nr",
        sequence=str(record.seq),
        entrez_query="Mammalia[Organism]",   # ο περιορισμός στα θηλαστικά
        hitlist_size=POSES_EMFANISEIS,
    )

    with open(XML_OUT, "w") as f:
        f.write(apotelesma.read())
    apotelesma.close()

    print(f"[✓] Αποθηκεύτηκε: {XML_OUT}")


if __name__ == "__main__":
    trekse_blast()
    print("\nΤρέξε τώρα το q2_analyze.py")
