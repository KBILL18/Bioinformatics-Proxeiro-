"""
Ερώτημα 3 -- Βήμα 1: αναζήτηση BLASTP

Στέλνει την άγνωστη αλληλουχία στους servers του NCBI και ΑΠΟΘΗΚΕΥΕΙ
το αποτέλεσμα σε αρχείο XML.

ΣΗΜΑΝΤΙΚΟ: τρέξε το ΜΙΑ φορά. Η αναζήτηση αργεί (1-5 λεπτά) γιατί
μπαίνει σε ουρά στους servers του NCBI. Μόλις δημιουργηθούν τα αρχεία
XML, δούλευε με το q3_report.py που διαβάζει από αυτά και είναι ακαριαίο.
"""

import os
from Bio.Blast import NCBIWWW

# Η άγνωστη αλληλουχία που δίνει η εκφώνηση (21 αμινοξέα)
SEQ = "GIVEQCCTSICSLYQLENYCN"


def trekse_blast(onoma_arxeiou, entrez_query=None):
    """Τρέχει BLASTP και σώζει το αποτέλεσμα σε XML.

    entrez_query: προαιρετικό φίλτρο οργανισμού, π.χ. "Bacteria[Organism]"
    """
    if os.path.exists(onoma_arxeiou):
        print(f"[!] Το {onoma_arxeiou} υπάρχει ήδη -- παραλείπεται.")
        print("    Σβήσε το αν θέλεις να ξανατρέξει η αναζήτηση.")
        return

    print(f"[>] Αναζήτηση BLASTP -> {onoma_arxeiou}")
    if entrez_query:
        print(f"    φίλτρο οργανισμού: {entrez_query}")
    print("    Περίμενε, αυτό μπορεί να πάρει αρκετά λεπτά...")

    # Παράμετροι βελτιστοποιημένες για ΜΙΚΡΕΣ αλληλουχίες.
    # Με τις προεπιλογές (word_size=3, BLOSUM62) μια αλληλουχία 21
    # αμινοξέων συχνά δεν επιστρέφει τίποτα, γιατί οι default τιμές
    # είναι σχεδιασμένες για πρωτεΐνες εκατοντάδων καταλοίπων.
    parametroi = dict(
        program="blastp",
        database="nr",
        sequence=SEQ,
        word_size=2,          # μικρότερο "παράθυρο" αρχικής ανίχνευσης
        matrix_name="PAM30",  # κατάλληλος πίνακας για κοντές αλληλουχίες
        gapcosts="9 1",       # τα κόστη που συνοδεύουν τον PAM30
        expect=10,            # χαλαρό κατώφλι E-value
        hitlist_size=50,
    )
    if entrez_query:
        parametroi["entrez_query"] = entrez_query

    apotelesma = NCBIWWW.qblast(**parametroi)

    with open(onoma_arxeiou, "w") as f:
        f.write(apotelesma.read())
    apotelesma.close()

    print(f"[✓] Αποθηκεύτηκε: {onoma_arxeiou}\n")


if __name__ == "__main__":
    print(f"Αλληλουχία: {SEQ}  ({len(SEQ)} αμινοξέα)\n")

    # (α) Ελεύθερη αναζήτηση -- για την ταυτοποίηση της πρωτεΐνης
    trekse_blast("q3_result.xml")

    # (β) Περιορισμένη σε βακτήρια -- για την τελευταία ερώτηση
    trekse_blast("q3_bacteria.xml", entrez_query="Bacteria[Organism]")

    print("Τέλος. Τρέξε τώρα το q3_report.py")
