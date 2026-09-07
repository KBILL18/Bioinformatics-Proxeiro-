"""
Ερώτημα 3 -- Βήμα 2: ανάγνωση των αποτελεσμάτων

Διαβάζει τα XML που παρήγαγε το q3_blast.py και τυπώνει τα στοιχεία που
χρειάζονται για την αναφορά: όνομα πρωτεΐνης, οργανισμό, E-value, score
και ποσοστό ταυτότητας.

Τρέχει σε δευτερόλεπτα -- ξανατρέξ' το όσες φορές θέλεις.
"""

import os
import re
from Bio.Blast import NCBIXML

POSA_HITS = 10   # πόσα αποτελέσματα να τυπωθούν


def organismos(titlos):
    """Βγάζει το όνομα του οργανισμού από τον τίτλο του hit.

    Οι τίτλοι του NCBI έχουν τη μορφή:
        'Chain A, Insulin [Homo sapiens]'
    δηλαδή ο οργανισμός βρίσκεται μέσα στις τελευταίες αγκύλες.
    """
    vrethenta = re.findall(r"\[([^\[\]]+)\]", titlos)
    return vrethenta[-1] if vrethenta else "(άγνωστος)"


def anafora(onoma_arxeiou, perigrafi):
    print("=" * 78)
    print(perigrafi)
    print("=" * 78)

    if not os.path.exists(onoma_arxeiou):
        print(f"[!] Δεν βρέθηκε το {onoma_arxeiou}. Τρέξε πρώτα το q3_blast.py.\n")
        return

    with open(onoma_arxeiou) as f:
        egrafi = NCBIXML.read(f)

    if not egrafi.alignments:
        print("Καμία εμφάνιση. Δοκίμασε πιο χαλαρές παραμέτρους.\n")
        return

    print(f"Μήκος ερωτήματος: {egrafi.query_length} αμινοξέα")
    print(f"Βάση δεδομένων:   {egrafi.database}")
    print(f"Συνολικά hits:    {len(egrafi.alignments)}\n")

    for i, alignment in enumerate(egrafi.alignments[:POSA_HITS], start=1):
        hsp = alignment.hsps[0]          # το καλύτερο τμήμα στοίχισης
        tautotita = 100.0 * hsp.identities / hsp.align_length

        print(f"--- Hit {i} " + "-" * 64)
        print(f"  Τίτλος     : {alignment.title[:110]}")
        print(f"  Οργανισμός : {organismos(alignment.title)}")
        print(f"  Accession  : {alignment.accession}")
        print(f"  E-value    : {hsp.expect:.3g}")
        print(f"  Score      : {hsp.score}  (bits: {hsp.bits:.1f})")
        print(f"  Ταυτότητα  : {hsp.identities}/{hsp.align_length}  ({tautotita:.1f}%)")
        print()
        print(f"    Query {hsp.query_start:>4}  {hsp.query}")
        print(f"              {hsp.match}")
        print(f"    Sbjct {hsp.sbjct_start:>4}  {hsp.sbjct}")
        print()

    # Σύνοψη ειδών -- χρήσιμη για τον σχολιασμό
    eidi = []
    for alignment in egrafi.alignments[:POSA_HITS]:
        o = organismos(alignment.title)
        if o not in eidi:
            eidi.append(o)
    print(f"Διακριτοί οργανισμοί στα πρώτα {POSA_HITS}: {', '.join(eidi)}\n")


if __name__ == "__main__":
    anafora("q3_result.xml", "(α) ΕΛΕΥΘΕΡΗ ΑΝΑΖΗΤΗΣΗ -- ταυτοποίηση της πρωτεΐνης")
    anafora("q3_bacteria.xml", "(β) ΠΕΡΙΟΡΙΣΜΕΝΗ ΣΕ BACTERIA")
