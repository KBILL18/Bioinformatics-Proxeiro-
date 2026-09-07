"""
Ερώτημα 2 -- Βήμα 1 (Fetching)

Ανακτά την αλληλουχία της ανθρώπινης πρωτεΐνης Cytochrome C από το NCBI
με τη βιβλιοθήκη Entrez της Biopython και την αποθηκεύει σε μορφή FASTA.

Τρέχει σε δευτερόλεπτα.
"""

import os
import warnings
from Bio import Entrez, SeqIO

warnings.filterwarnings("ignore")   # σιωπή στην προειδοποίηση για το email

# Το NCBI συνιστά (δεν επιβάλλει) διεύθυνση επικοινωνίας.
# Αν θέλεις, βάλε την εδώ. Αν το αφήσεις None, δουλεύει κανονικά.
Entrez.email = None

ARXEIO = "human_cytc.fasta"

# Το accession της εκφώνησης είναι το UniProt P99999. Το NCBI το
# αναγνωρίζει, αλλά κρατάμε και το αντίστοιχο RefSeq ως εφεδρεία.
IDS = ["P99999", "NP_061820.1"]


def kateveso_anthropini():
    if os.path.exists(ARXEIO):
        print(f"[!] Το {ARXEIO} υπάρχει ήδη -- παραλείπεται η λήψη.")
        return SeqIO.read(ARXEIO, "fasta")

    for pid in IDS:
        print(f"[>] Δοκιμή με ID: {pid}")
        try:
            handle = Entrez.efetch(
                db="protein", id=pid, rettype="fasta", retmode="text"
            )
            record = SeqIO.read(handle, "fasta")
            handle.close()
        except Exception as e:
            print(f"    απέτυχε ({e.__class__.__name__}) -- δοκιμάζω το επόμενο")
            continue

        SeqIO.write(record, ARXEIO, "fasta")
        print(f"[✓] Αποθηκεύτηκε: {ARXEIO}\n")
        return record

    raise SystemExit("Δεν κατέβηκε καμία αλληλουχία. Έλεγξε τη σύνδεση.")


if __name__ == "__main__":
    record = kateveso_anthropini()

    print("=" * 70)
    print("ΑΝΘΡΩΠΙΝΗ CYTOCHROME C")
    print("=" * 70)
    print(f"ID          : {record.id}")
    print(f"Περιγραφή   : {record.description}")
    print(f"Μήκος       : {len(record.seq)} αμινοξέα")
    print()
    # τυπώνουμε σε γραμμές των 60, όπως στη μορφή FASTA
    seq = str(record.seq)
    for i in range(0, len(seq), 60):
        print(f"  {i+1:>4}  {seq[i:i+60]}")
    print()
    print("Τρέξε τώρα το q2_blast.py")
