"""
Ερώτημα 2 -- Βήματα 3, 4 και 5

Βήμα 3 (Data Collection): διαβάζει τα hits του BLAST, κρατά ένα ανά είδος
                          και ανακτά τις πλήρεις αλληλουχίες μέσω Entrez.
Βήμα 4 (Alignment)      : στοιχίζει καθεμία με την ανθρώπινη, ανά ζεύγη,
                          με τον PairwiseAligner.
Βήμα 5 (Identity)       : υπολογίζει το ποσοστό ταυτότητας θέση προς θέση.

Επιπλέον εντοπίζει τις απολύτως συντηρημένες περιοχές, που ζητά η
τελευταία ερώτηση της εκφώνησης.

Τρέχει σε λίγα δευτερόλεπτα μετά την πρώτη φορά (οι αλληλουχίες
αποθηκεύονται τοπικά).
"""

import os
import re
import warnings
from Bio import Align, Entrez, SeqIO
from Bio.Align import substitution_matrices
from Bio.Blast import NCBIXML

warnings.filterwarnings("ignore")

Entrez.email = None          # όπως και στο q2_fetch.py -- προαιρετικό

XML_IN = "q2_blast.xml"
FASTA_ANTHROPOU = "human_cytc.fasta"
CACHE = "q2_omologes.fasta"          # τοπικό αντίγραφο των αλληλουχιών
CSV_OUT = "q2_taftotites.csv"
ALN_OUT = "q2_stoixiseis.txt"

ANA_TAXI = 2                 # πόσα είδη το πολύ από κάθε τάξη θηλαστικών
MIKOS_MIN, MIKOS_MAX = 100, 110   # αποδεκτό μήκος καταχώρησης

# ------------------------------------------------------------------
# Αντιστοίχιση γένους -> τάξης θηλαστικών.
# Χρειάζεται επειδή το BLAST κατατάσσει κατά ομοιότητα: τα πρώτα
# αποτελέσματα είναι όλα πρωτεύοντα. Για να μελετήσουμε τη διατήρηση
# ΣΕ ΟΛΑ τα θηλαστικά, επιλέγουμε ρητά αντιπροσώπους από κάθε τάξη.
# ------------------------------------------------------------------
TAXEIS = {
    "Primates": ["Pan", "Gorilla", "Pongo", "Macaca", "Papio", "Lemur",
                 "Callithrix", "Otolemur", "Microcebus", "Saimiri"],
    "Rodentia": ["Mus", "Rattus", "Cavia", "Marmota", "Sciurus",
                 "Microtus", "Jaculus", "Heterocephalus"],
    "Lagomorpha": ["Oryctolagus", "Ochotona"],
    "Carnivora": ["Ailuropoda", "Crocuta", "Hyaena", "Prionailurus",
                  "Nyctereutes", "Lynx", "Enhydra"],
    "Cetartiodactyla": ["Bos", "Ovis", "Hippopotamus", "Camelus", "Oryx",
                        "Phacochoerus", "Bubalus", "Tursiops",
                        "Balaenoptera", "Phocoena", "Inia", "Mesoplodon"],
    "Perissodactyla": ["Equus", "Ceratotherium", "Tapirus"],
    "Chiroptera": ["Myotis", "Desmodus", "Artibeus", "Pteropus",
                   "Rhinolophus", "Miniopterus", "Plecotus"],
    "Proboscidea": ["Elephas"],
    "Sirenia": ["Trichechus"],
    "Hyracoidea": ["Heterohyrax"],
    "Tubulidentata": ["Orycteropus"],
    "Macroscelidea": ["Rhynchocyon", "Elephantulus"],
    "Cingulata": ["Dasypus"],
    "Pholidota": ["Manis", "Smutsia"],
    "Eulipotyphla": ["Erinaceus", "Suncus", "Galemys"],
    "Dermoptera": ["Cynocephalus"],
    "Diprotodontia": ["Macropus", "Phascolarctos", "Trichosurus",
                      "Vombatus", "Petaurus"],
    "Dasyuromorphia": ["Dasyurus", "Sarcophilus", "Sminthopsis"],
    "Didelphimorphia": ["Monodelphis"],
}

GENOS_SE_TAXI = {g: t for t, gene in TAXEIS.items() for g in gene}


# ============================================================
# Βοηθητικά
# ============================================================

def organismos(titlos):
    """Εξάγει το όνομα του οργανισμού από τον τίτλο του hit."""
    vrethenta = re.findall(r"\[([^\[\]]+)\]", titlos)
    return vrethenta[-1].strip() if vrethenta else None


def taxi_tou(onoma):
    """Επιστρέφει την τάξη θηλαστικών στην οποία ανήκει το είδος, ή None
    αν το γένος δεν αναγνωρίζεται (π.χ. τεχνητές κατασκευές)."""
    if not onoma:
        return None
    genos = onoma.split()[0]
    return GENOS_SE_TAXI.get(genos)


# ============================================================
# Βήμα 3 -- Data Collection
# ============================================================

def epilogi_hits():
    """Επιλέγει αντιπροσώπους από κάθε τάξη θηλαστικών, εφαρμόζοντας
    και φίλτρο μήκους ώστε να αποκλειστούν οι επεκταμένες ή
    προβλεπόμενες καταχωρήσεις που θα αλλοίωναν τα ποσοστά."""
    with open(XML_IN) as f:
        egrafi = NCBIXML.read(f)

    print(f"Συνολικά hits στο XML: {len(egrafi.alignments)}")

    epilegmena = []
    ana_taxi = {}
    idi_eidame = set()
    aporrifthikan_mikos = 0

    for alignment in egrafi.alignments:
        onoma = organismos(alignment.title)
        taxi = taxi_tou(onoma)
        if taxi is None:
            continue
        if onoma in idi_eidame:
            continue

        # Φίλτρο μήκους: το κυτόχρωμα c είναι ~105 κατάλοιπα. Καταχωρήσεις
        # πολύ μεγαλύτερου μήκους είναι αυτόματες προβλέψεις ή κατασκευές
        # και παράγουν τεχνητά χαμηλά ποσοστά ταυτότητας.
        if not (MIKOS_MIN <= alignment.length <= MIKOS_MAX):
            aporrifthikan_mikos += 1
            continue

        if ana_taxi.get(taxi, 0) >= ANA_TAXI:
            continue

        idi_eidame.add(onoma)
        ana_taxi[taxi] = ana_taxi.get(taxi, 0) + 1
        epilegmena.append((onoma, alignment.accession, taxi))

    print(f"Απορρίφθηκαν λόγω ασυνήθιστου μήκους: {aporrifthikan_mikos}")
    print(f"Επιλέχθηκαν {len(epilegmena)} είδη από {len(ana_taxi)} τάξεις\n")
    return epilegmena


def kateveso_allilouxies(epilegmena):
    """Ανακτά τις πλήρεις αλληλουχίες μέσω Entrez, με τοπική αποθήκευση
    ώστε να μην ξαναζητούνται σε κάθε εκτέλεση."""
    zitoumena = {acc.split(".")[0] for _, acc, _ in epilegmena}
    if os.path.exists(CACHE):
        yparxonta = {r.id.split(".")[0]: str(r.seq)
                     for r in SeqIO.parse(CACHE, "fasta")}
        # χρησιμοποιούμε το τοπικό αντίγραφο μόνο αν καλύπτει όλα όσα θέλουμε
        if zitoumena <= set(yparxonta):
            print(f"[!] Χρήση τοπικού αντιγράφου: {CACHE}")
            return yparxonta
        print(f"[!] Το {CACHE} δεν καλύπτει τη νέα επιλογή -- νέα λήψη.")

    accessions = [acc for _, acc, _ in epilegmena]
    print(f"[>] Ανάκτηση {len(accessions)} αλληλουχιών από το Entrez...")

    # Ένα ενιαίο αίτημα αντί για ένα ανά αλληλουχία -- πολύ ταχύτερο
    # και συμβατό με τα όρια ρυθμού του NCBI.
    handle = Entrez.efetch(
        db="protein", id=",".join(accessions), rettype="fasta", retmode="text"
    )
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    SeqIO.write(records, CACHE, "fasta")
    print(f"[✓] Αποθηκεύτηκαν στο {CACHE}\n")

    return {r.id.split(".")[0]: str(r.seq) for r in records}


# ============================================================
# Βήμα 4 -- Alignment
# ============================================================

def ftiakse_aligner():
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


# ============================================================
# Βήμα 5 -- Identity Calculation
# ============================================================

def ypologise_taftotita(stoixismeni_a, stoixismeni_b):
    """Συγκρίνει τις δύο ΣΤΟΙΧΙΣΜΕΝΕΣ αλληλουχίες χαρακτήρα προς χαρακτήρα
    και επιστρέφει (ταυτίσεις, μήκος στοίχισης, ποσοστό).

    Η σύγκριση γίνεται υποχρεωτικά ΜΕΤΑ τη στοίχιση: αν συγκρίναμε τις
    αρχικές αλληλουχίες, μία μόνο διαφορά μήκους θα μετατόπιζε όλες τις
    επόμενες θέσεις και το αποτέλεσμα θα ήταν λανθασμένο.
    """
    tautiseis = 0
    mikos = len(stoixismeni_a)

    for i in range(mikos):
        if stoixismeni_a[i] == stoixismeni_b[i] and stoixismeni_a[i] != "-":
            tautiseis += 1

    pososto = 100.0 * tautiseis / mikos if mikos else 0.0
    return tautiseis, mikos, pososto


# ============================================================
# Συντηρημένες περιοχές
# ============================================================

def synthrimenes_theseis(anthropini, stoixiseis):
    """Επιστρέφει τις θέσεις της ανθρώπινης αλληλουχίας στις οποίες ΟΛΑ
    τα είδη έχουν το ίδιο αμινοξύ."""
    syntirimeni = [True] * len(anthropini)

    for top, bot in stoixiseis:          # top = άνθρωπος, bot = άλλο είδος
        thesi_anthropou = 0
        for i in range(len(top)):
            if top[i] == "-":
                continue                  # κενό στον άνθρωπο -> δεν αντιστοιχεί θέση
            if bot[i] != top[i]:
                syntirimeni[thesi_anthropou] = False
            thesi_anthropou += 1

    return syntirimeni


def perioxes(syntirimeni, anthropini, elaxisto_mikos=3):
    """Ομαδοποιεί διαδοχικές συντηρημένες θέσεις σε περιοχές."""
    apotelesma = []
    arxi = None
    for i, ok in enumerate(syntirimeni + [False]):
        if ok and arxi is None:
            arxi = i
        elif not ok and arxi is not None:
            if i - arxi >= elaxisto_mikos:
                apotelesma.append((arxi + 1, i, anthropini[arxi:i]))
            arxi = None
    return apotelesma


# ============================================================
# Κύριο πρόγραμμα
# ============================================================

def main():
    if not os.path.exists(XML_IN):
        raise SystemExit(f"Δεν βρέθηκε το {XML_IN}. Τρέξε πρώτα το q2_blast.py.")

    anthropini = str(SeqIO.read(FASTA_ANTHROPOU, "fasta").seq)
    print(f"Ανθρώπινη αλληλουχία: {len(anthropini)} αμινοξέα\n")

    epilegmena = epilogi_hits()
    allilouxies = kateveso_allilouxies(epilegmena)
    aligner = ftiakse_aligner()

    apotelesmata = []
    stoixiseis = []

    with open(ALN_OUT, "w", encoding="utf-8") as f_aln:
        for onoma, acc, taxi in epilegmena:
            seq = allilouxies.get(acc.split(".")[0])
            if seq is None:
                print(f"[!] Λείπει η αλληλουχία για {acc} ({onoma})")
                continue

            aln = aligner.align(anthropini, seq)[0]
            top, bot = aln[0], aln[1]

            tautiseis, mikos, pososto = ypologise_taftotita(top, bot)
            apotelesmata.append((onoma, taxi, acc, len(seq), tautiseis, mikos, pososto))
            stoixiseis.append((top, bot))

            f_aln.write(f"{'='*70}\n{onoma}  [{taxi}]  ({acc})\n")
            f_aln.write(f"Ταυτότητα: {tautiseis}/{mikos} = {pososto:.2f}%\n\n")
            for i in range(0, mikos, 60):
                f_aln.write(f"Human  {top[i:i+60]}\n")
                f_aln.write("       " + "".join(
                    "|" if a == b and a != "-" else " "
                    for a, b in zip(top[i:i+60], bot[i:i+60])) + "\n")
                f_aln.write(f"Other  {bot[i:i+60]}\n\n")

    # --- Πίνακας αποτελεσμάτων -----------------------------------
    apotelesmata.sort(key=lambda x: -x[6])

    print("=" * 78)
    print("ΠΟΣΟΣΤΑ ΤΑΥΤΟΤΗΤΑΣ ΩΣ ΠΡΟΣ ΤΗΝ ΑΝΘΡΩΠΙΝΗ CYTOCHROME C")
    print("=" * 78)
    print(f"{'#':>3}  {'Οργανισμός':<30}{'Τάξη':<18}{'Accession':<14}"
          f"{'Ταυτ.':>9}{'%':>8}")
    print("-" * 84)
    for i, (onoma, taxi, acc, mk, tt, ml, ps) in enumerate(apotelesmata, 1):
        print(f"{i:>3}  {onoma[:29]:<30}{taxi[:17]:<18}{acc:<14}"
              f"{tt:>4}/{ml:<4}{ps:>8.2f}")
    print("-" * 84)

    if apotelesmata:
        posostas = [x[6] for x in apotelesmata]
        print(f"Μέγιστο: {max(posostas):.2f}%   "
              f"Ελάχιστο: {min(posostas):.2f}%   "
              f"Μέσος όρος: {sum(posostas)/len(posostas):.2f}%\n")
        print(f"Μεγαλύτερη απόκλιση: {apotelesmata[-1][0]} "
              f"[{apotelesmata[-1][1]}] ({apotelesmata[-1][6]:.2f}%)\n")

    # --- Συντηρημένες περιοχές -----------------------------------
    syntirimeni = synthrimenes_theseis(anthropini, stoixiseis)
    plithos = sum(syntirimeni)

    print("=" * 78)
    print("ΑΠΟΛΥΤΩΣ ΣΥΝΤΗΡΗΜΕΝΕΣ ΠΕΡΙΟΧΕΣ (100% σε όλα τα είδη)")
    print("=" * 78)
    print(f"Συντηρημένες θέσεις: {plithos}/{len(anthropini)} "
          f"({100*plithos/len(anthropini):.1f}%)\n")
    print("Συνεχόμενες περιοχές μήκους >= 3 καταλοίπων:")
    for arxi, telos, akolouthia in perioxes(syntirimeni, anthropini):
        print(f"  θέσεις {arxi:>3}-{telos:<3}  ({telos-arxi+1:>2} κατάλοιπα)  {akolouthia}")

    # --- Αποθήκευση σε CSV ---------------------------------------
    with open(CSV_OUT, "w", encoding="utf-8") as f:
        f.write("organismos,taxi,accession,mikos,tautiseis,mikos_stoixisis,pososto\n")
        for row in apotelesmata:
            f.write(f'"{row[0]}","{row[1]}",{row[2]},{row[3]},'
                    f'{row[4]},{row[5]},{row[6]:.2f}\n')

    print(f"\n[✓] {CSV_OUT}      -- πίνακας για την αναφορά")
    print(f"[✓] {ALN_OUT}  -- οι πλήρεις στοιχίσεις")


if __name__ == "__main__":
    main()
