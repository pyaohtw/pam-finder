import re
import streamlit as st

st.set_page_config(page_title="gRNA PAM Finder", layout="centered")

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def clean_dna(text: str) -> str:
    """Keep only ACGT, uppercase, strip everything else (whitespace, numbers, FASTA headers)."""
    lines = [ln for ln in text.splitlines() if not ln.strip().startswith(">")]
    return re.sub(r"[^ACGT]", "", "".join(lines).upper())


def find_all(hay: str, needle: str):
    """Return all start indices of needle in hay (including overlaps)."""
    idxs, start = [], 0
    while True:
        i = hay.find(needle, start)
        if i == -1:
            break
        idxs.append(i)
        start = i + 1
    return idxs


def parse_grna_input(text: str):
    """
    Parse pasted gRNA input.
    Default: two columns (ID <tab/space> gRNA) per row.
    Fallback: single-column alternating lines (ID, then gRNA).
    Header rows (id / grna) are skipped.
    Returns list of (id, grna_raw).
    """
    raw_lines = [ln.strip() for ln in text.splitlines() if ln.strip()]

    def is_header(tokens):
        low = [t.lower() for t in tokens]
        return set(low) <= {"id", "grna", "guide", "sequence", "seq", "spacer"}

    # split each line into tokens
    tokenized = [re.split(r"[\t ]+", ln) for ln in raw_lines]
    multi = sum(1 for t in tokenized if len(t) >= 2)

    rows = []
    if multi >= max(1, len(tokenized) // 2):
        # two-column mode: first token = ID, last token = gRNA
        for toks in tokenized:
            if len(toks) < 2 or is_header(toks):
                continue
            rows.append((toks[0], toks[-1]))
    else:
        # single-column alternating mode: pair (ID, gRNA)
        flat = [t[0] for t in tokenized if not is_header(t)]
        for i in range(0, len(flat) - 1, 2):
            rows.append((flat[i], flat[i + 1]))
    return rows


def find_pam(spacer_raw: str, sense: str, antisense: str, pam_len: int = 7):
    """
    Search spacer on both strands. Return (grna_out, pam_out).
    - exactly one hit with full PAM -> (spacer, PAM)
    - not found / not enough downstream bases -> ("N/A", "N/A")
    - more than one hit -> ("N/A", "N/A (multiple matches)")
    """
    spacer = re.sub(r"[^ACGT]", "", spacer_raw.upper())
    if not spacer:
        return "N/A", "N/A"

    hits = []  # list of (strand_string, start_index)
    for strand in (sense, antisense):
        for i in find_all(strand, spacer):
            hits.append((strand, i))

    if len(hits) == 0:
        return "N/A", "N/A"
    if len(hits) > 1:
        return "N/A", "N/A (multiple matches)"

    strand, i = hits[0]
    pam = strand[i + len(spacer): i + len(spacer) + pam_len]
    if len(pam) < pam_len:
        return "N/A", "N/A"  # PAM not found (too close to end)
    return spacer, pam


st.title("gRNA PAM Finder")
st.caption("Finds the 7 nt PAM immediately 3' of each gRNA spacer, on either strand.")

col_left, col_right = st.columns(2)

with col_left:
    dna_text = st.text_area(
        "Target DNA",
        height=260,
        placeholder="Paste target DNA sequence here (FASTA header, spaces, and line breaks are ignored).",
    )

with col_right:
    grna_text = st.text_area(
        "gRNA list (paste from Excel: ID + gRNA columns)",
        height=260,
        placeholder="ID\tgRNA\n1934079\tGCACCCCAGGTCCCCATGCCTC\n1934111\tGGTACTCCTTGTTGTTGCCCTC",
    )

pam_len = st.number_input("PAM length (nt)", min_value=1, max_value=30, value=7, step=1)

if st.button("Find PAMs", type="primary"):
    sense = clean_dna(dna_text)
    rows = parse_grna_input(grna_text)

    if not sense:
        st.error("Please provide a valid target DNA sequence.")
    elif not rows:
        st.error("No gRNA rows detected. Paste an ID column and a gRNA column.")
    else:
        antisense = revcomp(sense)
        out_lines = ["ID\tgRNA\tPAM"]
        for gid, gseq in rows:
            grna_out, pam_out = find_pam(gseq, sense, antisense, int(pam_len))
            out_lines.append(f"{gid}\t{grna_out}\t{pam_out}")

        result = "\n".join(out_lines)
        st.subheader("Result (tab-delimited — copy/paste into Excel)")
        st.code(result, language="text")
        st.download_button("Download as .tsv", result, file_name="pam_results.tsv", mime="text/tab-separated-values")
