
import streamlit as st
import pandas as pd
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from io import StringIO

# ---------------------------------
# PAGE CONFIGURATION
# ---------------------------------

st.set_page_config(
    page_title="BLAST Automation Pipeline",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 BLAST Automation Pipeline")

st.write(
    "An automated bioinformatics tool for performing "
    "DNA and protein sequence similarity analysis."
)

st.divider()

# ---------------------------------
# SIDEBAR
# ---------------------------------

st.sidebar.header("About the Project")

st.sidebar.write(
    "This application uses Python, Biopython and "
    "NCBI BLAST to identify similar biological sequences."
)

st.sidebar.info(
    "DNA sequences use BLASTN. "
    "Protein sequences use BLASTP."
)

# ---------------------------------
# SEQUENCE INPUT
# ---------------------------------

st.header("1. Sequence Input")

sequence_type = st.radio(
    "Select sequence type:",
    ["DNA", "Protein"],
    horizontal=True
)

input_method = st.radio(
    "Choose input method:",
    ["Paste Sequence", "Upload FASTA File"],
    horizontal=True
)

sequence = ""

if input_method == "Paste Sequence":

    sequence = st.text_area(
        "Enter your biological sequence:",
        height=150,
        placeholder="Example: ATGCGTACGTAGCTAG..."
    )

else:

    uploaded_file = st.file_uploader(
        "Upload a FASTA file",
        type=["fasta", "fa", "fna", "faa", "txt"]
    )

    if uploaded_file is not None:

        try:
            content = uploaded_file.getvalue().decode("utf-8")
            record = SeqIO.read(
                StringIO(content),
                "fasta"
            )
            sequence = str(record.seq).upper()

            st.success(
                f"Successfully loaded: {record.id}"
            )

        except Exception as e:
            st.error(
                f"Could not read FASTA file: {e}"
            )

# ---------------------------------
# EMAIL INPUT
# ---------------------------------

st.header("2. NCBI Settings")

email = st.text_input(
    "Enter your email address:",
    placeholder="yourname@example.com",
    help="NCBI requests an email address for identification."
)

# ---------------------------------
# SEQUENCE VALIDATION
# ---------------------------------

def validate_sequence(sequence, sequence_type):

    sequence = "".join(sequence.split()).upper()

    if not sequence:
        return False, "Please enter a sequence."

    if sequence_type == "DNA":

        valid_characters = set(
            "ACGTNRYKMSWBDHV"
        )

    else:

        valid_characters = set(
            "ACDEFGHIKLMNPQRSTVWYBXZJUO"
        )

    if not set(sequence).issubset(valid_characters):

        return False, (
            "Invalid characters detected. "
            "Please check your sequence."
        )

    return True, sequence


# ---------------------------------
# BLAST ANALYSIS FUNCTION
# ---------------------------------

def run_blast(sequence, sequence_type, email):

    Entrez.email = email

    if sequence_type == "DNA":

        program = "blastn"
        database = "nt"

    else:

        program = "blastp"
        database = "nr"

    result_handle = NCBIWWW.qblast(
        program,
        database,
        sequence,
        hitlist_size=10,
        format_type="XML"
    )

    blast_record = NCBIXML.read(result_handle)

    result_handle.close()

    results = []

    for alignment in blast_record.alignments:

        if not alignment.hsps:
            continue

        # Select the highest-scoring HSP
        hsp = max(
            alignment.hsps,
            key=lambda x: x.bits
        )

        percentage_identity = (
            hsp.identities / hsp.align_length
        ) * 100

        query_coverage = (
            (hsp.query_end - hsp.query_start + 1)
            / len(sequence)
        ) * 100

        results.append({

            "Accession": alignment.accession,

            "Description": alignment.title,

            "Identity (%)": round(
                percentage_identity, 2
            ),

            "Query Coverage (%)": round(
                query_coverage, 2
            ),

            "E-value": hsp.expect,

            "Bit Score": round(
                hsp.bits, 2
            ),

            "Alignment Length": hsp.align_length

        })

    return results


# ---------------------------------
# RUN BLAST BUTTON
# ---------------------------------

st.divider()

st.header("3. Run Sequence Analysis")

if st.button(
    "🔍 Run BLAST Analysis",
    type="primary",
    use_container_width=True
):

    is_valid, validation_result = validate_sequence(
        sequence,
        sequence_type
    )

    if not is_valid:

        st.error(validation_result)

    elif not email or "@" not in email:

        st.error(
            "Please enter a valid email address."
        )

    else:

        sequence = validation_result

        st.success(
            f"Sequence validated successfully! "
            f"Length: {len(sequence)} bases/amino acids."
        )

        try:

            with st.spinner(
                "Submitting sequence to NCBI BLAST... "
                "This may take several minutes."
            ):

                results = run_blast(
                    sequence,
                    sequence_type,
                    email
                )

            st.divider()

            st.header("4. BLAST Results")

            if results:

                df = pd.DataFrame(results)

                st.success(
                    f"Found {len(results)} matching results!"
                )

                col1, col2, col3 = st.columns(3)

                col1.metric(
                    "Query Length",
                    len(sequence)
                )

                col2.metric(
                    "Matches Found",
                    len(results)
                )

                col3.metric(
                    "Best Identity",
                    f"{df['Identity (%)'].max()}%"
                )

                st.subheader("Similarity Results")

                st.dataframe(
                    df,
                    use_container_width=True,
                    hide_index=True
                )

                csv = df.to_csv(index=False)

                st.download_button(
                    label="⬇ Download Results as CSV",
                    data=csv,
                    file_name="blast_results.csv",
                    mime="text/csv"
                )

            else:

                st.warning(
                    "No significant sequence matches were returned."
                )

        except Exception as e:

            st.error(
                f"BLAST analysis failed: {e}"
            )

            st.info(
                "Check your internet connection and sequence. "
                "NCBI may also be temporarily unavailable."
            )


# ---------------------------------
# FOOTER
# ---------------------------------

st.divider()

st.caption(
    "BLAST Automation Pipeline | "
    "Developed using Python, Biopython and Streamlit"
)
