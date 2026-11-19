
import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from io import StringIO


# ==========================================
# 1. PAGE CONFIGURATION
# ==========================================

st.set_page_config(
    page_title="BLAST Automation Pipeline",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 BLAST Automation Pipeline")

st.write(
    "An automated bioinformatics application for "
    "DNA and protein sequence similarity analysis "
    "using Python, Biopython and NCBI BLAST."
)

st.divider()


# ==========================================
# 2. SIDEBAR
# ==========================================

st.sidebar.title("About the Project")

st.sidebar.write(
    "This application automates biological sequence "
    "similarity searching using NCBI BLAST."
)

st.sidebar.info(
    "DNA sequences → BLASTN\n\n"
    "Protein sequences → BLASTP"
)

st.sidebar.markdown("---")

st.sidebar.write("**Technologies Used**")

st.sidebar.write(
    "- Python\n"
    "- Biopython\n"
    "- NCBI BLAST\n"
    "- Streamlit\n"
    "- Pandas\n"
)


# ==========================================
# 3. SEQUENCE INPUT
# ==========================================

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
        placeholder="Enter DNA or protein sequence here..."
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
                f"FASTA file loaded successfully: {record.id}"
            )

        except Exception as e:

            st.error(
                f"Unable to read FASTA file: {e}"
            )


# ==========================================
# 4. NCBI SETTINGS
# ==========================================

st.header("2. NCBI Settings")

email = st.text_input(
    "Enter your email address:",
    placeholder="yourname@example.com",
    help="An email address is required to identify your NCBI request."
)


# ==========================================
# 5. SEQUENCE VALIDATION
# ==========================================

def validate_sequence(sequence, sequence_type):

    sequence = "".join(sequence.split()).upper()

    if not sequence:

        return False, "Please enter a biological sequence."

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
            "Please check your DNA or protein sequence."
        )

    return True, sequence


# ==========================================
# 6. BLAST ANALYSIS FUNCTION
# ==========================================

def run_blast(sequence, sequence_type, email):

    # Configure NCBI identification details

    NCBIWWW.email = email
    NCBIWWW.tool = "BLAST_Automation_Pipeline"

    # Select the appropriate BLAST algorithm

    if sequence_type == "DNA":

        program = "blastn"
        database = "nt"

    else:

        program = "blastp"
        database = "nr"

    # Submit sequence to NCBI BLAST

    result_handle = NCBIWWW.qblast(
        program,
        database,
        sequence,
        hitlist_size=10,
        expect=10.0,
        format_type="XML"
    )

    try:

        # Parse returned XML results

        blast_record = NCBIXML.read(result_handle)

    finally:

        result_handle.close()

    results = []

    # Extract alignment information

    for alignment in blast_record.alignments:

        if not alignment.hsps:
            continue

        # Select the highest-scoring alignment

        hsp = max(
            alignment.hsps,
            key=lambda x: x.bits
        )

        # Calculate percentage identity

        percentage_identity = (
            hsp.identities / hsp.align_length
        ) * 100

        # Calculate query coverage

        query_coverage = (
            (hsp.query_end - hsp.query_start + 1)
            / len(sequence)
        ) * 100

        # Store extracted results

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


# ==========================================
# 7. RUN BLAST ANALYSIS
# ==========================================

st.divider()

st.header("3. Run Sequence Analysis")

if st.button(
    "🔍 Run BLAST Analysis",
    type="primary",
    use_container_width=True
):

    # Validate the input sequence

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
            f"Length: {len(sequence)}"
        )

        # Submit BLAST request

        try:

            with st.spinner(
                "Connecting to NCBI BLAST. "
                "Please wait while your sequence is analyzed..."
            ):

                results = run_blast(
                    sequence,
                    sequence_type,
                    email
                )

            # Display results

            st.divider()

            st.header("4. BLAST Results")

            if results:

                df = pd.DataFrame(results)

                st.success(
                    f"Successfully retrieved {len(results)} matches!"
                )

                # Display summary metrics

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

                # Display results table

                st.subheader("Sequence Similarity Results")

                st.dataframe(
                    df,
                    use_container_width=True,
                    hide_index=True
                )

                # Download CSV results

                csv = df.to_csv(index=False)

                st.download_button(
                    label="⬇ Download Results as CSV",
                    data=csv,
                    file_name="blast_results.csv",
                    mime="text/csv"
                )

            else:

                st.warning(
                    "No matching sequences were returned by BLAST."
                )

        except Exception as e:

            st.error(
                f"BLAST analysis failed: {e}"
            )

            st.info(
                "Possible causes include an NCBI server issue, "
                "network problem, rejected request or invalid input. "
                "Please wait and try again."
            )


# ==========================================
# 8. FOOTER
# ==========================================

st.divider()

st.caption(
    "BLAST Automation Pipeline | "
    "Developed using Python, Biopython, NCBI BLAST and Streamlit"
)
