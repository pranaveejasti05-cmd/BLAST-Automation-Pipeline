
import streamlit as st
import pandas as pd

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

from io import StringIO
import time


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
    "DNA and protein sequence similarity analysis."
)

st.markdown(
    "Developed using Python, Biopython, NCBI BLAST "
    "and Streamlit."
)

st.divider()


# ==========================================
# 2. SIDEBAR
# ==========================================

st.sidebar.title("About the Project")

st.sidebar.write(
    """
    The BLAST Automation Pipeline simplifies
    biological sequence similarity searching.

    Users can submit DNA or protein sequences,
    perform BLAST searches and view matching
    sequences in an organized format.
    """
)

st.sidebar.info(
    "BLASTN is used for DNA sequences.\n\n"
    "BLASTP is used for protein sequences."
)

st.sidebar.markdown("---")

st.sidebar.subheader("Technologies Used")

st.sidebar.write(
    "- Python\n"
    "- Biopython\n"
    "- NCBI BLAST\n"
    "- Streamlit\n"
    "- Pandas"
)


# ==========================================
# 3. SEQUENCE INPUT
# ==========================================

st.header("1. Sequence Input")

sequence_type = st.radio(
    "Select your sequence type:",
    ["DNA", "Protein"],
    horizontal=True
)

input_method = st.radio(
    "Choose your input method:",
    ["Paste Sequence", "Upload FASTA File"],
    horizontal=True
)

sequence = ""


if input_method == "Paste Sequence":

    sequence = st.text_area(
        "Enter your biological sequence:",
        height=160,
        placeholder=(
            "Example:\n"
            "ATGCGTACGTAGCTAG..."
        )
    )

else:

    uploaded_file = st.file_uploader(
        "Upload your FASTA file",
        type=["fasta", "fa", "fna", "faa", "txt"]
    )

    if uploaded_file is not None:

        try:

            content = uploaded_file.getvalue().decode("utf-8")

            record = SeqIO.read(
                StringIO(content),
                "fasta"
            )

            sequence = str(record.seq)

            st.success(
                f"FASTA file loaded successfully: {record.id}"
            )

        except Exception as error:

            st.error(
                f"Unable to read FASTA file: {error}"
            )


# ==========================================
# 4. NCBI SETTINGS
# ==========================================

st.header("2. NCBI Settings")

email = st.text_input(
    "Enter your email address:",
    placeholder="yourname@example.com",
    help="Enter an email address for NCBI identification."
)


# ==========================================
# 5. SEQUENCE VALIDATION
# ==========================================

def validate_sequence(sequence, sequence_type):

    # Remove spaces and convert to uppercase
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
            "Invalid characters found in the sequence. "
            "Please check your DNA or protein input."
        )

    return True, sequence


# ==========================================
# 6. BLAST ANALYSIS FUNCTION
# ==========================================

def run_blast(sequence, sequence_type, email):

    # Configure NCBI identification
    NCBIWWW.email = email
    NCBIWWW.tool = "BLAST_Automation_Pipeline"

    # Select appropriate BLAST program
    if sequence_type == "DNA":

        program = "blastn"
        database = "nt"

    else:

        program = "blastp"
        database = "nr"

    # Submit the sequence to NCBI BLAST
    result_handle = NCBIWWW.qblast(
        program,
        database,
        sequence,
        hitlist_size=10,
        format_type="XML"
    )

    # Parse returned XML results
    try:

        blast_record = NCBIXML.read(result_handle)

    finally:

        result_handle.close()

    results = []

    # Extract information from each alignment
    for alignment in blast_record.alignments:

        if not alignment.hsps:
            continue

        # Select the HSP with the highest bit score
        hsp = max(
            alignment.hsps,
            key=lambda x: x.bits
        )

        # Calculate percentage identity
        percentage_identity = (
            hsp.identities / hsp.align_length
        ) * 100

        # Calculate query coverage
        aligned_query_length = sum(
            1 for character in hsp.query
            if character != "-"
        )

        query_coverage = (
            aligned_query_length / len(sequence)
        ) * 100

        # Store results
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

    # Validate sequence
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

        # Display selected algorithm
        if sequence_type == "DNA":

            st.info(
                "Selected algorithm: BLASTN | Database: nt"
            )

        else:

            st.info(
                "Selected algorithm: BLASTP | Database: nr"
            )

        # Submit BLAST request
        try:

            with st.spinner(
                "Connecting to NCBI BLAST... "
                "Your sequence is being analyzed. "
                "Please wait."
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
                    f"Analysis completed! "
                    f"{len(results)} matching results found."
                )

                # Summary metrics
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

                # Results table
                st.subheader("Sequence Similarity Results")

                st.dataframe(
                    df,
                    use_container_width=True,
                    hide_index=True
                )

                # CSV download
                csv_data = df.to_csv(
                    index=False
                ).encode("utf-8")

                st.download_button(
                    label="⬇ Download Results as CSV",
                    data=csv_data,
                    file_name="BLAST_Results.csv",
                    mime="text/csv",
                    use_container_width=True
                )

            else:

                st.warning(
                    "No matching sequences were returned "
                    "for this search."
                )

        except Exception as error:

            st.error(
                f"BLAST analysis failed: {error}"
            )

            st.warning(
                "NCBI may be temporarily unavailable, "
                "or the request may have been rejected."
            )

            st.info(
                "Wait a few minutes before trying again. "
                "You can also test the sequence directly "
                "on the official NCBI BLAST website."
            )


# ==========================================
# 8. PROJECT INFORMATION
# ==========================================

st.divider()

st.header("5. About BLAST Parameters")

with st.expander("Understand the results"):

    st.write(
        """
        **Percentage Identity:** Measures the percentage
        of identical residues in an alignment.

        **Query Coverage:** Indicates how much of the
        input sequence is included in the alignment.

        **E-value:** Estimates the number of alignments
        expected by chance. Lower values generally
        indicate more significant matches.

        **Bit Score:** Measures alignment quality.
        Higher scores generally indicate better alignments.

        **Accession Number:** Identifies the matching
        database sequence.
        """
    )


# ==========================================
# 9. FOOTER
# ==========================================

st.divider()

st.caption(
    "BLAST Automation Pipeline | "
    "Computational Biology Project | "
    "Python + Biopython + Streamlit"
)
