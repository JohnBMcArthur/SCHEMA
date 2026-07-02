"""
Automated sequence retrieval and alignment utilities.
Uses EBI MUSCLE API for web-based alignment (no local installation required).
"""

import tempfile
import os
import requests
import time
import random
from Bio import AlignIO
from xml.etree import ElementTree

from utils.config import BLAST_PARAMS, EBI_API_PARAMS, ALPHAFOLD_PARAMS
from utils.sequence_utils import calculate_identity, extract_uniprot_id, clean_sequence_id
from utils.temp_file_manager import temp_file_manager


def blast_search_sequences(query_sequence, num_sequences=10, min_identity=0.5, max_identity=0.9,
                          progress_callback=None):
    """
    Search for similar sequences using EBI BLAST against AlphaFold database.
    Returns sequences with min_identity <= identity < max_identity, randomly selected
    to ensure diversity (no more than max_identity between any pair).
    
    Args:
        query_sequence: Protein sequence string
        num_sequences: Number of sequences to retrieve (default: 10)
        min_identity: Minimum identity for BLAST hits (default: 0.5 = 50%)
        max_identity: Maximum identity for sequence selection (default: 0.9 = 90%)
        progress_callback: Optional callback for progress updates
        
    Returns:
        tuple: (list of (sequence_id, sequence) tuples, tsv_data_string or None, tsv_error_string or None, best_hit_info)
               - sequences: List including query as first entry
               - tsv_data: Raw TSV data from BLAST (None if XML fallback was used)
               - tsv_error: Error message if TSV retrieval failed (None otherwise)
               - best_hit: Tuple (uniprot_id, identity_pct) for best match (highest identity, regardless of filtering) or None
    """
    if progress_callback:
        progress_callback(0.1, "Submitting BLAST query to EBI...")
    
    try:
        # Use EBI BLAST API (more reliable, provides progress updates)
        url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
        data = {
            "email": "johnbmcarthur@gmail.com",
            "program": "blastp",
            "matrix": "BLOSUM62",
            "alignments": 1000,
            "scores": 1000,
            "exp": "1e-5",  # Must be string, and must be one of: 1e-200, 1e-100, 1e-50, 1e-10, 1e-5, 1e-4, 1e-3...
            "filter": "F",
            "gapalign": "true",
            "compstats": "F",
            "align": 0,  # Don't request alignments in initial result (we'll parse from XML)
            "stype": "protein",
            "sequence": query_sequence,
            "database": "afdb"  # AlphaFold database
        }
        
        response = requests.post(url, data=data, timeout=30)
        if not response.ok or "<" in response.text:
            raise Exception(f"Failed to submit BLAST job: {response.text}")
        
        job_id = response.text.strip()
        
        if progress_callback:
            progress_callback(0.2, f"BLAST job submitted (ID: {job_id[:8]}...), waiting for results...")
        
        # Poll for job completion
        status_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}"
        start_time = time.time()
        max_wait_time = EBI_API_PARAMS['max_wait_time']
        poll_interval = EBI_API_PARAMS['poll_interval']
        
        for i in range(int(max_wait_time / poll_interval)):
            status_response = requests.get(status_url, timeout=30)
            if status_response.status_code != 200:
                raise Exception(f"Failed to check job status: {status_response.status_code}")
            
            status = status_response.text.strip()
            elapsed = time.time() - start_time
            
            if progress_callback:
                elapsed_min = int(elapsed // 60)
                elapsed_sec = int(elapsed % 60)
                progress_callback(0.2 + (i / (max_wait_time / poll_interval)) * 0.5, 
                                f"BLAST status: {status} (elapsed: {elapsed_min}m {elapsed_sec}s)")
            
            if status == "FINISHED":
                break
            elif status in ("RUNNING", "PENDING", "QUEUED"):
                time.sleep(poll_interval)
            else:
                raise Exception(f"EBI BLAST job failed with status: {status}")
        else:
            raise TimeoutError("EBI BLAST job did not finish in time")
        
        if progress_callback:
            progress_callback(0.7, "BLAST completed, retrieving results...")
        
        # Small delay to ensure results are fully available
        time.sleep(2)
        
        # Get BLAST results in TSV format first (more reliable than XML)
        # TSV format: query acc, subject acc, identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
        tsv_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/tsv"
        tsv_response = requests.get(tsv_url, timeout=30)
        
        # Store TSV data for debugging
        tsv_data = None
        tsv_error = None
        
        if tsv_response.status_code == 200 and tsv_response.text.strip():
            tsv_data = tsv_response.text  # Store the TSV data
            
            # Parse TSV results
            if progress_callback:
                progress_callback(0.8, "Parsing BLAST results...")
            
            sequences = []
            seen_ids = set()
            
            # Add query sequence first
            sequences.append(("query", query_sequence))
            seen_ids.add("query")
            
            # Track best hit for structure selection (regardless of filtering)
            best_hit = None
            best_hit_identity = 0.0
            best_hit_uniprot_id = None
            
            # First, collect candidate hits (up to max_candidates to process)
            candidate_hits = []
            lines = tsv_response.text.strip().split('\n')
            min_identity_pct = min_identity * 100
            max_identity_pct = max_identity * 100
            
            # Try to detect TSV format - check first non-comment line
            header_line = None
            for line in lines:
                if line.strip() and not line.startswith('#'):
                    header_line = line
                    break
            
            # Parse based on format
            # Format 1: Standard EBI TSV (query acc, subject acc, identity, ...)
            # Format 2: Alternative format (Hit, DB, Accession, Description, ...)
            use_alt_format = False
            if header_line and 'Hit' in header_line and 'Accession' in header_line:
                use_alt_format = True
                if progress_callback:
                    progress_callback(0.81, f"Detected alternative TSV format (header: {header_line[:50]}...)")
            else:
                if progress_callback:
                    progress_callback(0.81, f"Using standard TSV format")
            
            for line in lines:
                if len(candidate_hits) >= BLAST_PARAMS['max_candidates']:
                    break
                
                if line.startswith('#') or not line.strip() or (header_line and line == header_line):
                    continue
                
                fields = line.split('\t')
                
                if use_alt_format:
                    # Alternative format: Hit, DB, Accession, Description, Organism, Length, Score(Bits), Identities(%), Positives(%), E()
                    # Note: Description may contain tabs, so we parse from the right (last columns are fixed)
                    if len(fields) < 10:
                        continue
                    
                    # Parse from the right since last columns are fixed
                    # E() is last column, Positives(%) is second to last, Identities(%) is third to last
                    try:
                        evalue_str = fields[-1].strip()  # Last column: E()
                        positives_str = fields[-2].rstrip('%').strip()  # Second to last: Positives(%)
                        identity_str = fields[-3].rstrip('%').strip()  # Third to last: Identities(%)
                        score_str = fields[-4].strip()  # Score(Bits)
                        length_str = fields[-5].strip()  # Length
                        
                        # Accession is always column 2 (third column, 0-indexed)
                        subject_acc = fields[2]  # Accession column
                        
                        # Parse identity percentage
                        identity_pct = float(identity_str)
                        
                        # Handle e-value format (might be "0.0" or scientific notation like "4.9e-111")
                        if evalue_str:
                            evalue = float(evalue_str)
                        else:
                            evalue = 1e5
                    except (ValueError, IndexError) as e:
                        # Skip if can't parse - likely malformed line
                        continue
                else:
                    # Standard format: query acc, subject acc, identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
                    if len(fields) < 12:
                        continue
                    subject_acc = fields[1]
                    identity_pct = float(fields[2])
                    evalue = float(fields[10])
                
                # Track best hit for structure selection (only check e-value, ignore identity filtering)
                if evalue <= 1e-5 and identity_pct > best_hit_identity:
                    uniprot_id_temp = extract_uniprot_id(subject_acc)
                    if uniprot_id_temp:
                        best_hit_identity = identity_pct
                        best_hit_uniprot_id = uniprot_id_temp
                        best_hit = (uniprot_id_temp, identity_pct)
                
                # Check e-value and identity range (min_identity <= identity < max_identity)
                # This applies to both formats
                if evalue <= 1e-5 and min_identity_pct <= identity_pct < max_identity_pct:
                    # Extract UniProt ID
                    uniprot_id = extract_uniprot_id(subject_acc)
                    seq_id = clean_sequence_id(uniprot_id)
                    
                    if seq_id not in seen_ids:
                        # Fetch sequence from UniProt
                        try:
                            uniprot_url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
                            uniprot_response = requests.get(uniprot_url, timeout=10)
                            if uniprot_response.status_code == 200:
                                fasta_lines = uniprot_response.text.split('\n')
                                seq_lines = [l for l in fasta_lines[1:] if l.strip()]
                                subject_seq = ''.join(seq_lines).upper()
                                
                                if len(subject_seq) > 0:
                                    candidate_hits.append((seq_id, subject_seq, identity_pct))
                                    seen_ids.add(seq_id)
                        except Exception as e:
                            # Log but continue - don't fail on individual sequence fetch errors
                            if progress_callback:
                                progress_callback(0.85, f"Warning: Could not fetch sequence for {uniprot_id}: {str(e)[:50]}")
                            continue
            
            # Randomly shuffle candidates
            random.shuffle(candidate_hits)
            
            # Select sequences non-redundantly
            max_identity_pct_check = max_identity * 100
            filtered_by_query_identity = []
            filtered_by_pairwise = []
            
            for seq_id, seq, identity_pct in candidate_hits:
                if len(sequences) >= num_sequences + 1:  # +1 for query
                    break
                
                # Verify identity with query is within range
                query_identity = calculate_identity(query_sequence, seq)
                if query_identity >= max_identity_pct_check:
                    filtered_by_query_identity.append((seq_id, query_identity))
                    continue
                
                # Check identity with previously selected sequences
                valid = True
                conflicting_seq = None
                for selected_id, selected_seq in sequences:
                    if selected_id == "query":
                        continue
                    pairwise_identity = calculate_identity(seq, selected_seq)
                    if pairwise_identity >= max_identity_pct_check:
                        valid = False
                        conflicting_seq = selected_id
                        break
                
                if valid:
                    sequences.append((seq_id, seq))
                else:
                    filtered_by_pairwise.append((seq_id, conflicting_seq))
            
            # Store debug info in a way that can be accessed (optional - for future use)
            # For now, the page-level debug will show the main filtering
            
            # Return sequences, TSV data, error info, and best hit for structure selection
            tsv_error_val = None  # TSV succeeded, so no error
            best_hit_val = best_hit if best_hit else None
            if progress_callback:
                progress_callback(1.0, f"Found {len(sequences)} sequence(s)")
            return sequences, tsv_data, tsv_error_val, best_hit_val
        else:
            # Fallback to XML if TSV fails
            result_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/xml"
            result_response = requests.get(result_url, timeout=30)
            
            if result_response.status_code != 200:
                # Try to get error details
                error_msg = result_response.text[:500] if result_response.text else "No error details"
                raise Exception(f"Failed to retrieve BLAST results: {result_response.status_code} - {error_msg}")
            
            # No TSV data available when using XML fallback
            tsv_data = None
            tsv_error = f"TSV retrieval failed: status_code={tsv_response.status_code}"
            if tsv_response.text:
                tsv_error += f", response_preview={tsv_response.text[:200]}"
            else:
                tsv_error += ", empty response"
            
            if progress_callback:
                progress_callback(0.8, "Parsing BLAST results...")
            
            # Track best hit for structure selection (regardless of filtering)
            best_hit = None
            best_hit_identity = 0.0
            best_hit_uniprot_id = None
            
            # Parse XML results - collect candidates first
            candidate_hits = []
            seen_ids = set()
            sequences = []  # Initialize sequences list for XML path
            min_identity_pct = min_identity * 100
            max_identity_pct = max_identity * 100
            
            try:
                root = ElementTree.fromstring(result_response.text)
                ns = {'ebi': 'http://www.ebi.ac.uk/schema'}
                
                for hit in root.findall('.//ebi:hit', ns):
                    if len(candidate_hits) >= BLAST_PARAMS['max_candidates']:
                        break
                    
                    uniprot_id = hit.attrib.get('ac')  # UniProt accession
                    if not uniprot_id:
                        continue
                    
                    align = hit.find('.//ebi:alignment', ns)
                    if align is not None:
                        identity_elem = align.find('ebi:identity', ns)
                        query_len_elem = align.find('ebi:queryLen', ns)
                        
                        if identity_elem is not None and query_len_elem is not None:
                            # Identity is a count, convert to percentage
                            num_identities = int(identity_elem.text)
                            query_len = int(query_len_elem.text)
                            identity_percent = (num_identities / query_len * 100) if query_len > 0 else 0
                            
                            # Track best hit for structure selection (ignore identity filtering)
                            uniprot_id_clean_temp = extract_uniprot_id(uniprot_id)
                            if uniprot_id_clean_temp and identity_percent > best_hit_identity:
                                best_hit_identity = identity_percent
                                best_hit_uniprot_id = uniprot_id_clean_temp
                                best_hit = (uniprot_id_clean_temp, identity_percent)
                            
                            # Check identity range (min_identity <= identity < max_identity)
                            if min_identity_pct <= identity_percent < max_identity_pct:
                                # Extract UniProt ID
                                uniprot_id_clean = extract_uniprot_id(uniprot_id)
                                seq_id = clean_sequence_id(uniprot_id_clean)
                                
                                # Fetch sequence from UniProt
                                try:
                                    if seq_id not in seen_ids:
                                        uniprot_url = f"https://www.uniprot.org/uniprot/{uniprot_id_clean}.fasta"
                                        uniprot_response = requests.get(uniprot_url, timeout=10)
                                        if uniprot_response.status_code == 200:
                                            fasta_lines = uniprot_response.text.split('\n')
                                            seq_lines = [l for l in fasta_lines[1:] if l.strip()]
                                            subject_seq = ''.join(seq_lines).upper()
                                            
                                            if len(subject_seq) > 0:
                                                candidate_hits.append((seq_id, subject_seq, identity_percent))
                                                seen_ids.add(seq_id)
                                except Exception:
                                    continue
                
                # Randomly shuffle candidates
                random.shuffle(candidate_hits)
                
                # Select sequences non-redundantly
                sequences = [("query", query_sequence)]  # Always include query
                max_identity_pct_check = max_identity * 100
                
                for seq_id, seq, identity_percent in candidate_hits:
                    if len(sequences) >= num_sequences + 1:  # +1 for query
                        break
                    
                    # Verify identity with query is within range
                    query_identity = calculate_identity(query_sequence, seq)
                    if query_identity >= max_identity_pct_check:
                        continue
                    
                    # Check identity with previously selected sequences
                    valid = True
                    for selected_id, selected_seq in sequences:
                        if selected_id == "query":
                            continue
                        pairwise_identity = calculate_identity(seq, selected_seq)
                        if pairwise_identity >= max_identity_pct_check:
                            valid = False
                            break
                    
                    if valid:
                        sequences.append((seq_id, seq))
            except Exception as e:
                raise Exception(f"Failed to parse BLAST XML: {str(e)}")
            
            # Ensure we always return at least the query sequence
            if not sequences:
                sequences = [("query", query_sequence)]
            
            if progress_callback:
                progress_callback(1.0, f"Found {len(sequences)} sequence(s) (XML fallback)")
            
            # Return sequences, TSV data, error info, and best hit for structure selection
            tsv_error_val = tsv_error if 'tsv_error' in locals() else None
            best_hit_val = best_hit if 'best_hit' in locals() else None
            return sequences, tsv_data, tsv_error_val, best_hit_val
        
    except Exception as e:
        # Return query sequence, None TSV, error message, and None best_hit on error
        if progress_callback:
            progress_callback(1.0, f"BLAST search failed: {str(e)}")
        return [("query", query_sequence)], None, str(e), None
        # Note: Could also re-raise: raise Exception(f"BLAST search failed: {str(e)}")


def align_sequences_ebi_muscle(sequences, output_file=None, progress_callback=None):
    """
    Align sequences using EBI MUSCLE REST API.
    
    Args:
        sequences: List of (id, sequence) tuples
        output_file: Optional output file path
        progress_callback: Optional callback for progress updates
        
    Returns:
        str: Path to alignment file
        
    Raises:
        ValueError: If less than 2 sequences are provided
    """
    # Validate that we have at least 2 sequences
    if not sequences or len(sequences) < 2:
        raise ValueError(
            f"At least 2 sequences are required for alignment. "
            f"Found {len(sequences) if sequences else 0} sequence(s). "
            f"Please check your BLAST search parameters or try a different query sequence."
        )
    
    if progress_callback:
        progress_callback(0.2, f"Preparing {len(sequences)} sequences for alignment...")
    
    try:
        # Prepare sequences in FASTA format
        fasta_data = ""
        for seq_id, sequence in sequences:
            # Clean sequence ID for FASTA format
            clean_id = seq_id.replace(' ', '_').replace('|', '_')
            fasta_data += f">{clean_id}\n{sequence}\n"
        
        if progress_callback:
            progress_callback(0.4, "Submitting alignment to EBI MUSCLE API...")
        
        # EBI MUSCLE REST API endpoint
        url = "https://www.ebi.ac.uk/Tools/services/rest/muscle/run"
        
        # Submit alignment job
        response = requests.post(
            url,
            data={
                'email': EBI_API_PARAMS['email'],
                'sequence': fasta_data,
                'format': 'fasta',
                'stype': 'protein'
            },
            headers={'Content-Type': 'application/x-www-form-urlencoded'},
            timeout=EBI_API_PARAMS['alignment_timeout']
        )
        
        if response.status_code != 200:
            raise Exception(f"EBI API error: {response.status_code} - {response.text}")
        
        job_id = response.text.strip()
        
        if progress_callback:
            progress_callback(0.6, f"Job submitted (ID: {job_id[:8]}...), waiting for results...")
        
        # Poll for results (EBI API is usually fast, but we poll to be safe)
        status_url = f"https://www.ebi.ac.uk/Tools/services/rest/muscle/status/{job_id}"
        result_url = f"https://www.ebi.ac.uk/Tools/services/rest/muscle/result/{job_id}/aln-fasta"
        
        max_attempts = int(EBI_API_PARAMS['alignment_timeout'] / 5)  # Based on poll interval
        attempt = 0
        
        while attempt < max_attempts:
            status_response = requests.get(status_url, timeout=30)
            
            if status_response.status_code == 200:
                status = status_response.text.strip()
                
                if status == "FINISHED":
                    if progress_callback:
                        progress_callback(0.9, "Alignment complete, downloading results...")
                    
                    # Get alignment result
                    result_response = requests.get(result_url, timeout=30)
                    
                    if result_response.status_code == 200:
                        alignment_text = result_response.text
                        
                        # Save to file
                        if output_file is None:
                            temp_output = tempfile.NamedTemporaryFile(
                                mode='w', suffix='.aln', delete=False
                            )
                            output_file = temp_output.name
                            temp_output.close()
                        
                        # Write alignment (EBI returns FASTA format)
                        with open(output_file, 'w') as f:
                            f.write(alignment_text)
                        
                        # Convert to SCHEMA format (space-separated)
                        schema_file = convert_fasta_to_schema_format(output_file)
                        
                        if progress_callback:
                            progress_callback(1.0, "Alignment saved")
                        
                        return schema_file
                    else:
                        raise Exception(f"Failed to retrieve alignment: {result_response.status_code}")
                
                elif status == "ERROR" or status == "NOT_FOUND":
                    raise Exception(f"Alignment job failed with status: {status}")
                
                # Still running, wait and retry
                if progress_callback:
                    progress_callback(0.6 + (attempt / max_attempts) * 0.3, 
                                    f"Alignment in progress... (attempt {attempt + 1}/{max_attempts})")
                
                time.sleep(5)
                attempt += 1
            else:
                raise Exception(f"Failed to check job status: {status_response.status_code}")
        
        raise Exception("Alignment timed out after 5 minutes")
        
    except requests.exceptions.Timeout:
        raise Exception("Alignment API request timed out. Please try again.")
    except requests.exceptions.RequestException as e:
        raise Exception(f"Network error during alignment: {str(e)}")
    except Exception as e:
        raise Exception(f"Alignment failed: {str(e)}")


def convert_fasta_to_schema_format(fasta_file):
    """
    Convert FASTA alignment to SCHEMA format (space-separated: ID SEQUENCE).
    
    Args:
        fasta_file: Path to FASTA alignment file
        
    Returns:
        str: Path to SCHEMA-formatted alignment file
    """
    try:
        # Read FASTA alignment
        alignment = AlignIO.read(fasta_file, "fasta")
        
        # Create SCHEMA format file
        schema_file = fasta_file.replace('.fasta', '_schema.txt').replace('.aln', '_schema.txt')
        if not schema_file.endswith('_schema.txt'):
            schema_file = fasta_file + '_schema.txt'
        
        with open(schema_file, 'w') as f:
            for record in alignment:
                seq = str(record.seq)
                seq_id = record.id
                # Clean up sequence ID (remove special characters that might cause issues)
                clean_id = seq_id.replace('|', '_').replace(' ', '_').replace('.', '_')
                f.write(f"{clean_id} {seq}\n")
        
        return schema_file
        
    except Exception as e:
        raise Exception(f"Alignment conversion failed: {str(e)}")


def download_alphafold_structure(uniprot_id, progress_callback=None):
    """
    Download AlphaFold structure for a UniProt ID.
    Uses version discovery to find the latest available model.
    
    Args:
        uniprot_id: UniProt ID (e.g., "P12345")
        progress_callback: Optional callback for progress updates
        
    Returns:
        tuple: (uniprot_id, pdb_file_path) or (None, None) if not found
    """
    if progress_callback:
        progress_callback(0.1, f"Downloading AlphaFold structure for {uniprot_id}...")
    
    try:
        # Discover AlphaFold model version
        cached_version = None
        
        if progress_callback:
            progress_callback(0.2, f"Discovering AlphaFold model version...")
        
        for version in range(ALPHAFOLD_PARAMS['max_version'], ALPHAFOLD_PARAMS['min_version'] - 1, -1):
            alphafold_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{version}.pdb"
            try:
                # Use HEAD request first to check if file exists (faster)
                head_response = requests.head(alphafold_url, timeout=10)
                if head_response.status_code == 200:
                    cached_version = version
                    if progress_callback:
                        progress_callback(0.5, f"Found AlphaFold model version v{version}")
                    break
            except Exception:
                # Try next version
                continue
        
        if cached_version is None:
            if progress_callback:
                progress_callback(1.0, f"No AlphaFold structure found for {uniprot_id} (tried v{ALPHAFOLD_PARAMS['min_version']}-v{ALPHAFOLD_PARAMS['max_version']})")
            return None, None
        
        # Download the structure
        alphafold_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{cached_version}.pdb"
        
        if progress_callback:
            progress_callback(0.7, f"Downloading structure (v{cached_version})...")
        
        get_response = requests.get(alphafold_url, timeout=ALPHAFOLD_PARAMS['timeout'])
        
        if get_response.status_code == 200 and get_response.content.strip().startswith(b'HEADER'):
            # Save to temp file
            temp_pdb = tempfile.NamedTemporaryFile(mode='wb', suffix='.pdb', delete=False)
            temp_pdb.write(get_response.content)
            temp_pdb.close()
            
            if progress_callback:
                progress_callback(1.0, f"AlphaFold structure {uniprot_id} downloaded (v{cached_version})")
            
            return uniprot_id, temp_pdb.name
        else:
            if progress_callback:
                progress_callback(1.0, f"AlphaFold download failed for {uniprot_id} (v{cached_version}): {get_response.status_code}")
            return None, None
            
    except Exception as e:
        if progress_callback:
            progress_callback(1.0, f"AlphaFold structure download failed: {str(e)}")
        return None, None


def find_best_structure_for_query(query_sequence, sequences, progress_callback=None):
    """
    Find the best AlphaFold structure for the query sequence.
    
    Strategy:
    1. If query has a UniProt ID, try that first
    2. Otherwise, use the sequence with highest identity to query
    3. Fall back to first non-query sequence
    
    Args:
        query_sequence: Query protein sequence
        sequences: List of (seq_id, sequence) tuples from BLAST
        progress_callback: Optional callback for progress updates
        
    Returns:
        tuple: (uniprot_id, pdb_file_path) or (None, None) if not found
    """
    if not sequences or len(sequences) < 2:
        if progress_callback:
            progress_callback(1.0, "No sequences available for structure lookup")
        return None, None
    
    # Try query sequence first (if it has a UniProt ID format)
    # Query is typically at index 0
    if sequences[0][0] != "query":
        # If first sequence is not query, check if it's a valid UniProt ID
        first_id = sequences[0][0]
        if first_id and not first_id.startswith("query"):
            uniprot_id = first_id.split('_')[0]
            result = download_alphafold_structure(uniprot_id, progress_callback)
            if result[1]:  # If structure found
                return result
    
    # Find sequence with highest identity to query
    best_match = None
    best_identity = 0.0
    
    for seq_id, seq in sequences[1:]:  # Skip query
        identity = calculate_identity(query_sequence, seq)
        if identity > best_identity:
            best_identity = identity
            best_match = (seq_id, seq)
    
    if best_match:
        uniprot_id = best_match[0].split('_')[0]  # Get clean UniProt ID
        return download_alphafold_structure(uniprot_id, progress_callback)
    
    # Fall back to first non-query sequence
    if len(sequences) > 1:
        uniprot_id = sequences[1][0].split('_')[0]
        return download_alphafold_structure(uniprot_id, progress_callback)
    
    return None, None


def find_alphafold_structure(sequence, progress_callback=None):
    """
    DEPRECATED: Use download_alphafold_structure() with UniProt IDs from BLAST search instead.
    This function is kept for backward compatibility but should not be used.
    
    Find AlphaFold structure for a sequence using EBI BLAST against AlphaFold database.
    AlphaFold structures are complete (no missing residues) unlike crystal structures.
    
    Args:
        sequence: Protein sequence string
        progress_callback: Optional callback for progress updates
        
    Returns:
        tuple: (uniprot_id, pdb_file_path) or (None, None) if not found
    """
    if progress_callback:
        progress_callback(0.1, "Submitting BLAST job to EBI AlphaFold database...")
    
    try:
        # Submit EBI BLAST job against AlphaFold database
        url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
        data = {
            "email": "johnbmcarthur@gmail.com",
            "program": "blastp",
            "matrix": "BLOSUM62",
            "alignments": 1000,
            "scores": 1000,
            "exp": 10,
            "filter": "F",
            "gapalign": "true",
            "compstats": "F",
            "align": 0,
            "stype": "protein",
            "sequence": sequence,
            "database": "afdb"  # AlphaFold database
        }
        
        response = requests.post(url, data=data, timeout=30)
        if not response.ok or "<" in response.text:
            raise Exception(f"Failed to submit BLAST job: {response.text}")
        
        job_id = response.text.strip()
        
        if progress_callback:
            progress_callback(0.2, f"BLAST job submitted (ID: {job_id[:8]}...), waiting for results...")
        
        # Poll for job completion
        status_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}"
        start_time = time.time()
        max_wait_time = 1800  # 30 minutes
        poll_interval = 5  # seconds
        
        for i in range(int(max_wait_time / poll_interval)):
            status_response = requests.get(status_url, timeout=30)
            if status_response.status_code != 200:
                raise Exception(f"Failed to check job status: {status_response.status_code}")
            
            status = status_response.text.strip()
            elapsed = time.time() - start_time
            
            if progress_callback:
                elapsed_min = int(elapsed // 60)
                elapsed_sec = int(elapsed % 60)
                progress_callback(0.2 + (i / (max_wait_time / poll_interval)) * 0.3, 
                                f"BLAST status: {status} (elapsed: {elapsed_min}m {elapsed_sec}s)")
            
            if status == "FINISHED":
                break
            elif status in ("RUNNING", "PENDING", "QUEUED"):
                time.sleep(poll_interval)
            else:
                raise Exception(f"EBI BLAST job failed with status: {status}")
        else:
            raise TimeoutError("EBI BLAST job did not finish in time")
        
        if progress_callback:
            progress_callback(0.6, "BLAST completed, retrieving results...")
        
        # Get BLAST results in XML format
        result_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/xml"
        result_response = requests.get(result_url, timeout=30)
        
        if result_response.status_code != 200:
            raise Exception(f"Failed to retrieve BLAST results: {result_response.status_code}")
        
        if progress_callback:
            progress_callback(0.7, "Parsing BLAST XML results...")
        
        # Parse XML results
        hits = []
        try:
            root = ElementTree.fromstring(result_response.text)
            ns = {'ebi': 'http://www.ebi.ac.uk/schema'}
            
            for hit in root.findall('.//ebi:hit', ns):
                uniprot_id = hit.attrib.get('ac')  # UniProt accession
                afdb_id = hit.attrib.get('id')  # AlphaFold DB ID (e.g., AF-A0A2P6VIL6-F1)
                length = int(hit.attrib.get('length', 0))
                
                align = hit.find('.//ebi:alignment', ns)
                if align is not None:
                    identity_elem = align.find('ebi:identity', ns)
                    if identity_elem is not None:
                        identity = float(identity_elem.text)
                        
                        # Filter by >50% identity
                        if identity > 50.0 and uniprot_id:
                            hits.append({
                                'uniprot_id': uniprot_id,
                                'afdb_id': afdb_id,
                                'length': length,
                                'identity': identity
                            })
        except Exception as e:
            raise Exception(f"Failed to parse BLAST XML: {str(e)}")
        
        if not hits:
            if progress_callback:
                progress_callback(1.0, "No hits found with >50% identity")
            return None, None
        
        # Sort by identity (highest first) and get best hit
        hits.sort(key=lambda x: x['identity'], reverse=True)
        best_hit = hits[0]
        best_uniprot_id = best_hit['uniprot_id']
        
        if progress_callback:
            progress_callback(0.8, f"Best hit: {best_uniprot_id} ({best_hit['identity']:.1f}% identity)")
            progress_callback(0.85, f"Downloading AlphaFold structure for {best_uniprot_id}...")
        
        # Download AlphaFold structure with version discovery
        # Try versions from v20 down to v1 to find the latest available
        max_version = 20
        downloaded = False
        
        for version in range(max_version, 0, -1):
            alphafold_url = f"https://alphafold.ebi.ac.uk/files/AF-{best_uniprot_id}-F1-model_v{version}.pdb"
            
            # Use HEAD request first to check if file exists (faster)
            try:
                head_response = requests.head(alphafold_url, timeout=10)
                if head_response.status_code == 200:
                    # File exists, download it
                    get_response = requests.get(alphafold_url, timeout=30)
                    if get_response.status_code == 200 and get_response.content.strip().startswith(b'HEADER'):
                        # Save to temp file
                        temp_pdb = tempfile.NamedTemporaryFile(mode='wb', suffix='.pdb', delete=False)
                        temp_pdb.write(get_response.content)
                        temp_pdb.close()
                        
                        if progress_callback:
                            progress_callback(1.0, f"AlphaFold structure {best_uniprot_id} downloaded (v{version})")
                        
                        return best_uniprot_id, temp_pdb.name
            except Exception:
                # Try next version
                continue
        
        if progress_callback:
            progress_callback(1.0, f"AlphaFold structure not available for {best_uniprot_id}")
        return None, None
            
    except Exception as e:
        if progress_callback:
            progress_callback(1.0, f"AlphaFold search failed: {str(e)}")
        return None, None


# Keep old function name for backward compatibility, but redirect to AlphaFold
def find_pdb_structure(sequence, progress_callback=None):
    """
    Find structure for a sequence. Now uses AlphaFold by default.
    Redirects to find_alphafold_structure for complete structures.
    """
    return find_alphafold_structure(sequence, progress_callback)


# Removed create_pdb_parent_alignment() - PDB alignment is now handled automatically 
# by schemarecomb's PDBStructure.renumber() method. No manual alignment needed.
