"""
File upload and validation utilities for SCHEMA-RASPP Streamlit app.
"""

import os
import tempfile
from pathlib import Path
from typing import Optional, Tuple

from utils.config import MAX_UPLOAD_BYTES


def _safe_upload_basename(name: str) -> str:
    """Return a single path segment safe for writing under a temp directory."""
    base = Path(name or "upload").name
    if base in (".", "..") or not base:
        return "upload.bin"
    return base


def _check_upload_size(uploaded_file) -> Tuple[bool, Optional[str]]:
    size = getattr(uploaded_file, "size", None)
    if size is None:
        try:
            pos = uploaded_file.tell()
            uploaded_file.seek(0, os.SEEK_END)
            size = uploaded_file.tell()
            uploaded_file.seek(pos)
        except Exception:
            return True, None
    if size and size > MAX_UPLOAD_BYTES:
        mb = MAX_UPLOAD_BYTES / (1024 * 1024)
        return False, f"File exceeds maximum upload size ({mb:.0f} MB)"
    return True, None


def validate_pdb_file(uploaded_file):
    """
    Validate PDB file format.
    
    Args:
        uploaded_file: Streamlit UploadedFile object
        
    Returns:
        tuple: (is_valid, error_message)
    """
    if uploaded_file is None:
        return False, "No file uploaded"

    ok, err = _check_upload_size(uploaded_file)
    if not ok:
        return False, err
    
    # Check file extension
    if not uploaded_file.name.lower().endswith('.pdb'):
        return False, "File must have .pdb extension"
    
    # Read first few lines to check format
    try:
        content = uploaded_file.read()
        uploaded_file.seek(0)  # Reset file pointer
        
        # Check for PDB header indicators
        content_str = content.decode('utf-8', errors='ignore')
        lines = content_str.split('\n')[:10]
        
        # Look for common PDB record types
        pdb_indicators = ['ATOM', 'HETATM', 'HEADER', 'TITLE', 'REMARK']
        has_pdb_content = any(
            any(line.startswith(indicator) for indicator in pdb_indicators)
            for line in lines
        )
        
        if not has_pdb_content and len(lines) > 0:
            return False, "File does not appear to be a valid PDB file"
        
        return True, None
    except Exception as e:
        return False, f"Error reading file: {str(e)}"


def validate_msa_file(uploaded_file):
    """
    Validate multiple sequence alignment file format.
    
    Args:
        uploaded_file: Streamlit UploadedFile object
        
    Returns:
        tuple: (is_valid, error_message)
    """
    if uploaded_file is None:
        return False, "No file uploaded"

    ok, err = _check_upload_size(uploaded_file)
    if not ok:
        return False, err
    
    # Check file extension
    valid_extensions = ['.txt', '.msa', '.fasta', '.fa', '.aln']
    if not any(uploaded_file.name.lower().endswith(ext) for ext in valid_extensions):
        return False, f"File must have one of these extensions: {', '.join(valid_extensions)}"
    
    # Read content to check format
    try:
        content = uploaded_file.read()
        uploaded_file.seek(0)  # Reset file pointer
        
        content_str = content.decode('utf-8', errors='ignore')
        lines = [line.strip() for line in content_str.split('\n') if line.strip()]
        
        if len(lines) < 2:
            return False, "MSA file must contain at least 2 sequences"
        
        # Check for sequence identifiers (lines starting with > for FASTA or sequence names)
        has_identifiers = any(
            line.startswith('>') or 
            (not line.startswith('-') and any(c.isalpha() for c in line[:20]))
            for line in lines[:10]
        )
        
        if not has_identifiers:
            return False, "MSA file does not appear to contain sequence identifiers"
        
        return True, None
    except Exception as e:
        return False, f"Error reading file: {str(e)}"


def validate_crossover_file(uploaded_file):
    """
    Validate crossover points file format.
    
    Args:
        uploaded_file: Streamlit UploadedFile object
        
    Returns:
        tuple: (is_valid, error_message)
    """
    if uploaded_file is None:
        return False, "No file uploaded"

    ok, err = _check_upload_size(uploaded_file)
    if not ok:
        return False, err
    
    # Check file extension
    if not uploaded_file.name.lower().endswith('.txt'):
        return False, "Crossover file must have .txt extension"
    
    # Read content to check format
    try:
        content = uploaded_file.read()
        uploaded_file.seek(0)  # Reset file pointer
        
        content_str = content.decode('utf-8', errors='ignore')
        lines = [line.strip() for line in content_str.split('\n') if line.strip()]
        
        if len(lines) == 0:
            return False, "Crossover file is empty"
        
        # Check if lines contain numbers (crossover positions)
        has_numbers = False
        for line in lines[:10]:
            try:
                # Try to parse as numbers
                numbers = [int(x) for x in line.split() if x.isdigit()]
                if numbers:
                    has_numbers = True
                    break
            except:
                continue
        
        if not has_numbers:
            return False, "Crossover file does not appear to contain valid crossover positions"
        
        return True, None
    except Exception as e:
        return False, f"Error reading file: {str(e)}"


def validate_contact_file(uploaded_file):
    """Validate SCHEMA contact file (space-separated residue pairs)."""
    if uploaded_file is None:
        return False, "No file uploaded"

    ok, err = _check_upload_size(uploaded_file)
    if not ok:
        return False, err

    if not uploaded_file.name.lower().endswith(".txt"):
        return False, "Contact file must have .txt extension"

    try:
        content = uploaded_file.read()
        uploaded_file.seek(0)
        lines = [
            line.strip()
            for line in content.decode("utf-8", errors="ignore").splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]
        if not lines:
            return False, "Contact file is empty"
        parsed = 0
        for line in lines[:50]:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                parsed += 1
        if parsed == 0:
            return False, "Contact file does not contain valid residue pair lines"
        return True, None
    except Exception as e:
        return False, f"Error reading file: {str(e)}"


def save_uploaded_file(uploaded_file, directory=None):
    """
    Save uploaded file to temporary directory.
    
    Args:
        uploaded_file: Streamlit UploadedFile object
        directory: Optional directory path (defaults to temp directory)
        
    Returns:
        str: Path to saved file
    """
    if uploaded_file is None:
        return None
    
    if directory is None:
        directory = tempfile.gettempdir()
    
    # Create directory if it doesn't exist
    os.makedirs(directory, exist_ok=True)
    
    # Save file with sanitized name (prevent path traversal)
    safe_name = _safe_upload_basename(uploaded_file.name)
    file_path = os.path.join(directory, safe_name)
    with open(file_path, 'wb') as f:
        f.write(uploaded_file.getbuffer())
    
    return file_path
