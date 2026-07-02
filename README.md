# SCHEMA-RASPP Streamlit Application

A web-based interface for the SCHEMA-RASPP protein recombination library design tool.

## Overview

This Streamlit application provides an intuitive interface for protein engineers to design recombination libraries using structure-guided approaches. The app bundles the SCHEMA-RASPP modules and provides:

- **SCHEMA Energy Calculation**: Calculate disruption energies from protein structures
- **RASPP Library Design**: Find optimal crossover points for library design
- **Results Visualization**: Interactive plots and export capabilities

## Installation

1. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Run the application:**
   ```bash
   streamlit run app.py
   ```

3. **Open your browser** to the URL shown in the terminal (typically `http://localhost:8501`)

## Usage

### Page 1: SCHEMA Energy Calculation

1. Upload a PDB structure file
2. Upload a multiple sequence alignment (MSA) file
3. Optionally upload a PDB-parent alignment file
4. Set contact distance threshold (default: 5.0 Å)
5. Click "Calculate SCHEMA Contacts"
6. Optionally upload crossover points file and calculate energies

### Page 2: RASPP Library Design

1. Use SCHEMA contacts from Page 1, or upload MSA and contact files
2. Set number of crossovers and minimum fragment diversity
3. Click "Run RASPP Algorithm"
4. View optimal crossover designs

### Page 3: Results Visualization

1. View RASPP curves and energy distributions
2. Compare different designs
3. Export results as CSV, JSON, or text files

## File Formats

### PDB File
Standard Protein Data Bank format (.pdb)

### MSA File
Multiple sequence alignment in ALN format:
```
SEQ1  MKTAYIAKQR...
SEQ2  MKTAYIAKQR...
```

### Contact File
Space-separated contact positions:
```
i j ri rj
```

### Crossover File
Space-separated crossover positions:
```
10 20 30
```

## Example Files

Example files are provided in the `examples/` directory:
- `1G68.pdb` - Example PDB structure
- `lac-msa.txt` - Example MSA file
- `PSE4-1G68.txt` - Example crossover points

## Dependencies

- streamlit >= 1.28.0
- numpy >= 1.24.0
- matplotlib >= 3.7.0
- plotly >= 5.14.0
- pandas >= 2.0.0

## Features

### Automated Workflow
- **Single Sequence Input**: Enter a protein sequence and the app will:
  - Automatically find similar sequences via EBI BLAST against AlphaFold database
  - Select diverse sequences (randomly, with <90% identity between pairs)
  - Align sequences using EBI MUSCLE API (no local installation needed)
  - Find best matching AlphaFold structure automatically
  - Calculate SCHEMA contacts
- **Multi-Fragment Testing**: Automatically test fragment counts from 5-20 to find optimal library designs
- **Progress Tracking**: Real-time progress updates with time estimates

### Manual Workflow
- Upload your own PDB, MSA, and crossover files
- Full control over all parameters

### Code Organization
- **Centralized Configuration**: All constants and session state keys in `utils/config.py`
- **Session Management**: Proper initialization and cleanup via `utils/session_manager.py`
- **Temp File Management**: Automatic cleanup of temporary files
- **Error Handling**: Comprehensive error handling throughout
- **Progress Feedback**: Time estimates for long-running operations

## Notes

- The SCHEMA-RASPP modules are bundled in the `schema_raspp/` package
- Original SCHEMA-RASPP code is from: https://github.com/mattasmith/SCHEMA-RASPP
- Alignment uses EBI MUSCLE REST API (no local MUSCLE installation required)
- BLAST searches use EBI BLAST API against AlphaFold database (more reliable than NCBI)
- Structures are downloaded from AlphaFold database (complete structures, no missing residues)
- Sequence selection ensures diversity (<90% identity between selected sequences)

## References

- Voigt, C. et al., "Protein building blocks preserved by recombination," Nature Structural Biology 9(7):553-558 (2002)
- Endelman, J. et al., "Site-directed protein recombination as a shortest-path problem," Protein Engineering, Design & Selection 17(7):589-594 (2005)

## License

GPL-3.0 (inherited from SCHEMA-RASPP)
