"""
Visualization functions for SCHEMA-RASPP results.
Creates plots for energy maps, RASPP curves, and library results.
"""

import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

try:
    from schema_raspp import rasppcurve, schema
except ImportError:
    import rasppcurve, schema


def _fragments_from_crossovers(crossovers, parent_length):
    """
    Get fragment (start, end) list from crossover positions.
    Uses schema.getFragments when available, else a simple fallback.
    """
    try:
        return schema.getFragments(crossovers, 'A' * parent_length)
    except Exception:
        fragments = []
        start = 0
        for crossover in sorted(crossovers):
            fragments.append((start, crossover))
            start = crossover
        fragments.append((start, parent_length))
        return fragments


def plot_contact_map(contacts, num_residues=None, title="SCHEMA Contact Map"):
    """
    Create a contact map visualization.
    
    Args:
        contacts: List of contact tuples (i, j, ri, rj)
        num_residues: Number of residues (if None, inferred from contacts)
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    if num_residues is None:
        max_res = max(max(c[0], c[1]) for c in contacts) if contacts else 0
        num_residues = max_res + 1
    
    # Create contact matrix
    matrix = np.zeros((num_residues, num_residues))
    for (i, j, ri, rj) in contacts:
        if i < num_residues and j < num_residues:
            matrix[i, j] = 1
            matrix[j, i] = 1
    
    # Create Plotly heatmap
    fig = go.Figure(data=go.Heatmap(
        z=matrix,
        colorscale='Blues',
        showscale=True,
        colorbar=dict(title="Contact")
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="Residue Position",
        yaxis_title="Residue Position",
        width=700,
        height=700
    )
    
    return fig


def plot_energy_distribution(energy_results, title="SCHEMA Energy Distribution"):
    """
    Plot distribution of SCHEMA disruption energies.
    
    Args:
        energy_results: List of dictionaries with 'disruption' key
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    disruptions = [r['disruption'] for r in energy_results if r.get('disruption') is not None]
    
    if not disruptions:
        return None
    
    fig = go.Figure()
    
    # Histogram
    fig.add_trace(go.Histogram(
        x=disruptions,
        nbinsx=30,
        name="Energy Distribution",
        marker_color='steelblue'
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="SCHEMA Disruption Energy",
        yaxis_title="Frequency",
        width=700,
        height=500
    )
    
    return fig


def plot_raspp_curve(raspp_results, parents, bin_width=1.0, title="RASPP Curve"):
    """
    Generate RASPP curve showing average energy vs average mutation.
    
    Args:
        raspp_results: List of RASPP result tuples (avg_energy, crossovers, l_min, l_max)
        parents: List of parent sequences
        bin_width: Width of mutation bins
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    try:
        # Use rasppcurve module if available
        curve_data = rasppcurve.curve(raspp_results, parents, bin_width)
        
        if not curve_data:
            return None
        
        # Extract data
        avg_mutations = [d[0] for d in curve_data]
        avg_energies = [d[1] for d in curve_data]
        counts = [d[2] for d in curve_data]
        
        # Create Plotly scatter plot
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=avg_mutations,
            y=avg_energies,
            mode='lines+markers',
            name="RASPP Curve",
            marker=dict(size=8, color='steelblue'),
            line=dict(width=2, color='steelblue'),
            text=[f"Count: {c}" for c in counts],
            hovertemplate="Avg Mutations: %{x}<br>Avg Energy: %{y}<br>Count: %{text}<extra></extra>"
        ))
        
        fig.update_layout(
            title=title,
            xaxis_title="Average Mutations",
            yaxis_title="Average SCHEMA Energy",
            width=800,
            height=600
        )
        
        return fig
        
    except Exception as e:
        # Fallback: simple scatter plot
        avg_energies = [r[0] for r in raspp_results]
        # Estimate mutations from fragment lengths
        mutations = []
        for r in raspp_results:
            _, crossovers, l_min, l_max = r
            # Rough estimate based on fragment lengths
            avg_mut = (l_min + l_max) / 2.0
            mutations.append(avg_mut)
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=mutations,
            y=avg_energies,
            mode='markers',
            name="RASPP Results",
            marker=dict(size=8, color='steelblue')
        ))
        
        fig.update_layout(
            title=title,
            xaxis_title="Estimated Mutations",
            yaxis_title="Average SCHEMA Energy",
            width=800,
            height=600
        )
        
        return fig


def plot_crossover_comparison(raspp_results, title="Crossover Comparison"):
    """
    Compare different crossover designs.
    
    Args:
        raspp_results: List of RASPP result dictionaries
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    if not raspp_results:
        return None
    
    # Extract data
    energies = [r['average_energy'] for r in raspp_results]
    num_fragments = [r['num_fragments'] for r in raspp_results]
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=num_fragments,
        y=energies,
        mode='markers',
        name="Designs",
        marker=dict(
            size=10,
            color=energies,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="Energy")
        ),
        text=[f"Crossovers: {len(r['crossovers'])}" for r in raspp_results],
        hovertemplate="Fragments: %{x}<br>Energy: %{y}<br>%{text}<extra></extra>"
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="Number of Fragments",
        yaxis_title="Average SCHEMA Energy",
        width=700,
        height=500
    )
    
    return fig


def plot_fragment_lengths(crossovers, parent_length, title="Fragment Length Distribution"):
    """
    Visualize fragment lengths for a given crossover design.
    
    Args:
        crossovers: List of crossover positions
        parent_length: Length of parent sequence
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    fragments = _fragments_from_crossovers(crossovers, parent_length)
    lengths = [end - start for start, end in fragments]
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=[f"Fragment {i+1}" for i in range(len(lengths))],
        y=lengths,
        marker_color='steelblue',
        text=lengths,
        textposition='outside'
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="Fragment",
        yaxis_title="Length (residues)",
        width=700,
        height=400
    )
    
    return fig


def create_results_dataframe(raspp_results):
    """
    Convert RASPP results to a pandas DataFrame.
    
    Args:
        raspp_results: List of RASPP result dictionaries
        
    Returns:
        pandas.DataFrame: DataFrame with results
    """
    data = []
    for i, result in enumerate(raspp_results):
        data.append({
            'Design': i + 1,
            'Average Energy': result.get('average_energy', 0),
            'Number of Crossovers': len(result.get('crossovers', [])),
            'Number of Fragments': result.get('num_fragments', 0),
            'Min Fragment Length': result.get('min_fragment_length', 0),
            'Max Fragment Length': result.get('max_fragment_length', 0),
            'Crossovers': ', '.join(map(str, result.get('crossovers', [])))
        })
    
    return pd.DataFrame(data)


def plot_sequence_with_crossovers(crossovers, parent_length, parent_sequence=None, 
                                  title="Fragmentation Sites", contacts=None):
    """
    Create a linear sequence view with crossover positions marked.
    Shows fragments as colored segments with crossover sites as vertical lines.
    
    Args:
        crossovers: List of crossover positions (1-based)
        parent_length: Length of parent sequence
        parent_sequence: Optional sequence string for display
        contacts: Optional list of contacts for overlay
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    fragments = _fragments_from_crossovers(crossovers, parent_length)
    
    fig = go.Figure()
    
    # Color each fragment differently
    colors = px.colors.qualitative.Set3
    y_pos = 0
    
    for i, (start, end) in enumerate(fragments):
        fragment_length = end - start
        fig.add_trace(go.Bar(
            x=[fragment_length],
            y=[y_pos],
            orientation='h',
            base=start,
            marker_color=colors[i % len(colors)],
            name=f'Fragment {i+1}',
            text=[f'F{i+1}'],
            textposition='inside',
            hovertemplate=f'Fragment {i+1}<br>Position: {start}-{end}<br>Length: {fragment_length} residues<extra></extra>'
        ))
    
    # Add crossover markers
    for i, crossover in enumerate(sorted(crossovers)):
        fig.add_vline(
            x=crossover,
            line_dash="dash",
            line_color="red",
            line_width=2,
            annotation_text=f"C{i+1}",
            annotation_position="top",
            annotation_font_size=10
        )
    
    # Add contact lines if provided
    if contacts:
        for i, j, _, _ in contacts[:200]:  # Limit for performance
            if i < parent_length and j < parent_length:
                fig.add_shape(
                    type="line",
                    x0=i, x1=j,
                    y0=-0.1, y1=-0.1,
                    line=dict(color="rgba(173, 216, 230, 0.3)", width=1),  # lightblue with 0.3 opacity
                    layer="below"
                )
    
    fig.update_layout(
        title=title,
        xaxis_title="Sequence Position",
        yaxis=dict(visible=False, range=[-0.5, 0.5]),
        height=200,
        showlegend=True,
        hovermode='x unified'
    )
    
    return fig


def plot_crossover_frequency_heatmap(all_crossovers, sequence_length, 
                                     title="Crossover Position Frequency"):
    """
    Create a heatmap showing how often each position is used as a crossover.
    Useful for multi-fragment testing to see preferred fragmentation sites.
    
    Args:
        all_crossovers: List of lists, each containing crossover positions for a design
        sequence_length: Length of sequence
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    # Count frequency of each position
    position_counts = [0] * sequence_length
    
    for crossovers in all_crossovers:
        for crossover in crossovers:
            if 0 <= crossover < sequence_length:
                position_counts[crossover] += 1
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=[position_counts],
        colorscale='YlOrRd',
        showscale=True,
        colorbar=dict(title="Frequency"),
        y=['Crossover Frequency'],
        hovertemplate='Position: %{x}<br>Frequency: %{z}<extra></extra>'
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="Sequence Position",
        yaxis_title="",
        height=150,
        width=800
    )
    
    return fig


def plot_contact_map_with_crossovers(contacts, crossovers, num_residues, 
                                     title="Contact Map with Fragmentation Sites"):
    """
    Overlay crossover positions on contact map.
    Shows which contacts are disrupted by fragmentation.
    
    Args:
        contacts: List of contact tuples
        crossovers: List of crossover positions
        num_residues: Number of residues
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    # Create contact matrix
    matrix = np.zeros((num_residues, num_residues))
    for (i, j, ri, rj) in contacts:
        if i < num_residues and j < num_residues:
            matrix[i, j] = 1
            matrix[j, i] = 1
    
    fig = go.Figure(data=go.Heatmap(
        z=matrix,
        colorscale='Blues',
        showscale=True,
        colorbar=dict(title="Contact"),
        hovertemplate='Residue %{x} ↔ Residue %{y}<br>Contact: %{z}<extra></extra>'
    ))
    
    # Add crossover lines
    for i, crossover in enumerate(crossovers):
        if 0 <= crossover < num_residues:
            # Vertical line
            fig.add_shape(
                type="line",
                x0=crossover, x1=crossover,
                y0=0, y1=num_residues,
                line=dict(color="red", width=2, dash="dash"),
                layer="above"
            )
            # Horizontal line
            fig.add_shape(
                type="line",
                x0=0, x1=num_residues,
                y0=crossover, y1=crossover,
                line=dict(color="red", width=2, dash="dash"),
                layer="above"
            )
    
    fig.update_layout(
        title=title,
        xaxis_title="Residue Position",
        yaxis_title="Residue Position",
        width=700,
        height=700
    )
    
    return fig


def plot_energy_landscape_with_crossovers(avg_energies, crossovers, 
                                         sequence_length, title="Energy Landscape"):
    """
    Plot energy contribution along sequence with crossover positions marked.
    Shows where fragmentation minimizes disruption.
    
    Args:
        avg_energies: List of (i, j, energy) tuples or energy matrix
        crossovers: List of crossover positions
        sequence_length: Length of sequence
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    # Calculate per-position energy (sum of contacts involving each position)
    position_energies = [0.0] * sequence_length
    
    if isinstance(avg_energies, list) and len(avg_energies) > 0:
        if isinstance(avg_energies[0], tuple) and len(avg_energies[0]) == 3:
            # List of (i, j, energy) tuples
            for i, j, energy in avg_energies:
                if i < sequence_length:
                    position_energies[i] += energy
                if j < sequence_length:
                    position_energies[j] += energy
        else:
            # Try to treat as matrix
            try:
                for i in range(min(sequence_length, len(avg_energies))):
                    if isinstance(avg_energies[i], (list, np.ndarray)):
                        for j in range(min(sequence_length, len(avg_energies[i]))):
                            position_energies[i] += float(avg_energies[i][j])
                            if i != j:
                                position_energies[j] += float(avg_energies[i][j])
            except Exception:
                pass
    
    fig = go.Figure()
    
    # Energy landscape
    fig.add_trace(go.Scatter(
        x=list(range(sequence_length)),
        y=position_energies,
        mode='lines',
        name='Energy',
        line=dict(color='steelblue', width=2),
        fill='tozeroy',
        fillcolor='rgba(70, 130, 180, 0.3)',
        hovertemplate='Position: %{x}<br>Energy: %{y:.2f}<extra></extra>'
    ))
    
    # Crossover markers
    for i, crossover in enumerate(crossovers):
        if 0 <= crossover < sequence_length:
            fig.add_vline(
                x=crossover,
                line_dash="dash",
                line_color="red",
                line_width=2,
                annotation_text=f"XO{i+1}",
                annotation_position="top",
                annotation_font_size=10
            )
    
    fig.update_layout(
        title=title,
        xaxis_title="Sequence Position",
        yaxis_title="Energy Contribution",
        height=400,
        hovermode='x unified'
    )
    
    return fig


def plot_fragment_energy_contribution(crossovers, contacts, parent_length,
                                     title="Fragment Energy Contribution"):
    """
    Calculate and visualize energy contribution of each fragment.
    Shows which fragments contribute most to disruption.
    
    Args:
        crossovers: List of crossover positions
        contacts: List of contact tuples
        parent_length: Length of parent sequence
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    fragments = _fragments_from_crossovers(crossovers, parent_length)
    
    # Calculate energy per fragment
    fragment_energies = []
    fragment_lengths = []
    
    for start, end in fragments:
        length = end - start
        fragment_lengths.append(length)
        
        # Count contacts within this fragment
        fragment_contacts = [
            c for c in contacts 
            if start <= c[0] < end and start <= c[1] < end
        ]
        fragment_energies.append(len(fragment_contacts))
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=[f"Fragment {i+1}" for i in range(len(fragments))],
        y=fragment_energies,
        marker_color='steelblue',
        text=fragment_energies,
        textposition='outside',
        hovertemplate='Fragment %{x}<br>Contacts: %{y}<br>Length: %{customdata} residues<extra></extra>',
        customdata=fragment_lengths
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="Fragment",
        yaxis_title="Number of Contacts",
        height=400
    )
    
    return fig


def plot_multiple_designs_comparison(designs_list, parent_length,
                                     title="Design Comparison"):
    """
    Show multiple fragmentation designs overlaid for comparison.
    Each design shown as a different colored set of fragments.
    
    Args:
        designs_list: List of tuples (crossovers, energy) or dicts with 'crossovers' and 'average_energy'
        parent_length: Length of parent sequence
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    fig = go.Figure()
    
    colors = px.colors.qualitative.Set3
    y_offset = 0
    
    for design_idx, design in enumerate(designs_list):
        # Handle both tuple and dict formats
        if isinstance(design, tuple):
            crossovers, energy = design
        elif isinstance(design, dict):
            crossovers = design.get('crossovers', [])
            energy = design.get('average_energy', 0)
        else:
            continue
        
        fragments = _fragments_from_crossovers(crossovers, parent_length)
        
        for frag_idx, (start, end) in enumerate(fragments):
            fig.add_trace(go.Bar(
                x=[end - start],
                y=[y_offset],
                base=start,
                orientation='h',
                marker_color=colors[design_idx % len(colors)],
                name=f'Design {design_idx+1}',
                text=[f'D{design_idx+1}F{frag_idx+1}'],
                textposition='inside',
                hovertemplate=f'Design {design_idx+1}, Fragment {frag_idx+1}<br>'
                             f'Energy: {energy:.2f}<br>'
                             f'Position: {start}-{end}<br>'
                             f'Length: {end-start} residues<extra></extra>',
                showlegend=(frag_idx == 0)  # Only show legend for first fragment
            ))
        
        y_offset += 1
    
    # Add crossover markers for all designs
    for design_idx, design in enumerate(designs_list):
        if isinstance(design, tuple):
            crossovers, _ = design
        elif isinstance(design, dict):
            crossovers = design.get('crossovers', [])
        else:
            continue
            
        for crossover in crossovers:
            fig.add_vline(
                x=crossover,
                line_dash="dot",
                line_color=colors[design_idx % len(colors)],
                opacity=0.5,
                line_width=1
            )
    
    fig.update_layout(
        title=title,
        xaxis_title="Sequence Position",
        yaxis_title="Design",
        height=max(300, len(designs_list) * 50),
        showlegend=True,
        hovermode='x unified'
    )
    
    return fig


def plot_sequence_with_secondary_structure(crossovers, secondary_structure, 
                                           parent_length, title="Fragmentation with Secondary Structure"):
    """
    Show crossover positions relative to secondary structure elements.
    Requires secondary structure prediction or PDB data.
    
    Args:
        crossovers: List of crossover positions
        secondary_structure: List of secondary structure assignments (H, E, C, T)
        parent_length: Length of parent sequence
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    fig = go.Figure()
    
    # Map secondary structure to colors
    ss_colors = {
        'H': 'red',      # Helix
        'E': 'blue',     # Sheet
        'C': 'gray',     # Coil
        'T': 'orange',   # Turn
        'S': 'green',    # Bend
        'G': 'purple',   # 3-helix
        'I': 'cyan',     # Pi helix
    }
    
    # Ensure secondary structure matches length
    if len(secondary_structure) != parent_length:
        secondary_structure = (secondary_structure * ((parent_length // len(secondary_structure)) + 1))[:parent_length]
    
    # Plot secondary structure as colored bars
    for i, ss in enumerate(secondary_structure):
        fig.add_trace(go.Bar(
            x=[1],
            y=[i],
            base=i,
            orientation='h',
            marker_color=ss_colors.get(ss, 'gray'),
            name=ss,
            showlegend=False,
            hovertemplate=f'Position: {i}<br>SS: {ss}<extra></extra>'
        ))
    
    # Add crossover markers
    for i, crossover in enumerate(crossovers):
        if 0 <= crossover < parent_length:
            fig.add_vline(
                x=crossover,
                line_dash="dash",
                line_color="black",
                line_width=3,
                annotation_text=f"XO{i+1}",
                annotation_position="top",
                annotation_font_size=10
            )
    
    fig.update_layout(
        title=title,
        xaxis_title="Sequence Position",
        yaxis_title="Secondary Structure",
        height=200,
        xaxis=dict(range=[0, parent_length])
    )
    
    return fig


def plot_interactive_sequence_viewer(crossovers, parent_sequence, 
                                    contacts=None, title="Interactive Sequence Viewer"):
    """
    Create an interactive, zoomable sequence viewer with crossover positions.
    Allows users to explore different regions of the sequence.
    
    Args:
        crossovers: List of crossover positions
        parent_sequence: Parent sequence string
        contacts: Optional list of contacts for overlay
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    fig = go.Figure()
    
    seq_length = len(parent_sequence)
    
    # Add sequence as text annotations (show every 10th residue for readability)
    positions = []
    residues = []
    for i in range(0, seq_length, max(1, seq_length // 100)):  # Show up to 100 residues
        positions.append(i)
        residues.append(parent_sequence[i])
    
    fig.add_trace(go.Scatter(
        x=positions,
        y=[0] * len(positions),
        mode='text',
        text=residues,
        textfont=dict(size=10, family="monospace"),
        showlegend=False,
        hovertemplate='Position: %{x}<br>Residue: %{text}<extra></extra>'
    ))
    
    # Add crossover markers
    for i, crossover in enumerate(crossovers):
        if 0 <= crossover < seq_length:
            fig.add_vline(
                x=crossover,
                line_dash="dash",
                line_color="red",
                line_width=2,
                annotation_text=f"Pos {crossover}",
                annotation_position="top",
                annotation_font_size=10
            )
    
    # Add contact lines if provided
    if contacts:
        for i, j, _, _ in contacts[:100]:  # Limit for performance
            if i < seq_length and j < seq_length:
                fig.add_shape(
                    type="line",
                    x0=i, x1=j,
                    y0=-0.5, y1=-0.5,
                    line=dict(color="rgba(173, 216, 230, 0.3)", width=1),  # lightblue with 0.3 opacity
                    layer="below"
                )
    
    fig.update_layout(
        title=title,
        xaxis_title="Sequence Position",
        yaxis=dict(visible=False, range=[-1, 1]),
        height=300,
        hovermode='x unified'
    )
    
    return fig


def plot_fragment_length_vs_energy(crossovers, contacts, parent_length, 
                                   avg_energy, title="Fragment Length vs Energy"):
    """
    Scatter plot showing fragment lengths vs their energy contribution.
    Helps identify optimal fragment sizes.
    
    Args:
        crossovers: List of crossover positions
        contacts: List of contact tuples
        parent_length: Length of parent sequence
        avg_energy: Average energy for this design
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    fragments = _fragments_from_crossovers(crossovers, parent_length)
    
    fragment_lengths = []
    fragment_energies = []
    
    for start, end in fragments:
        length = end - start
        fragment_lengths.append(length)
        
        # Calculate contacts in fragment
        frag_contacts = [
            c for c in contacts
            if start <= c[0] < end and start <= c[1] < end
        ]
        fragment_energies.append(len(frag_contacts))
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=fragment_lengths,
        y=fragment_energies,
        mode='markers+text',
        marker=dict(size=12, color='steelblue'),
        text=[f"F{i+1}" for i in range(len(fragments))],
        textposition="top center",
        hovertemplate='Fragment %{text}<br>Length: %{x} residues<br>Contacts: %{y}<extra></extra>'
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="Fragment Length (residues)",
        yaxis_title="Number of Contacts",
        height=500
    )
    
    return fig


def compute_crossover_position_counts(all_crossovers, sequence_length):
    """
    Count how often each sequence index appears as a crossover (0-based indices).

    Args:
        all_crossovers: List of lists of crossover positions (0-based, same as RASPP results).
        sequence_length: Full sequence length (indices must satisfy 0 <= i < sequence_length).

    Returns:
        dict[int, int]: position (0-based) -> count. Empty dict if none in range.
    """
    position_counts = {}
    for crossovers in all_crossovers:
        for crossover in crossovers:
            if 0 <= crossover < sequence_length:
                position_counts[crossover] = position_counts.get(crossover, 0) + 1
    return position_counts


def plot_crossover_distribution(
    all_crossovers,
    sequence_length,
    title="Crossover Position Distribution",
    selected_positions_1based=None,
):
    """
    Bar chart showing distribution of crossover positions at each site.
    Shows frequency at each individual position (not binned).
    Identifies preferred fragmentation regions.
    
    Args:
        all_crossovers: List of lists, each containing crossover positions
        sequence_length: Length of sequence
        title: Plot title
        selected_positions_1based: Applied crossover positions (1-based); bars shown in red
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    position_counts = compute_crossover_position_counts(all_crossovers, sequence_length)
    if not position_counts:
        return None
    
    # Sort positions and get frequencies
    positions = sorted(position_counts.keys())
    frequencies = [position_counts[pos] for pos in positions]
    max_freq = max(frequencies)
    label_threshold = 0.1 * max_freq
    # Label sequence positions whose frequency is strictly above 10% of the peak bar
    bar_text = [str(pos) if freq > label_threshold else "" for pos, freq in zip(positions, frequencies)]

    selected_0based = {
        int(p) - 1
        for p in (selected_positions_1based or [])
        if 1 <= int(p) <= int(sequence_length)
    }
    bar_colors = [
        "#d32f2f" if pos in selected_0based else "steelblue"
        for pos in positions
    ]

    fig = go.Figure()

    # Use bar chart instead of histogram to show each site individually
    fig.add_trace(go.Bar(
        x=positions,
        y=frequencies,
        marker_color=bar_colors,
        opacity=0.85,
        hovertemplate='Position: %{x}<br>Frequency: %{y}<extra></extra>',
        width=0.8,  # Narrow bars to show individual sites
        text=bar_text,
        textposition="outside",
        cliponaxis=False,
    ))

    fig.update_layout(
        title=title,
        xaxis_title="Sequence Position",
        yaxis_title="Frequency",
        height=520,
        margin=dict(t=90, l=60, r=40, b=60),
        xaxis=dict(
            type='linear',
            range=[0, sequence_length],
        ),
    )
    
    return fig


def plot_3d_fragment_energy(all_designs, contacts, parent_length,
                            title="3D Fragment Energy Visualization"):
    """
    Create a 3D scatter plot showing fragment start position, end position, and energy.
    
    Args:
        all_designs: List of design dictionaries, each with 'crossovers' key
        contacts: List of contact tuples (i, j, ri, rj)
        parent_length: Length of parent sequence
        title: Plot title
        
    Returns:
        plotly.graph_objects.Figure: 3D scatter plot figure
    """
    if not all_designs or not contacts:
        return None
    
    # Collect all fragments from all designs
    fragment_data = []  # List of (start, end, energy, design_idx, fragment_idx)
    
    for design_idx, design in enumerate(all_designs):
        crossovers = design.get('crossovers', [])
        if not crossovers:
            continue
        
        fragments = _fragments_from_crossovers(crossovers, parent_length)
        
        # Calculate energy for each fragment
        for frag_idx, (start, end) in enumerate(fragments):
            # Count contacts that involve residues within this fragment
            # This includes contacts within the fragment and contacts crossing boundaries
            fragment_contacts = [
                c for c in contacts
                if (start <= c[0] < end or start <= c[1] < end)
            ]
            
            # Calculate energy as number of contacts involving this fragment
            # Weight by how much of the contact is within the fragment
            energy = 0.0
            for i, j, ri, rj in fragment_contacts:
                # Count full contact if both residues are in fragment
                if start <= i < end and start <= j < end:
                    energy += 1.0
                # Count partial contact if only one residue is in fragment
                elif start <= i < end or start <= j < end:
                    energy += 0.5
            
            fragment_data.append({
                'start': start,
                'end': end,
                'energy': energy,
                'design_idx': design_idx,
                'fragment_idx': frag_idx,
                'length': end - start
            })
    
    if not fragment_data:
        return None
    
    # Extract data for plotting
    starts = [f['start'] for f in fragment_data]
    ends = [f['end'] for f in fragment_data]
    energies = [f['energy'] for f in fragment_data]
    lengths = [f['length'] for f in fragment_data]
    design_indices = [f['design_idx'] for f in fragment_data]
    fragment_indices = [f['fragment_idx'] for f in fragment_data]
    
    # Create 3D scatter plot
    fig = go.Figure(data=go.Scatter3d(
        x=starts,
        y=ends,
        z=energies,
        mode='markers',
        marker=dict(
            size=5,
            color=energies,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="Fragment Energy"),
            line=dict(width=0.5, color='rgba(0,0,0,0.1)')
        ),
        text=[f"Design {d+1}, Fragment {f+1}<br>Start: {s}, End: {e}<br>Energy: {en:.2f}<br>Length: {l}" 
              for d, f, s, e, en, l in zip(design_indices, fragment_indices, starts, ends, energies, lengths)],
        hovertemplate='%{text}<extra></extra>'
    ))
    
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='Fragment Start Position',
            yaxis_title='Fragment End Position',
            zaxis_title='Fragment Energy',
            xaxis=dict(range=[0, parent_length]),
            yaxis=dict(range=[0, parent_length]),
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.2)
            )
        ),
        height=700,
        width=900
    )
    
    return fig


def format_msa_with_crossovers_html(alignment_data, crossover_positions, 
                                    max_sequences=10, chars_per_line=80):
    """
    Format multiple sequence alignment with crossover markers (carets) as HTML.
    
    Based on SCHEMA-RASPP: crossovers are 1-based, and a crossover at position X means:
    - Fragment before ends at position X-1 (1-based)
    - Fragment after starts at position X (1-based)
    - The caret should be placed BETWEEN positions X-1 and X, which is at position X-1 (0-based)
    - Or we can show it at position X (1-based) to indicate where the new fragment starts
    
    Args:
        alignment_data: List of (name, sequence) tuples or list of sequences
        crossover_positions: List of crossover positions (1-based)
        max_sequences: Maximum number of sequences to display (default: 10)
        chars_per_line: Characters to display per line (default: 80)
        
    Returns:
        str: HTML formatted alignment
    """
    if not alignment_data:
        return None
    
    # Handle different input formats
    if isinstance(alignment_data[0], tuple):
        # List of (name, sequence) tuples
        seq_names = [name for name, _ in alignment_data[:max_sequences]]
        sequences = [seq for _, seq in alignment_data[:max_sequences]]
    else:
        # List of sequences
        seq_names = [f"Sequence {i+1}" for i in range(min(len(alignment_data), max_sequences))]
        sequences = alignment_data[:max_sequences]
    
    if not sequences:
        return None
    
    alignment_length = len(sequences[0])
    
    # Find max name length for proper alignment
    max_name_len = max(len(name) for name in seq_names) if seq_names else 15
    
    # Convert crossover positions to 0-based for internal use
    # A crossover at position X (1-based) means the new fragment starts at X
    # So we mark position X-1 (0-based) as the end of previous fragment
    # For display, we'll place the caret at position X (1-based) to show where new fragment starts
    crossovers_1based = sorted(set([p for p in crossover_positions if 1 <= p <= alignment_length]))
    crossovers_0based = [p - 1 for p in crossovers_1based]  # Position where new fragment starts (0-based)
    
    # Create HTML with proper table structure for alignment
    html_parts = ['<div style="font-family: monospace; font-size: 12px; line-height: 1.4;">']
    html_parts.append('<h4>Multiple Sequence Alignment with Crossover Markers</h4>')
    html_parts.append('<div style="overflow-x: auto;">')
    html_parts.append('<table style="border-collapse: collapse; font-family: monospace; font-size: 12px;">')
    
    # Split alignment into chunks for display
    num_chunks = (alignment_length + chars_per_line - 1) // chars_per_line
    
    for chunk_idx in range(num_chunks):
        start_pos = chunk_idx * chars_per_line
        end_pos = min(start_pos + chars_per_line, alignment_length)
        chunk_length = end_pos - start_pos
        
        # Position numbers row - properly formatted with spacing
        pos_nums_row = ['<tr style="background-color: #f0f0f0;">']
        pos_nums_row.append(f'<td style="padding-right: 10px; text-align: right; width: {max_name_len}ch;">&nbsp;</td>')
        pos_nums_row.append('<td style="padding: 0;">')
        pos_nums_cells = []
        for i in range(chunk_length):
            pos = start_pos + i + 1  # 1-based position
            # Format position number with proper spacing
            if pos % 10 == 0:
                pos_str = str(pos)
            else:
                pos_str = str(pos % 10)
            pos_nums_cells.append(f'<span style="display: inline-block; width: 1ch; text-align: center; color: #666;">{pos_str}</span>')
        pos_nums_row.append(''.join(pos_nums_cells))
        pos_nums_row.append('</td></tr>')
        html_parts.append(''.join(pos_nums_row))
        
        # Caret row - place carets at crossover positions
        caret_row = ['<tr style="background-color: #fff0f0;">']
        caret_row.append(f'<td style="padding-right: 10px; text-align: right; width: {max_name_len}ch;">&nbsp;</td>')
        caret_row.append('<td style="padding: 0;">')
        caret_cells = []
        for i in range(chunk_length):
            pos_0based = start_pos + i
            pos_1based = pos_0based + 1
            # Place caret at position X (1-based) to show where new fragment starts
            if pos_1based in crossovers_1based:
                caret_cells.append('<span style="display: inline-block; width: 1ch; text-align: center; color: red; font-weight: bold;">^</span>')
            else:
                caret_cells.append('<span style="display: inline-block; width: 1ch; text-align: center;">&nbsp;</span>')
        caret_row.append(''.join(caret_cells))
        caret_row.append('</td></tr>')
        html_parts.append(''.join(caret_row))
        
        # Each sequence row - properly aligned
        for name, seq in zip(seq_names, sequences):
            seq_chunk = seq[start_pos:end_pos]
            seq_row = ['<tr>']
            # Name column with fixed width
            seq_row.append(f'<td style="padding-right: 10px; text-align: right; width: {max_name_len}ch; font-weight: bold;">{name}</td>')
            # Sequence column with character-by-character alignment
            seq_row.append('<td style="padding: 0; letter-spacing: 0;">')
            seq_cells = []
            for i, char in enumerate(seq_chunk):
                pos_0based = start_pos + i
                pos_1based = pos_0based + 1
                # Highlight position X (1-based) where new fragment starts
                if pos_1based in crossovers_1based:
                    seq_cells.append(f'<span style="display: inline-block; width: 1ch; text-align: center; background-color: #ffcccc;">{char}</span>')
                else:
                    seq_cells.append(f'<span style="display: inline-block; width: 1ch; text-align: center;">{char}</span>')
            seq_row.append(''.join(seq_cells))
            seq_row.append('</td></tr>')
            html_parts.append(''.join(seq_row))
        
        # Add spacing between chunks
        if chunk_idx < num_chunks - 1:
            html_parts.append('<tr><td colspan="2" style="height: 10px;">&nbsp;</td></tr>')
    
    html_parts.append('</table>')
    html_parts.append('</div>')
    html_parts.append('</div>')
    return ''.join(html_parts)


def compute_msa_column_conservation_diversity(sequences):
    """
    Per-column statistics for an MSA (equal-length sequences).

    Returns:
        positions: 1-based column indices
        conservation_pct: frequency of the most abundant residue × 100
        diversity_score: count of distinct characters in the column (including gap '-')
    """
    from collections import Counter

    if not sequences:
        return [], [], []
    width = len(sequences[0])
    if any(len(s) != width for s in sequences):
        raise ValueError("All sequences in the MSA must have the same length")

    nseq = len(sequences)
    positions = []
    conservation_pct = []
    diversity_score = []

    for col in range(width):
        column_chars = [seq[col] for seq in sequences]
        ctr = Counter(column_chars)
        diversity_score.append(len(ctr))
        max_cnt = max(ctr.values())
        conservation_pct.append(100.0 * max_cnt / nseq)
        positions.append(col + 1)

    return positions, conservation_pct, diversity_score


def compute_consensus_segment_midpoints(
    positions,
    conservation_pct,
    diversity_score,
    conservation_min_pct,
    diversity_max,
):
    """
    For each pair of consecutive alignment columns (x1, x2), if BOTH columns satisfy
    conservation > conservation_min_pct AND diversity < diversity_max, record the
    midpoint x = avg(x1, x2) at y = 100% (conservation axis).

    Returns:
        mid_x: list of midpoint positions (1-based scale, may be fractional)
        pair_labels: list of strings like "12–13" for hover
    """
    mid_x = []
    pair_labels = []
    n = len(positions)
    for j in range(n - 1):
        if (
            conservation_pct[j] > conservation_min_pct
            and conservation_pct[j + 1] > conservation_min_pct
            and diversity_score[j] < diversity_max
            and diversity_score[j + 1] < diversity_max
        ):
            mid_x.append((positions[j] + positions[j + 1]) / 2.0)
            pair_labels.append(f"{positions[j]}–{positions[j + 1]}")
    return mid_x, pair_labels


def allowed_full_schema_crossovers_from_segment_marks(
    parents,
    conservation_min_pct,
    diversity_max,
):
    """
    SCHEMA 1-based crossover indices where a crossover may occur at graph segment marks.

    For each consecutive column pair (p, p+1) passing the same filters as the green dots,
    the crossover between columns is identified by SCHEMA 1-based position p+1 (second column).
    """
    positions, cons, div = compute_msa_column_conservation_diversity(parents)
    allowed = set()
    for j in range(len(positions) - 1):
        if (
            cons[j] > conservation_min_pct
            and cons[j + 1] > conservation_min_pct
            and div[j] < diversity_max
            and div[j + 1] < diversity_max
        ):
            allowed.add(positions[j + 1])
    return allowed


def allowed_collapsed_starts_from_graph_segment_marks(
    parents,
    collapsed_parents,
    identical_sites,
    conservation_min_pct,
    diversity_max,
):
    """
    0-based collapsed indices where a fragment may start if the implied full SCHEMA crossover
    matches one of the segment-mark boundaries from allowed_full_schema_crossovers_from_segment_marks.
    """
    try:
        from schema_raspp.raspp import translate_collapsed_indices
    except ImportError:
        from raspp import translate_collapsed_indices

    full_allowed = allowed_full_schema_crossovers_from_segment_marks(
        parents, conservation_min_pct, diversity_max
    )
    if not full_allowed:
        return set()

    L = len(collapsed_parents[0])
    sites_sorted = sorted(identical_sites)
    allowed = set()
    for s in range(1, L):
        xb = [s + 1]
        translate_collapsed_indices(xb, list(sites_sorted))
        if xb[0] in full_allowed:
            allowed.add(s)
    return allowed


def plot_msa_conservation_diversity(
    sequences,
    title="MSA conservation and diversity",
    segment_conservation_min_pct=None,
    segment_diversity_max=None,
):
    """
    Dual-axis line plot: conservation % (left) and diversity score (right) vs alignment position.

    If segment_conservation_min_pct and segment_diversity_max are set, adds green dots at
    y=100% on the conservation axis at x = average of consecutive positions (x1, x2) where
    both columns exceed the conservation threshold and are below the diversity threshold.

    Hover shows position and the y-value for each trace.
    """
    try:
        from plotly.subplots import make_subplots
    except ImportError:
        make_subplots = None

    positions, conservation_pct, diversity_score = compute_msa_column_conservation_diversity(
        sequences
    )
    if not positions:
        return None

    fig = make_subplots(specs=[[{"secondary_y": True}]])

    fig.add_trace(
        go.Scatter(
            x=positions,
            y=conservation_pct,
            mode="lines",
            name="Conservation %",
            line=dict(color="#1f77b4", width=1.5),
            hovertemplate=(
                "Position: %{x}<br>"
                "Conservation: %{y:.2f}%<extra></extra>"
            ),
        ),
        secondary_y=False,
    )
    fig.add_trace(
        go.Scatter(
            x=positions,
            y=diversity_score,
            mode="lines",
            name="Diversity score",
            line=dict(color="#ff7f0e", width=1.5),
            hovertemplate=(
                "Position: %{x}<br>"
                "Diversity: %{y}<extra></extra>"
            ),
        ),
        secondary_y=True,
    )

    if segment_conservation_min_pct is not None and segment_diversity_max is not None:
        mid_x, pair_labels = compute_consensus_segment_midpoints(
            positions,
            conservation_pct,
            diversity_score,
            segment_conservation_min_pct,
            segment_diversity_max,
        )
        if mid_x:
            fig.add_trace(
                go.Scatter(
                    x=mid_x,
                    y=[100.0] * len(mid_x),
                    mode="markers",
                    name="Consensus segment midpoints",
                    marker=dict(size=11, color="#2ca02c", line=dict(width=1, color="white")),
                    customdata=pair_labels,
                    hovertemplate=(
                        "x (midpoint): %{x:.2f}<br>"
                        "Conservation: 100% (marker)<br>"
                        "Column pair: %{customdata}<extra></extra>"
                    ),
                ),
                secondary_y=False,
            )

    fig.update_layout(
        title=title,
        hovermode="x unified",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=60, r=60, t=80, b=50),
    )
    fig.update_xaxes(title_text="Alignment position (1-based)")
    fig.update_yaxes(title_text="Conservation %", secondary_y=False, range=[0, 105])
    fig.update_yaxes(title_text="Diversity score (distinct residues)", secondary_y=True)

    return fig
